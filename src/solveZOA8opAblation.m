function result = solveZOA8opAblation(prob, opts)
%SOLVEZOA8OPABLATION  MOZOA with individual components switched off, for
%   the ablation study requested during review (run_ablation_study.m).
%
%   result = solveZOA8opAblation(prob, opts)
%
%   This is solveZOA8op.m with three additional opts fields that each
%   disable one structural component of MOZOA, so the contribution of
%   each can be isolated:
%
%     opts.variant  'full'            -- MOZOA exactly as in the paper
%                                        (identical to solveZOA8op.m)
%                   'noArchive'       -- the external Pareto archive is
%                                        replaced by a single running-best
%                                        "elite" solution per objective
%                                        combination is NOT tracked;
%                                        instead only the CURRENT
%                                        population is used for
%                                        leader/partner selection and for
%                                        computing hvCurve, i.e. there is
%                                        no persistent archive across
%                                        iterations beyond the population
%                                        itself. Acceptance falls back to
%                                        a simple per-individual Pareto
%                                        dominance check against that
%                                        individual's own previous
%                                        objective (no archive-entry
%                                        signal is available).
%                   'noRecombination' -- operators O6 (archive OX) and O7
%                                        (leader OX) are removed from the
%                                        operator set; the remaining six
%                                        operators (O1,O2,O3,O4,O5,O8)
%                                        are drawn with their relative
%                                        weights from opts.fixedProbs
%                                        renormalised to sum to 1 (or
%                                        uniformly among the six if
%                                        opts.opSelectMode='random')
%                   'uniformSelection' -- identical structure to 'full'
%                                        (archive + all 8 operators) but
%                                        opSelectMode is forced to
%                                        'random' regardless of what is
%                                        passed in, i.e. this is exactly
%                                        solveZOA8op.m with
%                                        opSelectMode='random' -- kept
%                                        here as a named ablation
%                                        condition for convenience so all
%                                        four conditions can be driven
%                                        from one function with one
%                                        switch.
%
%   All other opts fields (popSize, numIter, Amax, p1, fixedProbs,
%   verbose, seed) match solveZOA8op.m exactly, and 'full' with default
%   opts reproduces solveZOA8op.m's behaviour bit-for-bit (same RNG
%   sequence given the same seed), so this file can be used as a
%   drop-in replacement for solveZOA8op.m in run_ablation_study.m without
%   changing any other paper result.
%
%   result fields: same as solveZOA8op.m, plus result.variant echoing
%   the condition that was run, and (for 'full'/'uniformSelection' only,
%   since these are the only conditions with a real persistent archive)
%   result.archiveSize, matching solveZOA8op.m.

if nargin < 2, opts = struct; end
def = struct('popSize',100,'numIter',600,'Amax',100,'p1',0.60, ...
             'opSelectMode','fixed', ...
             'fixedProbs',[0.10 0.10 0.20 0.05 0.10 0.20 0.20 0.05], ...
             'variant','full', ...
             'verbose',true,'seed',[]);
opts = mergeOptsAblation(def, opts);
if ~isempty(opts.seed), rng(opts.seed); end

validVariants = {'full','noArchive','noRecombination','uniformSelection'};
if ~ismember(opts.variant, validVariants)
    error('solveZOA8opAblation:badVariant', ...
        'opts.variant must be one of: %s', strjoin(validVariants, ', '));
end

% 'uniformSelection' is structurally identical to 'full' except the
% operator-selection mode; force it here so the rest of the function
% only needs to branch on the archive/recombination conditions.
if strcmp(opts.variant, 'uniformSelection')
    opts.opSelectMode = 'random';
end

% ---- build the active operator index list and matching probabilities ----
if strcmp(opts.variant, 'noRecombination')
    activeOps = [1 2 3 4 5 8];                 % O6, O7 removed
    baseProbs = opts.fixedProbs([1 2 3 4 5 8]);
    baseProbs = baseProbs / sum(baseProbs);     % renormalise to sum 1
else
    activeOps = 1:8;
    baseProbs = opts.fixedProbs;
    if abs(sum(baseProbs) - 1) > 1e-9
        error('solveZOA8opAblation:fixedProbs', ...
            'opts.fixedProbs must sum to 1 (got %.6f).', sum(baseProbs));
    end
end
nAct = numel(activeOps);
useArchive = ~strcmp(opts.variant, 'noArchive');

n    = prob.n;
NP   = opts.popSize;
cust = 2:(n+1);

tStart = tic;

% ---------- hybrid initialisation (identical to solveZOA8op.m) ----------
pop = cell(NP,1);
popObj = zeros(NP,3);
half = round(NP/2);
pop{1} = nearestNeighborChrom(prob);
pop{2} = clarkeWrightChrom(prob);
for i = 3:half
    base = nearestNeighborChrom(prob);
    pop{i} = segmentRelocation(base);
end
for i = half+1:NP
    pop{i} = cust(randperm(n));
end
for i = 1:NP, popObj(i,:) = evaluate(pop{i}, prob); end

refPoint = max(popObj,[],1) * 1.1 + 1;

if useArchive
    A = struct('obj',{},'chrom',{});
    for i = 1:NP
        cand.obj = popObj(i,:); cand.chrom = pop{i};
        A = archiveUpdate(A, cand, opts.Amax);
    end
else
    % No persistent archive: 'A' here is recomputed EVERY iteration from
    % only the current population's non-dominated subset, so nothing
    % survives across iterations beyond what is currently in the
    % population -- this is the deliberate ablation of the archive
    % component. leader/partner selection below draws from this
    % transient set instead of a true external archive.
    A = paretoSubsetAblation(pop, popObj);
end

hvCurve  = zeros(opts.numIter,1);
archSize = zeros(opts.numIter,1);
% Track operator usage statistics throughout the run (used by
% run_operator_usage_analysis.m; harmless overhead for the plain
% ablation study, which ignores these fields).
opTried    = zeros(1,8);
opAccepted = zeros(1,8);
opArchived = zeros(1,8);

for t = 1:opts.numIter
    if useArchive
        archiveObjs = vertcat(A.obj);
    else
        A = paretoSubsetAblation(pop, popObj);
        archiveObjs = vertcat(A.obj);
    end

    for i = 1:NP
        phase2 = (rand < opts.p1);

        % ----- operator selection (restricted to activeOps) -----
        if strcmp(opts.opSelectMode,'random')
            a = activeOps(randi(nAct));
        else
            pick = find(rand <= cumsum(baseProbs), 1, 'first');
            if isempty(pick), pick = nAct; end
            a = activeOps(pick);
        end
        opTried(a) = opTried(a) + 1;

        % ----- apply operator -----
        opCtx = struct();
        opCtx.dist = prob.dist;
        opCtx.Wi   = 2 - 2*(t/opts.numIter);
        if ~isempty(A)
            AO = archiveObjs;
            rngO = max(AO,[],1) - min(AO,[],1); rngO(rngO==0) = 1;
            score = sum((AO - min(AO,[],1))./rngO, 2);
            [~, li] = min(score);
            opCtx.leader  = A(li).chrom;
            opCtx.partner = A(randi(numel(A))).chrom;
        end
        newChrom = applyOperator(pop{i}, a, opCtx);
        newObj = evaluate(newChrom, prob);

        cand.obj = newObj; cand.chrom = newChrom;
        if useArchive
            [A, entered] = archiveUpdate(A, cand, opts.Amax);
            archiveObjs = vertcat(A.obj);
        else
            % No persistent archive: "entered" is redefined as "this
            % candidate is non-dominated within the current population
            % plus itself", which is the closest available analogue of
            % archive entry when there is no archive to enter.
            entered = ~any(all(popObj <= newObj, 2) & any(popObj < newObj, 2));
        end
        if entered, opArchived(a) = opArchived(a) + 1; end

        accept = false;
        if phase2
            accept = entered;
        else
            if entered || softParetoImproveAblation(newObj, popObj(i,:))
                accept = true;
            end
        end
        if accept
            pop{i} = newChrom;
            popObj(i,:) = newObj;
            opAccepted(a) = opAccepted(a) + 1;
        end
    end

    if useArchive
        HVnorm = hypervolume(archiveObjs, refPoint) / max(prod(refPoint), eps);
        archSize(t) = numel(A);
    else
        currentFront = paretoSubsetAblation(pop, popObj);
        HVnorm = hypervolume(vertcat(currentFront.obj), refPoint) / max(prod(refPoint), eps);
        archSize(t) = numel(currentFront);
    end
    hvCurve(t) = HVnorm;

    if opts.verbose && (mod(t,50)==0 || t==1)
        fprintf('  Ablation[%s] iter %4d/%d  |A|=%3d  HVnorm=%.4f\n', ...
            opts.variant, t, opts.numIter, archSize(t), HVnorm);
    end
end

% ---------- output ----------
if useArchive
    archiveObjs = vertcat(A.obj);
    [uObj, ia] = unique(round(archiveObjs,6), 'rows', 'stable');
    result.paretoChrom = {A(ia).chrom};
else
    finalFront = paretoSubsetAblation(pop, popObj);
    archiveObjs = vertcat(finalFront.obj);
    [uObj, ia] = unique(round(archiveObjs,6), 'rows', 'stable');
    result.paretoChrom = {finalFront(ia).chrom};
end
result.paretoObj = uObj;
result.routes = cellfun(@(c) decodeRoutes(c,prob), result.paretoChrom, ...
                        'UniformOutput', false);
result.archiveSize = archSize;
result.hvCurve = hvCurve;
result.refPoint = refPoint;
result.cpuTime = toc(tStart);
result.variant = opts.variant;
result.opTried = opTried;
result.opAccepted = opAccepted;
result.opArchived = opArchived;
result.opts = opts;
end

% =========================================================================
function A = paretoSubsetAblation(pop, popObj)
% Non-dominated subset of the CURRENT population only (no memory of past
% iterations) -- used to emulate "no archive" while still allowing
% leader/partner selection and an HV curve to be computed each iteration.
n = size(popObj,1);
keep = true(n,1);
for i = 1:n
    if ~keep(i), continue; end
    for j = 1:n
        if i~=j && keep(j)
            if all(popObj(j,:) <= popObj(i,:)) && any(popObj(j,:) < popObj(i,:))
                keep(i) = false; break;
            end
        end
    end
end
idx = find(keep);
A = struct('obj', num2cell(popObj(idx,:),2)', 'chrom', pop(idx)');
end

function tf = softParetoImproveAblation(newObj, curObj)
tf = all(newObj <= curObj) && any(newObj < curObj);
end

function o = mergeOptsAblation(def, in)
o = def; fn = fieldnames(in);
for i = 1:numel(fn), o.(fn{i}) = in.(fn{i}); end
end
