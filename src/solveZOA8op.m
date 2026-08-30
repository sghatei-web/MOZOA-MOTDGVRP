function result = solveZOA8op(prob, opts)
%SOLVMOZOA8OP  Multi-Operator Zebra Optimization Algorithm (MOZOA).
%
%   result = solveZOA8op(prob, opts)
%
%   The method combines eight discrete routing operators (five local moves,
%   two Order-Crossover recombination operators, and a double-bridge
%   perturbation), two-phase foraging/defence acceptance, and an external
%   Pareto archive with crowding-distance truncation. Operators are selected
%   using either the fixed schedule reported in the paper or a uniform-random
%   schedule for the controlled schedule comparison.
%
%   opts fields:
%     popSize      100     NP zebras
%     numIter      600     T_max iterations (paper experiments explicitly
%                          set this value to 100 through their run scripts)
%     Amax         100     max Pareto-archive size
%     p1           0.60    probability of entering Phase-2 (defense/greedy
%                          acceptance), same zebra behavioural split as
%                          MOZOA
%     opSelectMode 'fixed' | 'random'
%                          'fixed'  (default): operators drawn with the
%                                   fixed probabilities opts.fixedProbs
%                                   (default: weighted toward the
%                                   intensification/recombination operators
%                                   O3/O6/O7, matching the non-uniform
%                                   default used in run_ablation_study.m)
%                          'random': operators drawn uniformly at random
%                                   every step
%     fixedProbs   1x8 vector, used only when opSelectMode='fixed', must sum
%                  to 1 (default: [0.10 0.10 0.20 0.05 0.10 0.20 0.20 0.05])
%     verbose      true
%     seed         []
%
%   result fields: same as solveMOZOA.m (paretoObj, paretoChrom, routes,
%   archiveSize, hvCurve, refPoint, cpuTime).

if nargin < 2, opts = struct; end
def = struct('popSize',100,'numIter',600,'Amax',100,'p1',0.60, ...
             'opSelectMode','fixed', ...
             'fixedProbs',[0.10 0.10 0.20 0.05 0.10 0.20 0.20 0.05], ...
             'verbose',true,'seed',[]);
opts = mergeOptsZOA8(def, opts);
if ~isempty(opts.seed), rng(opts.seed); end
if strcmp(opts.opSelectMode,'fixed') && abs(sum(opts.fixedProbs) - 1) > 1e-9
    error('solveZOA8op:fixedProbs', ...
        'opts.fixedProbs must sum to 1 (got %.6f).', sum(opts.fixedProbs));
end

n    = prob.n;
NP   = opts.popSize;
cust = 2:(n+1);
nAct = 8;

tStart = tic;

% ---------- hybrid initialisation (identical to solveMOZOA.m) -----------
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

% ---------- Pareto archive: reference point ----------
refPoint = max(popObj,[],1) * 1.1 + 1;

A = struct('obj',{},'chrom',{});
for i = 1:NP
    cand.obj = popObj(i,:); cand.chrom = pop{i};
    A = archiveUpdate(A, cand, opts.Amax);
end

hvCurve  = zeros(opts.numIter,1);
archSize = zeros(opts.numIter,1);

for t = 1:opts.numIter
    archiveObjs = vertcat(A.obj);

    for i = 1:NP
        phase2 = (rand < opts.p1);

        % ----- operator selection: no learning, no network -----
        if strcmp(opts.opSelectMode,'random')
            a = randi(nAct);
        else
            a = find(rand <= cumsum(opts.fixedProbs), 1, 'first');
            if isempty(a), a = nAct; end
        end

        % ----- apply operator (same MN-ZOA context as MOZOA) -----
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
        [A, entered] = archiveUpdate(A, cand, opts.Amax);
        archiveObjs = vertcat(A.obj);

        % ----- same two-phase zebra acceptance rule as MOZOA -----
        accept = false;
        if phase2
            accept = entered;
        else
            if entered || softParetoImproveZOA8(newObj, popObj(i,:))
                accept = true;
            end
        end
        if accept
            pop{i} = newChrom;
            popObj(i,:) = newObj;
        end
    end

    HVnorm = hypervolume(archiveObjs, refPoint) / max(prod(refPoint), eps);
    hvCurve(t)  = HVnorm;
    archSize(t) = numel(A);

    if opts.verbose && (mod(t,50)==0 || t==1)
        bo = min(archiveObjs,[],1);
        fprintf('  ZOA-8op iter %4d/%d  |A|=%3d  HVnorm=%.4f  bestD=%.1f bestT=%.2f bestF=%.2f\n', ...
            t, opts.numIter, numel(A), HVnorm, bo(1), bo(2), bo(3));
    end
end

% ---------- output ----------
archiveObjs = vertcat(A.obj);
[uObj, ia] = unique(round(archiveObjs,6), 'rows', 'stable');
result.paretoObj   = uObj;
result.paretoChrom = {A(ia).chrom};
result.routes = cellfun(@(c) decodeRoutes(c,prob), result.paretoChrom, ...
                        'UniformOutput', false);
result.archiveSize = archSize;
result.hvCurve = hvCurve;
result.refPoint = refPoint;
result.cpuTime = toc(tStart);
result.opts = opts;
end

% =========================================================================
function tf = softParetoImproveZOA8(newObj, curObj)
tf = all(newObj <= curObj) && any(newObj < curObj);
end

function o = mergeOptsZOA8(def, in)
o = def; fn = fieldnames(in);
for i = 1:numel(fn), o.(fn{i}) = in.(fn{i}); end
end
