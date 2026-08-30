function result = solveMOTDGVRPTracked(prob, opts)
%SOLVEMOTDGVRPTRACKED  Same algorithm as solveMOTDGVRP.m (NSGA-2/MLNSGA-2),
%   with an added per-generation hypervolume curve for convergence plots.
%
%   result = solveMOTDGVRPTracked(prob, opts)
%
%   This file is a drop-in copy of solveMOTDGVRP.m with exactly one
%   addition: at the end of every generation, the hypervolume of the
%   CURRENT population's non-dominated front is computed against a
%   caller-supplied reference point and stored in result.hvCurve, so
%   convergence behaviour (not just final-generation quality) can be
%   plotted and compared across algorithms on a shared axis. No other
%   line of the original algorithm is changed, and opts.refPoint is the
%   only new field; every other opts field and the paper's original
%   solveMOTDGVRP.m are unaffected by this file's existence.
%
%   opts fields: identical to solveMOTDGVRP.m, PLUS
%     refPoint   1x3 vector, REQUIRED. Use the same shared reference
%                point across every algorithm being compared on a given
%                instance (e.g. built once from the pooled objectives of
%                all algorithms and all runs, as in
%                run_convergence_analysis.m), so hvCurve values are
%                directly comparable between algorithms -- this mirrors
%                the shared-reference-point convention already used for
%                the paper's final-quality hypervolume comparison
%                (Table 16), applied here per-generation instead of only
%                at the end of the run.
%
%   result fields: identical to solveMOTDGVRP.m, PLUS
%     hvCurve    numGen x 1, normalised hypervolume
%                (raw HV / prod(refPoint)) of the current population's
%                non-dominated front at the end of each generation.

if nargin < 2, opts = struct; end
def = struct('popSize',100,'numGen',100,'pc',0.30,'pm',0.25, ...
             'useML',false,'lambda',2,'maxPatLen',4,'verbose',true,'seed',[], ...
             'refPoint',[]);
opts = mergeOptsTracked(def, opts);
if isempty(opts.refPoint)
    error('solveMOTDGVRPTracked:noRefPoint', ...
        'opts.refPoint is required (a 1x3 vector) so hvCurve is comparable across algorithms.');
end
if ~isempty(opts.seed), rng(opts.seed); end

n   = prob.n;
N   = opts.popSize;
cust = 2:(n+1);

tStart = tic;

pop = cell(N,1);
pop{1} = nearestNeighborChrom(prob);
pop{2} = clarkeWrightChrom(prob);
for i = 3:N
    pop{i} = cust(randperm(n));
end
if opts.useML
    nOpp = round(0.25*N);
    for i = 1:nOpp
        base = pop{randi(2 + min(5,N-2))};
        opp  = oblOpposite(base);
        if aggregateTracked(opp,prob) < aggregateTracked(base,prob)
            pop{i+2} = opp;
        end
    end
end

F = zeros(N,3);
for i = 1:N, F(i,:) = evaluate(pop{i}, prob); end

history = zeros(opts.numGen,3);
hvCurve = zeros(opts.numGen,1);
patterns = {};

for g = 1:opts.numGen
    [~, rank] = fastNonDominatedSort(F);
    cd = zeros(N,1);
    for r = unique(rank)'
        members = find(rank==r);
        cd(members) = crowdingDistance(F(members,:));
    end

    Q = cell(N,1);
    for i = 1:N
        par = tournamentSelect(rank, cd, 2);
        if rand < opts.pc
            child = onePointCrossover(pop{par(1)}, pop{par(2)});
        else
            child = pop{par(1)};
        end
        if rand < opts.pm
            if rand < 0.5
                child = shiftMutation(child);
            else
                child = twoOptMutation(child);
            end
        end
        Q{i} = child;
    end

    if opts.useML && ~isempty(patterns)
        half = round(N/2);
        for i = 1:half
            Q{i} = chromFromPattern(patterns, prob);
        end
    end

    FQ = zeros(N,3);
    for i = 1:N, FQ(i,:) = evaluate(Q{i}, prob); end

    R   = [pop; Q];
    FR  = [F; FQ];
    keep = environmentalSelect(FR, N);
    pop = R(keep);
    F   = FR(keep,:);

    if opts.useML
        [fr,~] = fastNonDominatedSort(F);
        elite = pop(fr{1});
        patterns = mineContiguousPatterns(elite, opts.lambda, opts.maxPatLen);
    end

    history(g,:) = min(F,[],1);

    % ---- NEW: track per-generation hypervolume of the current
    %      non-dominated front against the shared reference point ----
    [frontsHV,~] = fastNonDominatedSort(F);
    currentFront = F(frontsHV{1}, :);
    hvCurve(g) = hypervolume(currentFront, opts.refPoint) / max(prod(opts.refPoint), eps);

    if opts.verbose && (mod(g,10)==0 || g==1)
        fprintf('  gen %3d/%d  bestDist=%.1f  bestTime=%.2f  bestFuel=%.2f  HVnorm=%.4f\n', ...
            g, opts.numGen, history(g,1), history(g,2), history(g,3), hvCurve(g));
    end
end

[fronts,~] = fastNonDominatedSort(F);
first = fronts{1};
[uObj, ia] = unique(round(F(first,:),6), 'rows', 'stable');
result.paretoObj   = uObj;
result.paretoChrom = pop(first(ia));
result.routes = cellfun(@(c) decodeRoutes(c,prob), result.paretoChrom, ...
                        'UniformOutput', false);
result.cpuTime = toc(tStart);
result.history = history;
result.hvCurve = hvCurve;
result.opts    = opts;
end

% =========================================================================
function s = aggregateTracked(chrom, prob)
o = evaluate(chrom, prob);
s = o(1) + o(2)*50 + o(3);
end

function o = mergeOptsTracked(def, in)
o = def;
fn = fieldnames(in);
for i = 1:numel(fn), o.(fn{i}) = in.(fn{i}); end
end
