function result = solveMOTDGVRP(prob, opts)
%SOLVEMOTDGVRP  NSGA-2 and MLNSGA-2 for the time-dependent green VRP.
%
%   result = solveMOTDGVRP(prob, opts)
%
%   Implements Algorithm 1 (NSGA-2) and Algorithm 2 (MLNSGA-2) of
%   Nyako, Tayachi & Ben Abdelaziz (2025). Set opts.useML = true to enable
%   the machine-learning enhancements: Opposition-Based Learning in the
%   initial population (Section 5.2.1) and Apriori contiguous-pattern
%   injection during the search (Section 5.2.2).
%
%   opts fields (all optional; paper defaults shown):
%     popSize   100      population size
%     numGen    100      number of generations
%     pc        0.30     crossover probability
%     pm        0.25     mutation probability
%     useML     false    NSGA-2 (false) or MLNSGA-2 (true)
%     lambda    2        minimum support for frequent patterns (MLNSGA-2)
%     maxPatLen 4        maximum mined pattern length (MLNSGA-2)
%     verbose   true     print progress
%     seed      []       RNG seed for reproducibility
%
%   result fields:
%     paretoChrom : cell array of non-dominated chromosomes (first front)
%     paretoObj   : K x 3 objective matrix [distance time fuel]
%     routes      : cell array of decoded routes for each Pareto solution
%     cpuTime     : wall-clock seconds
%     history     : per-generation best (min) of each objective

if nargin < 2, opts = struct; end
def = struct('popSize',100,'numGen',100,'pc',0.30,'pm',0.25, ...
             'useML',false,'lambda',2,'maxPatLen',4,'verbose',true,'seed',[]);
opts = mergeOpts(def, opts);
if ~isempty(opts.seed), rng(opts.seed); end

n   = prob.n;
N   = opts.popSize;
cust = 2:(n+1);

tStart = tic;

% ---------- Initial population ----------
pop = cell(N,1);
% Seeds: nearest neighbour + Clarke-Wright, rest random
pop{1} = nearestNeighborChrom(prob);
pop{2} = clarkeWrightChrom(prob);
for i = 3:N
    pop{i} = cust(randperm(n));
end
% MLNSGA-2: inject opposition-based-learning opposites of the seeds
if opts.useML
    nOpp = round(0.25*N);
    for i = 1:nOpp
        base = pop{randi(2 + min(5,N-2))};   % oppose a seed/early individual
        opp  = oblOpposite(base);
        % keep whichever (base vs opposite) is better on aggregated objective
        if aggregate(opp,prob) < aggregate(base,prob)
            pop{i+2} = opp;
        end
    end
end

% ---------- Evaluate ----------
F = zeros(N,3);
for i = 1:N, F(i,:) = evaluate(pop{i}, prob); end

history = zeros(opts.numGen,3);
patterns = {};

% ---------- Evolutionary loop ----------
for g = 1:opts.numGen
    [~, rank] = fastNonDominatedSort(F);
    cd = zeros(N,1);
    for r = unique(rank)'
        members = find(rank==r);
        cd(members) = crowdingDistance(F(members,:));
    end

    % ---- generate offspring Q ----
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

    % MLNSGA-2: replace half of Q with pattern-derived individuals
    if opts.useML && ~isempty(patterns)
        half = round(N/2);
        for i = 1:half
            Q{i} = chromFromPattern(patterns, prob);
        end
    end

    FQ = zeros(N,3);
    for i = 1:N, FQ(i,:) = evaluate(Q{i}, prob); end

    % ---- combine R = P + Q and select N survivors ----
    R   = [pop; Q];
    FR  = [F; FQ];
    keep = environmentalSelect(FR, N);
    pop = R(keep);
    F   = FR(keep,:);

    % MLNSGA-2: update frequent-pattern set from the current first front
    if opts.useML
        [fr,~] = fastNonDominatedSort(F);
        elite = pop(fr{1});
        patterns = mineContiguousPatterns(elite, opts.lambda, opts.maxPatLen);
    end

    history(g,:) = min(F,[],1);
    if opts.verbose && (mod(g,10)==0 || g==1)
        fprintf('  gen %3d/%d  bestDist=%.1f  bestTime=%.2f  bestFuel=%.2f\n', ...
            g, opts.numGen, history(g,1), history(g,2), history(g,3));
    end
end

% ---------- Extract final Pareto front ----------
[fronts,~] = fastNonDominatedSort(F);
first = fronts{1};
% de-duplicate identical objective vectors
[uObj, ia] = unique(round(F(first,:),6), 'rows', 'stable');
result.paretoObj   = uObj;
result.paretoChrom = pop(first(ia));
result.routes = cellfun(@(c) decodeRoutes(c,prob), result.paretoChrom, ...
                        'UniformOutput', false);
result.cpuTime = toc(tStart);
result.history = history;
result.opts    = opts;
end

% =========================================================================
function s = aggregate(chrom, prob)
% Simple scalar surrogate (equal weights on normalised-ish scales) used only
% to decide OBL acceptance; not part of the multi-objective ranking.
o = evaluate(chrom, prob);
s = o(1) + o(2)*50 + o(3);   % rough commensurate scaling
end

function o = mergeOpts(def, in)
o = def;
fn = fieldnames(in);
for i = 1:numel(fn), o.(fn{i}) = in.(fn{i}); end
end
