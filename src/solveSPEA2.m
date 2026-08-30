function result = solveSPEA2(prob, opts)
%SOLVESPEA2  SPEA2 (Zitzler, Laumanns & Thiele, 2001) for the MOTDGVRP.
%
%   result = solveSPEA2(prob, opts)
%
%   Added to address the reviewer's request to compare MOZOA against
%   established multi-objective evolutionary algorithms beyond NSGA-2 and
%   MLNSGA-2. SPEA2 is implemented with the same chromosome representation,
%   objective evaluation, and variation operators as solveMOTDGVRP.m (one-
%   point crossover, shift/2-opt mutation), so the only methodological
%   difference from NSGA-2 in this codebase is SPEA2's environmental
%   selection mechanism: fitness by dominance STRENGTH (Eq. below) combined
%   with a k-th-nearest-neighbour density estimate, and a fixed-size external
%   archive maintained every generation (rather than NSGA-2's rank + crowding
%   distance on a combined parent+offspring population).
%
%   SPEA2 fitness (Zitzler et al. 2001):
%     S(i)     = |{j in P+Arc : i dominates j}|                (strength)
%     R(i)     = sum_{j dominates i} S(j)                       (raw fitness)
%     D(i)     = 1 / (sigma_i^k + 2)                             (density)
%     F(i)     = R(i) + D(i)                                     (lower=better)
%   where sigma_i^k is the distance (in objective space) to the k-th nearest
%   neighbour in P+Arc, k = sqrt(|P|+|Arc|) rounded to nearest integer.
%   Individuals with F(i) < 1 are non-dominated; the archive is filled with
%   these first, truncated by iterative nearest-neighbour removal if there
%   are too many, and padded with the best-F dominated individuals if there
%   are too few (Section 3.2 of Zitzler et al. 2001).
%
%   opts fields (mirrors solveMOTDGVRP.m for a fair, matched comparison):
%     popSize   100      population size
%     archiveSize 100    external archive size (paper default: same as popSize)
%     numGen    100      number of generations
%     pc        0.30     crossover probability
%     pm        0.25     mutation probability
%     verbose   true
%     seed      []
%
%   result fields: same as solveMOTDGVRP (paretoChrom, paretoObj, routes,
%   cpuTime, history).

if nargin < 2, opts = struct; end
def = struct('popSize',100,'archiveSize',100,'numGen',100,'pc',0.30, ...
             'pm',0.25,'verbose',true,'seed',[]);
opts = mergeOptsSPEA2(def, opts);
if ~isempty(opts.seed), rng(opts.seed); end

n    = prob.n;
N    = opts.popSize;
NA   = opts.archiveSize;
cust = 2:(n+1);

tStart = tic;

% ---------- initial population (same hybrid seeding as NSGA-2 driver) ----
pop = cell(N,1);
pop{1} = nearestNeighborChrom(prob);
pop{2} = clarkeWrightChrom(prob);
for i = 3:N
    pop{i} = cust(randperm(n));
end
popObj = zeros(N,3);
for i = 1:N, popObj(i,:) = evaluate(pop{i}, prob); end

archChrom = {};      % external archive chromosomes
archObj   = zeros(0,3);

history = zeros(opts.numGen, 3);

for g = 1:opts.numGen
    % ---- combine population + archive, compute SPEA2 fitness ----
    allChrom = [pop; archChrom(:)];
    allObj   = [popObj; archObj];
    M = size(allObj,1);

    F = spea2Fitness(allObj);

    % ---- environmental selection: build next archive ----
    nd = find(F < 1);                       % non-dominated within allChrom
    if numel(nd) <= NA
        keepIdx = nd;
        if numel(keepIdx) < NA
            % pad with best (lowest-F) dominated individuals
            dominatedIdx = setdiff(1:M, nd);
            [~, ord] = sort(F(dominatedIdx));
            need = NA - numel(keepIdx);
            keepIdx = [keepIdx(:); dominatedIdx(ord(1:min(need,numel(dominatedIdx))))'];
        end
    else
        % truncate by iterative nearest-neighbour removal (Zitzler et al. 2001,
        % truncation operator): repeatedly drop the individual with the
        % smallest distance to its nearest neighbour among the survivors
        keepIdx = nd(:);
        while numel(keepIdx) > NA
            sub = allObj(keepIdx,:);
            k = size(sub,1);
            Dmat = pdist2Local(sub, sub);
            Dmat(1:k+1:end) = inf;          % ignore self-distance
            Dsort = sort(Dmat, 2);
            % remove the individual whose nearest-neighbour distance vector
            % is lexicographically smallest (standard SPEA2 tie-break)
            [~, worst] = min(Dsort(:,1));
            % tie-break by next-nearest distances if needed
            tiedMask = abs(Dsort(:,1) - Dsort(worst,1)) < 1e-12;
            if sum(tiedMask) > 1
                tiedIdx = find(tiedMask);
                [~, sub2] = sortrows(Dsort(tiedIdx,2:end));
                worst = tiedIdx(sub2(1));
            end
            keepIdx(worst) = [];
        end
    end

    archChrom = allChrom(keepIdx);
    archObj   = allObj(keepIdx,:);
    archF     = F(keepIdx);

    history(g,:) = min(archObj,[],1);

    if opts.verbose && (mod(g,10)==0 || g==1)
        fprintf('  SPEA2 gen %4d/%d  |Arc|=%3d  bestD=%.1f bestT=%.2f bestF=%.2f\n', ...
            g, opts.numGen, size(archObj,1), history(g,1), history(g,2), history(g,3));
    end

    if g == opts.numGen, break; end   % skip breeding on the last generation

    % ---- mating selection: binary tournament on SPEA2 fitness (archive-only,
    %      standard SPEA2 convention -- lower F wins) ----
    newPop = cell(N,1);
    for i = 1:N
        p1 = binaryTournament(archChrom, archF);
        p2 = binaryTournament(archChrom, archF);
        if rand < opts.pc
            child = onePointCrossover(p1, p2);
        else
            child = p1;
        end
        if rand < opts.pm
            if rand < 0.5
                child = shiftMutation(child);
            else
                child = twoOptMutation(child);
            end
        end
        newPop{i} = child;
    end
    pop = newPop;
    popObj = zeros(N,3);
    for i = 1:N, popObj(i,:) = evaluate(pop{i}, prob); end
end

% ---------- output: non-dominated subset of the final archive ----------
[fronts, ~] = fastNonDominatedSort(archObj);
firstFront = fronts{1};
[uObj, ia] = unique(round(archObj(firstFront,:),6), 'rows', 'stable');
keptChrom = archChrom(firstFront);

result.paretoChrom = keptChrom(ia);
result.paretoObj = uObj;
result.routes = cellfun(@(c) decodeRoutes(c,prob), result.paretoChrom, ...
                        'UniformOutput', false);
result.cpuTime = toc(tStart);
result.history = history;
result.opts = opts;
end

% =========================================================================
function F = spea2Fitness(objs)
% SPEA2 raw fitness + density estimate (lower is better; <1 = non-dominated).
M = size(objs,1);
domCount = zeros(M,1);    % S(i): how many others i dominates
dominatedBy = cell(M,1);  % indices that dominate i

for i = 1:M
    for j = 1:M
        if i == j, continue; end
        if dominates(objs(i,:), objs(j,:))
            domCount(i) = domCount(i) + 1;
        end
    end
end
for i = 1:M
    for j = 1:M
        if i == j, continue; end
        if dominates(objs(j,:), objs(i,:))
            dominatedBy{i}(end+1) = j; %#ok<AGROW>
        end
    end
end
R = zeros(M,1);
for i = 1:M
    R(i) = sum(domCount(dominatedBy{i}));
end

k = max(1, round(sqrt(M)));
Dmat = pdist2Local(objs, objs);
Dmat(1:M+1:end) = inf;
Dsort = sort(Dmat, 2);
kk = min(k, size(Dsort,2));
sigma_k = Dsort(:,kk);
D = 1 ./ (sigma_k + 2);

F = R + D;
end

function D = pdist2Local(A, B)
% Euclidean pairwise distance without the Statistics Toolbox's pdist2.
% A: p x d, B: q x d  ->  D: p x q
sumA = sum(A.^2, 2);
sumB = sum(B.^2, 2);
D2 = sumA + sumB' - 2*(A*B');
D2(D2 < 0) = 0;             % guard tiny negative values from floating point
D = sqrt(D2);
end

function p = binaryTournament(chroms, F)
n = numel(chroms);
i = randi(n); j = randi(n);
if F(i) <= F(j), p = chroms{i}; else, p = chroms{j}; end
end

function o = mergeOptsSPEA2(def, in)
o = def; fn = fieldnames(in);
for i = 1:numel(fn), o.(fn{i}) = in.(fn{i}); end
end
