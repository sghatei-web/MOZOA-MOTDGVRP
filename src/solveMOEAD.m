function result = solveMOEAD(prob, opts)
%SOLVEMOEAD  MOEA/D (Zhang & Li, 2007) for the MOTDGVRP, Tchebycheff variant.
%
%   result = solveMOEAD(prob, opts)
%
%   Added alongside solveSPEA2.m to address the reviewer's request for
%   decomposition-based competitors. MOEA/D decomposes the 3-objective
%   problem into opts.popSize scalar Tchebycheff subproblems, each with its
%   own weight vector lambda_i (generated on the simplex via a systematic
%   grid, standard for 3 objectives), and each subproblem is optimised
%   cooperatively using only its T nearest neighbouring subproblems (by
%   weight-vector distance) for mating, following Algorithm 1 of Zhang & Li
%   (2007). Uses the same chromosome representation and variation operators
%   (onePointCrossover, shiftMutation, twoOptMutation) as the other
%   competitors in this codebase, so the comparison isolates the search
%   strategy rather than the encoding.
%
%   Tchebycheff scalarising function for subproblem i:
%     g(x | lambda_i, z*) = max_k { lambda_i(k) * |f_k(x) - z*_k| }
%   where z* is the current ideal point (component-wise best objective seen).
%
%   opts fields:
%     popSize   100      number of subproblems (= population size)
%     numGen    100      number of generations (matched to the other methods)
%     T         20       neighbourhood size (number of closest weight vectors)
%     pc        0.30     crossover probability
%     pm        0.25     mutation probability
%     delta     0.9      probability of selecting mates from the neighbourhood
%                         (vs. the whole population)
%     nr        2        max number of neighbour solutions replaced per child
%     verbose   true
%     seed      []
%
%   result fields: same as solveMOTDGVRP (paretoChrom, paretoObj, routes,
%   cpuTime, history). The returned Pareto set is the non-dominated subset of
%   the final population (the standard way to report a MOEA/D result).

if nargin < 2, opts = struct; end
def = struct('popSize',100,'numGen',100,'T',20,'pc',0.30,'pm',0.25, ...
             'delta',0.9,'nr',2,'verbose',true,'seed',[]);
opts = mergeOptsMOEAD(def, opts);
if ~isempty(opts.seed), rng(opts.seed); end

n    = prob.n;
N    = opts.popSize;
cust = 2:(n+1);

tStart = tic;

% ---------- weight vectors: simplex-lattice for 3 objectives, resized to N --
W = simplexWeights3(N);
N = size(W,1);   % simplex-lattice sizes only hit N approximately; use actual

% ---------- neighbourhoods: T closest weight vectors (Euclidean) ----------
Dw = pdist2Local(W, W);
[~, nbOrd] = sort(Dw, 2);
T = min(opts.T, N);
B = nbOrd(:, 1:T);     % N x T neighbour index list

% ---------- initial population ----------
pop = cell(N,1);
pop{1} = nearestNeighborChrom(prob);
if N > 1, pop{2} = clarkeWrightChrom(prob); end
for i = 3:N
    pop{i} = cust(randperm(n));
end
popObj = zeros(N,3);
for i = 1:N, popObj(i,:) = evaluate(pop{i}, prob); end

zStar = min(popObj,[],1);   % ideal point

history = zeros(opts.numGen, 3);

for g = 1:opts.numGen
    for i = 1:N
        % ---- mating pool: neighbourhood with prob delta, else whole pop ----
        if rand < opts.delta
            pool = B(i,:);
        else
            pool = 1:N;
        end
        idx = pool(randperm(numel(pool), min(2,numel(pool))));
        p1 = pop{idx(1)};
        if numel(idx) >= 2
            p2 = pop{idx(2)};
        else
            p2 = pop{idx(1)};
        end

        if rand < opts.pc
            child = onePointCrossover(p1, p2);
        else
            child = p1;
        end
        if rand < opts.pm
            if rand < 0.5, child = shiftMutation(child);
            else,          child = twoOptMutation(child); end
        end

        childObj = evaluate(child, prob);
        zStar = min(zStar, childObj);

        % ---- update neighbouring subproblems (bounded replacement, nr) ----
        replaceCount = 0;
        order = pool(randperm(numel(pool)));
        for jj = order
            gNew = tchebycheff(childObj, W(jj,:), zStar);
            gOld = tchebycheff(popObj(jj,:), W(jj,:), zStar);
            if gNew <= gOld
                pop{jj} = child;
                popObj(jj,:) = childObj;
                replaceCount = replaceCount + 1;
                if replaceCount >= opts.nr, break; end
            end
        end
    end

    history(g,:) = min(popObj,[],1);
    if opts.verbose && (mod(g,10)==0 || g==1)
        fprintf('  MOEA/D gen %4d/%d  bestD=%.1f bestT=%.2f bestF=%.2f\n', ...
            g, opts.numGen, history(g,1), history(g,2), history(g,3));
    end
end

% ---------- output: non-dominated subset of the final population ----------
[fronts, ~] = fastNonDominatedSort(popObj);
firstFront = fronts{1};
[uObj, ia] = unique(round(popObj(firstFront,:),6), 'rows', 'stable');
keptChrom = pop(firstFront);

result.paretoChrom = keptChrom(ia);
result.paretoObj = uObj;
result.routes = cellfun(@(c) decodeRoutes(c,prob), result.paretoChrom, ...
                        'UniformOutput', false);
result.cpuTime = toc(tStart);
result.history = history;
result.opts = opts;
end

% =========================================================================
function g = tchebycheff(obj, lambda, zStar)
lambda = max(lambda, 1e-6);   % avoid zero weight collapsing the scalarisation
g = max(lambda .* abs(obj - zStar));
end

function W = simplexWeights3(Napprox)
% Systematic simplex-lattice weight vectors for 3 objectives (Das & Dennis
% 1998 construction), sized to produce approximately Napprox vectors.
% Number of points for divisions H is C(H+2,2) = (H+1)(H+2)/2.
H = 1;
while (H+1)*(H+2)/2 < Napprox, H = H + 1; end
W = [];
for i = 0:H
    for j = 0:(H-i)
        k = H - i - j;
        W(end+1,:) = [i, j, k] / H; %#ok<AGROW>
    end
end
end

function D = pdist2Local(A, B)
sumA = sum(A.^2, 2);
sumB = sum(B.^2, 2);
D2 = sumA + sumB' - 2*(A*B');
D2(D2 < 0) = 0;
D = sqrt(D2);
end

function o = mergeOptsMOEAD(def, in)
o = def; fn = fieldnames(in);
for i = 1:numel(fn), o.(fn{i}) = in.(fn{i}); end
end
