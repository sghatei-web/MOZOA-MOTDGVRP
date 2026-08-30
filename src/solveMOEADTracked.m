function result = solveMOEADTracked(prob, opts)
%SOLVEMOEADTRACKED  Same algorithm as solveMOEAD.m, with an added
%   per-generation hypervolume curve for convergence plots.
%
%   result = solveMOEADTracked(prob, opts)
%
%   Drop-in copy of solveMOEAD.m with one addition: at the end of every
%   generation, the hypervolume of the current population's non-dominated
%   subset is computed against a caller-supplied shared reference point
%   and stored in result.hvCurve. No other line of the original algorithm
%   is changed.
%
%   opts fields: identical to solveMOEAD.m, PLUS
%     refPoint   1x3 vector, REQUIRED (see solveMOTDGVRPTracked.m for why
%                a shared reference point across algorithms matters here).
%
%   result fields: identical to solveMOEAD.m, PLUS
%     hvCurve    numGen x 1, normalised hypervolume of the population's
%                non-dominated subset at the end of each generation.

if nargin < 2, opts = struct; end
def = struct('popSize',100,'numGen',100,'T',20,'pc',0.30,'pm',0.25, ...
             'delta',0.9,'nr',2,'verbose',true,'seed',[],'refPoint',[]);
opts = mergeOptsMOEADTracked(def, opts);
if isempty(opts.refPoint)
    error('solveMOEADTracked:noRefPoint', ...
        'opts.refPoint is required (a 1x3 vector) so hvCurve is comparable across algorithms.');
end
if ~isempty(opts.seed), rng(opts.seed); end

n    = prob.n;
N    = opts.popSize;
cust = 2:(n+1);

tStart = tic;

W = simplexWeights3Tracked(N);
N = size(W,1);

Dw = pdist2LocalMOEADTracked(W, W);
[~, nbOrd] = sort(Dw, 2);
T = min(opts.T, N);
B = nbOrd(:, 1:T);

pop = cell(N,1);
pop{1} = nearestNeighborChrom(prob);
if N > 1, pop{2} = clarkeWrightChrom(prob); end
for i = 3:N
    pop{i} = cust(randperm(n));
end
popObj = zeros(N,3);
for i = 1:N, popObj(i,:) = evaluate(pop{i}, prob); end

zStar = min(popObj,[],1);

history = zeros(opts.numGen, 3);
hvCurve = zeros(opts.numGen, 1);

for g = 1:opts.numGen
    for i = 1:N
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

        replaceCount = 0;
        order = pool(randperm(numel(pool)));
        for jj = order
            gNew = tchebycheffTracked(childObj, W(jj,:), zStar);
            gOld = tchebycheffTracked(popObj(jj,:), W(jj,:), zStar);
            if gNew <= gOld
                pop{jj} = child;
                popObj(jj,:) = childObj;
                replaceCount = replaceCount + 1;
                if replaceCount >= opts.nr, break; end
            end
        end
    end

    history(g,:) = min(popObj,[],1);

    % ---- NEW: track per-generation hypervolume of the population's
    %      non-dominated subset ----
    [frontsHV,~] = fastNonDominatedSort(popObj);
    currentFront = popObj(frontsHV{1}, :);
    hvCurve(g) = hypervolume(currentFront, opts.refPoint) / max(prod(opts.refPoint), eps);

    if opts.verbose && (mod(g,10)==0 || g==1)
        fprintf('  MOEA/D gen %4d/%d  bestD=%.1f bestT=%.2f bestF=%.2f  HVnorm=%.4f\n', ...
            g, opts.numGen, history(g,1), history(g,2), history(g,3), hvCurve(g));
    end
end

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
result.hvCurve = hvCurve;
result.opts = opts;
end

% =========================================================================
function g = tchebycheffTracked(obj, lambda, zStar)
lambda = max(lambda, 1e-6);
g = max(lambda .* abs(obj - zStar));
end

function W = simplexWeights3Tracked(Napprox)
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

function D = pdist2LocalMOEADTracked(A, B)
sumA = sum(A.^2, 2);
sumB = sum(B.^2, 2);
D2 = sumA + sumB' - 2*(A*B');
D2(D2 < 0) = 0;
D = sqrt(D2);
end

function o = mergeOptsMOEADTracked(def, in)
o = def; fn = fieldnames(in);
for i = 1:numel(fn), o.(fn{i}) = in.(fn{i}); end
end
