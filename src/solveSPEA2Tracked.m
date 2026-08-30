function result = solveSPEA2Tracked(prob, opts)
%SOLVESPEA2TRACKED  Same algorithm as solveSPEA2.m, with an added
%   per-generation hypervolume curve for convergence plots.
%
%   result = solveSPEA2Tracked(prob, opts)
%
%   Drop-in copy of solveSPEA2.m with one addition: at the end of every
%   generation, the hypervolume of the current external archive (which is
%   always non-dominated by SPEA2's own construction) is computed against
%   a caller-supplied shared reference point and stored in
%   result.hvCurve. No other line of the original algorithm is changed.
%
%   opts fields: identical to solveSPEA2.m, PLUS
%     refPoint   1x3 vector, REQUIRED (see solveMOTDGVRPTracked.m for why
%                a shared reference point across algorithms matters here).
%
%   result fields: identical to solveSPEA2.m, PLUS
%     hvCurve    numGen x 1, normalised hypervolume of the archive at the
%                end of each generation.

if nargin < 2, opts = struct; end
def = struct('popSize',100,'archiveSize',100,'numGen',100,'pc',0.30, ...
             'pm',0.25,'verbose',true,'seed',[],'refPoint',[]);
opts = mergeOptsSPEA2Tracked(def, opts);
if isempty(opts.refPoint)
    error('solveSPEA2Tracked:noRefPoint', ...
        'opts.refPoint is required (a 1x3 vector) so hvCurve is comparable across algorithms.');
end
if ~isempty(opts.seed), rng(opts.seed); end

n    = prob.n;
N    = opts.popSize;
NA   = opts.archiveSize;
cust = 2:(n+1);

tStart = tic;

pop = cell(N,1);
pop{1} = nearestNeighborChrom(prob);
pop{2} = clarkeWrightChrom(prob);
for i = 3:N
    pop{i} = cust(randperm(n));
end
popObj = zeros(N,3);
for i = 1:N, popObj(i,:) = evaluate(pop{i}, prob); end

archChrom = {};
archObj   = zeros(0,3);

history = zeros(opts.numGen, 3);
hvCurve = zeros(opts.numGen, 1);

for g = 1:opts.numGen
    allChrom = [pop; archChrom(:)];
    allObj   = [popObj; archObj];
    M = size(allObj,1);

    F = spea2FitnessTracked(allObj);

    nd = find(F < 1);
    if numel(nd) <= NA
        keepIdx = nd;
        if numel(keepIdx) < NA
            dominatedIdx = setdiff(1:M, nd);
            [~, ord] = sort(F(dominatedIdx));
            need = NA - numel(keepIdx);
            keepIdx = [keepIdx(:); dominatedIdx(ord(1:min(need,numel(dominatedIdx))))'];
        end
    else
        keepIdx = nd(:);
        while numel(keepIdx) > NA
            sub = allObj(keepIdx,:);
            k = size(sub,1);
            Dmat = pdist2LocalTracked(sub, sub);
            Dmat(1:k+1:end) = inf;
            Dsort = sort(Dmat, 2);
            [~, worst] = min(Dsort(:,1));
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

    % ---- NEW: track per-generation hypervolume of the (already
    %      non-dominated, by SPEA2's own construction) archive ----
    [frontsHV,~] = fastNonDominatedSort(archObj);
    currentFront = archObj(frontsHV{1}, :);
    hvCurve(g) = hypervolume(currentFront, opts.refPoint) / max(prod(opts.refPoint), eps);

    if opts.verbose && (mod(g,10)==0 || g==1)
        fprintf('  SPEA2 gen %4d/%d  |Arc|=%3d  bestD=%.1f bestT=%.2f bestF=%.2f  HVnorm=%.4f\n', ...
            g, opts.numGen, size(archObj,1), history(g,1), history(g,2), history(g,3), hvCurve(g));
    end

    if g == opts.numGen, break; end

    newPop = cell(N,1);
    for i = 1:N
        p1 = binaryTournamentTracked(archChrom, archF);
        p2 = binaryTournamentTracked(archChrom, archF);
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
result.hvCurve = hvCurve;
result.opts = opts;
end

% =========================================================================
function F = spea2FitnessTracked(objs)
M = size(objs,1);
domCount = zeros(M,1);
dominatedBy = cell(M,1);

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
Dmat = pdist2LocalTracked(objs, objs);
Dmat(1:M+1:end) = inf;
Dsort = sort(Dmat, 2);
kk = min(k, size(Dsort,2));
sigma_k = Dsort(:,kk);
D = 1 ./ (sigma_k + 2);

F = R + D;
end

function D = pdist2LocalTracked(A, B)
sumA = sum(A.^2, 2);
sumB = sum(B.^2, 2);
D2 = sumA + sumB' - 2*(A*B');
D2(D2 < 0) = 0;
D = sqrt(D2);
end

function p = binaryTournamentTracked(chroms, F)
n = numel(chroms);
i = randi(n); j = randi(n);
if F(i) <= F(j), p = chroms{i}; else, p = chroms{j}; end
end

function o = mergeOptsSPEA2Tracked(def, in)
o = def; fn = fieldnames(in);
for i = 1:numel(fn), o.(fn{i}) = in.(fn{i}); end
end
