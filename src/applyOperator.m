function c = applyOperator(chrom, a, ctx)
%APPLYOPERATOR  Apply one of MOZOA's discrete operators (MN-ZOA enhanced).
%
%   c = applyOperator(chrom, a, ctx)
%
%   a = 1 : Swap                        (O1, mutation)
%   a = 2 : Insertion                   (O2, mutation)
%   a = 3 : 2-opt full local search     (O3, intensification)  *MN-ZOA*
%   a = 4 : 3-opt move                  (O4, diversification)
%   a = 5 : Segment relocation          (O5, mutation)
%   a = 6 : Order crossover w/ archive  (O6, recombination)    *MN-ZOA*
%   a = 7 : Order crossover w/ leader   (O7, recombination)    *MN-ZOA*
%   a = 8 : Double-bridge restart       (O8, perturbation)     *MN-ZOA*
%
%   The last three operators are adapted from the MN-ZOA / ADZOA design,
%   which strengthened the discrete Zebra algorithm by replacing destructive
%   positional moves with Order Crossover (OX) that preserves good
%   sub-sequences, a full 2-opt local search whose intensity follows an
%   adaptive weight, and a double-bridge kick for stalled individuals.
%
%   ctx (optional) supplies the recombination/intensification context:
%     .dist     N x N distance matrix (enables full 2-opt local search)
%     .partner  a chromosome to recombine with (an archive member; O6)
%     .leader   the current best/leader chromosome (O7)
%     .Wi       adaptive weight in [0,2] controlling 2-opt passes
%   When ctx or a needed field is missing, the operator degrades gracefully
%   to a mutation, so callers that pass only (chrom, a) still work.

if nargin < 3, ctx = struct; end

switch a
    case 1
        c = swapOp(chrom);
    case 2
        c = insertionOp(chrom);
    case 3
        c = twoOptOp(chrom, ctx);          % full local search if dist given
    case 4
        c = threeOptOp(chrom);
    case 5
        c = segmentRelocation(chrom);
    case 6
        c = crossoverOp(chrom, ctx, 'partner');
    case 7
        c = crossoverOp(chrom, ctx, 'leader');
    case 8
        c = doubleBridgeOp(chrom);
    otherwise
        c = chrom;
end
end

% ---- O1 Swap -------------------------------------------------------------
function c = swapOp(t)
n = numel(t); c = t;
if n < 2, return; end
idx = randperm(n,2);
c([idx(1) idx(2)]) = t([idx(2) idx(1)]);
end

% ---- O2 Insertion --------------------------------------------------------
function c = insertionOp(t)
n = numel(t);
if n < 2, c = t; return; end
idx = randperm(n,2); i = idx(1); j = idx(2);
city = t(i); c = t; c(i) = [];
if j > i, j = j - 1; end
c = [c(1:j), city, c(j+1:end)];
end

% ---- O3 2-opt: full local search when a distance matrix is available -----
function c = twoOptOp(t, ctx)
if isfield(ctx,'dist') && ~isempty(ctx.dist)
    if isfield(ctx,'Wi') && ~isempty(ctx.Wi)
        maxPasses = max(1, round(ctx.Wi * 5) + 1);   % adaptive (MN-ZOA)
    else
        maxPasses = 6;
    end
    c = twoOptLocalSearchTour(t, ctx.dist, maxPasses);
else
    n = numel(t);
    if n < 2, c = t; return; end
    i = randi(n); j = randi(n);
    if i > j, [i,j] = deal(j,i); end
    c = t; c(i:j) = t(j:-1:i);
end
end

% ---- O4 3-opt (two reversals for stronger diversification) --------------
function c = threeOptOp(t)
n = numel(t);
if n < 4, c = twoOptOp(t, struct()); return; end
cuts = sort(randperm(n-1, 3));
i = cuts(1); j = cuts(2); k = cuts(3);
A = t(1:i); B = t(i+1:j); C = t(j+1:k); D = t(k+1:end);
c = [A, fliplr(B), fliplr(C), D];
end

% ---- O6/O7 Order Crossover (OX) with an archive partner or the leader ----
function c = crossoverOp(t, ctx, which)
switch which
    case 'partner'
        if isfield(ctx,'partner') && ~isempty(ctx.partner)
            mate = ctx.partner;
        elseif isfield(ctx,'leader') && ~isempty(ctx.leader)
            mate = ctx.leader;
        else
            c = swapOp(t); return;
        end
    otherwise
        if isfield(ctx,'leader') && ~isempty(ctx.leader)
            mate = ctx.leader;
        elseif isfield(ctx,'partner') && ~isempty(ctx.partner)
            mate = ctx.partner;
        else
            c = swapOp(t); return;
        end
end
c = orderCrossoverTour(t, mate);
end

% ---- O8 Double-bridge perturbation (4-opt double bridge) ----------------
function c = doubleBridgeOp(t)
n = numel(t);
if n < 8, c = threeOptOp(t); return; end
p = sort(randperm(n-1, 3));
A = t(1:p(1)); B = t(p(1)+1:p(2)); C = t(p(2)+1:p(3)); D = t(p(3)+1:end);
c = [A, D, C, B];          % classic double-bridge reconnection
end

% =========================================================================
% Self-contained helpers (operate on a giant-tour permutation).
% =========================================================================

function child = orderCrossoverTour(parent1, parent2)
% Standard Order Crossover (OX): copy a contiguous block from parent1, fill
% the rest from parent2 in its relative order. Preserves good sub-sequences.
n = numel(parent1);
if n < 2, child = parent1; return; end
idx = sort(randperm(n, 2));
i = idx(1); j = idx(2);
child = zeros(1, n);
segment = parent1(i:j);
child(i:j) = segment;
mask = ismember(parent2, segment);
remaining = parent2(~mask);
fillPos = [j+1:n, 1:i-1];
child(fillPos) = remaining(1:numel(fillPos));
end

function outTour = twoOptLocalSearchTour(tour, D, maxPasses)
% Full 2-opt local search on a giant tour using the customer-indexed
% distance matrix D (vectorised over the second edge for speed).
if nargin < 3 || isempty(maxPasses), maxPasses = 30; end
n = numel(tour); szD = size(D);
improved = true; pass = 0;
while improved && pass < maxPasses
    improved = false; pass = pass + 1;
    for i = 1:n-1
        a = tour(i); b = tour(i+1);
        jIdx = (i+1):n;
        cc = tour(jIdx);
        dIdx = mod(jIdx, n) + 1;
        dd = tour(dIdx);
        oldEdges = D(a,b) + D(sub2ind(szD, cc, dd));
        newEdges = D(a,cc) + D(b,dd);
        deltas = newEdges - oldEdges;
        [minDelta, idxMin] = min(deltas);
        if minDelta < -1e-9
            j = jIdx(idxMin);
            tour(i+1:j) = fliplr(tour(i+1:j));
            improved = true;
        end
    end
end
outTour = tour;
end
