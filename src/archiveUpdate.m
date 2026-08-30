function [A, entered] = archiveUpdate(A, cand, Amax)
%ARCHIVEUPDATE  Insert a candidate into the external Pareto archive.
%
%   [A, entered] = archiveUpdate(A, cand, Amax)
%
%   Implements Pareto Archive Management (Section 4.4 of Qat'ei et al.):
%     (1) if cand is dominated by any archive member -> reject;
%     (2) remove archive members dominated by cand;
%     (3) add cand;
%     (4) if |A| > Amax, drop the minimum-crowding-distance member (Eq. 16).
%
%   A is a struct array with fields:
%       .obj   (1 x m objective vector)
%       .chrom (permutation)
%   cand has the same fields. entered is true iff cand was inserted.

entered = false;

% Empty archive: accept directly
if isempty(A)
    A = cand;
    entered = true;
    return;
end

objs = vertcat(A.obj);

% (1) dominated by an existing member?
for i = 1:size(objs,1)
    if dominates(objs(i,:), cand.obj) || isequal(objs(i,:), cand.obj)
        return;   % reject (not better than what we have)
    end
end

% (2) remove members dominated by cand
domByCand = false(numel(A),1);
for i = 1:numel(A)
    if dominates(cand.obj, A(i).obj)
        domByCand(i) = true;
    end
end
A(domByCand) = [];

% (3) add candidate
A(end+1) = cand;
entered = true;

% (4) truncate by crowding distance if oversized
if numel(A) > Amax
    objs = vertcat(A.obj);
    cd = archiveCrowding(objs);
    [~, worst] = min(cd);
    A(worst) = [];
end
end

% -------------------------------------------------------------------------
function cd = archiveCrowding(F)
% Crowding distance per Eq. (16): sum over objectives of normalised gaps
% between neighbouring solutions; boundary solutions get Inf.
[m, M] = size(F);
cd = zeros(m,1);
if m <= 2
    cd(:) = inf;
    return;
end
for k = 1:M
    [vals, idx] = sort(F(:,k));
    cd(idx(1))   = inf;
    cd(idx(end)) = inf;
    rng = vals(end) - vals(1);
    if rng == 0, rng = 1; end
    for i = 2:m-1
        cd(idx(i)) = cd(idx(i)) + (vals(i+1) - vals(i-1)) / rng;
    end
end
end
