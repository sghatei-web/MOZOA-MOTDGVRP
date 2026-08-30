function hv = hypervolume(F, ref)
%HYPERVOLUME  Hypervolume dominated by front F w.r.t. reference point ref.
%
%   hv = hypervolume(F, ref)
%
%   F   : K x m matrix of objective vectors (minimisation).
%   ref : 1 x m reference (nadir) point; must be dominated by all of F
%         (i.e. larger in every objective). Points beyond ref are clipped.
%
%   Supports m = 2 (exact sweep) and m = 3 (slice-based exact computation).
%   For m = 3 the front is sorted on f3 and the 2-D HV of the accumulated
%   projection is integrated across the third axis. Used for R_HV (Eq. 12)
%   and for reporting.

if isempty(F), hv = 0; return; end
m = size(F,2);

% keep only points that dominate the reference
keep = all(F < ref, 2);
F = F(keep,:);
if isempty(F), hv = 0; return; end

switch m
    case 2
        hv = hv2d(F, ref);
    case 3
        hv = hv3d(F, ref);
    otherwise
        error('hypervolume: only 2 or 3 objectives supported.');
end
end

% -------------------------------------------------------------------------
function hv = hv2d(F, ref)
% Exact 2-D hypervolume for minimisation.
% Take the non-dominated subset, sort by f1 ascending (=> f2 descending),
% sum rectangle areas.
F = paretoFilter(F);
[~, idx] = sort(F(:,1), 'ascend');
F = F(idx,:);
hv = 0;
prevF2 = ref(2);
for i = 1:size(F,1)
    width  = ref(1) - F(i,1);
    height = prevF2 - F(i,2);
    if width > 0 && height > 0
        hv = hv + width * height;
    end
    prevF2 = min(prevF2, F(i,2));
end
end

% -------------------------------------------------------------------------
function hv = hv3d(F, ref)
% Exact 3-D hypervolume by slab sweep along the third objective.
% Sort points by f3 ascending. Between consecutive z-levels, the dominated
% area equals the 2-D HV of all points whose f3 is <= the lower z-level.
F = paretoFilter(F);
[~, idx] = sort(F(:,3), 'ascend');
F = F(idx,:);
K = size(F,1);
hv = 0;
zlevels = [F(:,3); ref(3)];           % next-larger z for each point (last = ref)
for i = 1:K
    dz = zlevels(i+1) - F(i,3);
    if dz <= 0, continue; end
    proj = F(1:i, 1:2);               % points with f3 <= current level
    area = hv2d(proj, ref(1:2));
    hv = hv + area * dz;
end
end

% -------------------------------------------------------------------------
function P = paretoFilter(F)
% Return the non-dominated rows of F (minimisation).
n = size(F,1);
keep = true(n,1);
for i = 1:n
    if ~keep(i), continue; end
    for j = 1:n
        if i~=j && keep(j)
            if all(F(j,:) <= F(i,:)) && any(F(j,:) < F(i,:))
                keep(i) = false; break;
            end
        end
    end
end
P = F(keep,:);
end
