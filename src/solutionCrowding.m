function [cdVal, maxCd, rankVal] = solutionCrowding(obj, archiveObjs)
%SOLUTIONCROWDING  Crowding distance & non-dominated rank of one solution.
%
%   [cdVal, maxCd, rankVal] = solutionCrowding(obj, archiveObjs)
%
%   obj          : 1 x m objective vector of the solution of interest.
%   archiveObjs  : K x m matrix of current archive objective vectors.
%
%   Returns the crowding distance cdVal of obj within the combined set
%   (archive U {obj}), the maximum crowding distance maxCd over the archive
%   (used to normalise R_DIV, Eq. 13), and the non-dominated rank rankVal of
%   obj relative to the archive (1 = non-dominated). Boundary/extreme points
%   yield Inf crowding, which is capped to 1 for the bounded state feature.

if isempty(archiveObjs)
    cdVal = 1; maxCd = 1; rankVal = 1;
    return;
end

F = [archiveObjs; obj];
idxSelf = size(F,1);

cd = crowdingVec(F);
cdVal = cd(idxSelf);

cdArch = crowdingVec(archiveObjs);
finite = cdArch(~isinf(cdArch));
if isempty(finite), maxCd = 1; else, maxCd = max([finite; eps]); end

% non-dominated rank of obj vs archive
rankVal = 1;
for i = 1:size(archiveObjs,1)
    if dominates(archiveObjs(i,:), obj)
        rankVal = rankVal + 1;
    end
end
end

% -------------------------------------------------------------------------
function cd = crowdingVec(F)
[m, M] = size(F);
cd = zeros(m,1);
if m <= 2, cd(:) = inf; return; end
for k = 1:M
    [vals, idx] = sort(F(:,k));
    cd(idx(1)) = inf; cd(idx(end)) = inf;
    rng = vals(end) - vals(1); if rng == 0, rng = 1; end
    for i = 2:m-1
        cd(idx(i)) = cd(idx(i)) + (vals(i+1) - vals(i-1)) / rng;
    end
end
end
