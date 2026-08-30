function R = buildReferenceSet(allFronts)
%BUILDREFERENCESET  Pool multiple Pareto fronts and keep only the
%   non-dominated union, for use as the IGD+ reference set (igdPlus.m).
%
%   R = buildReferenceSet(allFronts)
%
%   allFronts : cell array, each cell a K_i x m matrix of objective
%               vectors (e.g. one cell per algorithm's best-run front, or
%               one cell per (algorithm, run) pair for a stricter
%               reference set built from every run rather than just the
%               best one).
%
%   Since the true Pareto front of MOTDGVRP is unknown, this project uses
%   the best-known approximation: the non-dominated subset of every
%   front pooled together, following the same "pool everything, then
%   filter" convention already used to build the shared hypervolume
%   reference point (Section 5).
%
%   USAGE
%     R = buildReferenceSet({mozoaFront, nsga2Front, mlnsga2Front, ...});

pooled = vertcat(allFronts{:});
if isempty(pooled)
    R = pooled;
    return;
end
pooled = unique(round(pooled, 6), 'rows', 'stable');

n = size(pooled, 1);
keep = true(n, 1);
for i = 1:n
    if ~keep(i), continue; end
    for j = 1:n
        if i ~= j && keep(j)
            if all(pooled(j,:) <= pooled(i,:)) && any(pooled(j,:) < pooled(i,:))
                keep(i) = false;
                break;
            end
        end
    end
end
R = pooled(keep, :);
end
