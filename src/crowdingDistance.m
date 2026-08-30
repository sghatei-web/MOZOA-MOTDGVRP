function cd = crowdingDistance(F)
%CROWDINGDISTANCE  NSGA-II crowding distance for one front.
%
%   cd = crowdingDistance(F)
%   F : m x M objective matrix for the individuals of a single front.
%   cd: m x 1 crowding distances (Inf for boundary solutions).

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
