function newTour = segmentRelocation(tour)
%SEGMENTRELOCATION  Operator O5 of MOZOA: relocate a contiguous block.
%
%   newTour = segmentRelocation(tour)
%
%   Removes a contiguous block of l in [1, floor(n/5)] cities and reinserts
%   it at another random position. This produces larger structural changes
%   than swap/insertion and helps Pareto-front diversification (Section 4.2,
%   operator O5 of Qat'ei et al.).

n = numel(tour);
if n < 3, newTour = tour; return; end
L = max(1, floor(n/5));
l = randi(L);
i = randi(n - l + 1);
block = tour(i:i+l-1);
rest = tour;
rest(i:i+l-1) = [];
j = randi(numel(rest) + 1);
newTour = [rest(1:j-1), block, rest(j:end)];
end
