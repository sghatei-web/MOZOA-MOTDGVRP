function chrom = nearestNeighborChrom(prob)
%NEARESTNEIGHBORCHROM  Nearest-neighbour giant tour over all customers.
%
%   chrom = nearestNeighborChrom(prob)
%
%   Builds a single permutation of customers (indices 2..n+1) by repeatedly
%   moving to the nearest unvisited customer, starting from the depot. Used
%   to seed the initial population (Section 5.1, Solomon 1987 heuristic).

n = prob.n;
unvisited = 2:(n+1);
chrom = zeros(1,n);
current = 1;                      % depot
for k = 1:n
    d = prob.dist(current, unvisited);
    [~, idx] = min(d);
    nxt = unvisited(idx);
    chrom(k) = nxt;
    current = nxt;
    unvisited(idx) = [];
end
end
