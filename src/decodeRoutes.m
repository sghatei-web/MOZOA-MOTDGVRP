function routes = decodeRoutes(chrom, prob)
%DECODEROUTES  Split a customer permutation ("giant tour") into routes.
%
%   routes = decodeRoutes(chrom, prob)
%
%   The chromosome is a permutation of customer indices (2..n+1, i.e. the
%   depot at index 1 is excluded). Customers are appended to the current
%   route until adding the next one would exceed capacity Q, at which point
%   a new route (vehicle) starts. This is the standard split used with the
%   route-as-sequence representation described in Section 5.1 of the paper.
%
%   OUTPUT  routes : cell array, each cell a row vector of 1-based node
%                    indices WITHOUT depot bookends (depot is added by the
%                    evaluator). Number of routes <= prob.NV for feasibility.

chrom = chrom(:)';
routes = {};
cur = [];
load = 0;
for c = chrom
    q = prob.demand(c);
    if load + q > prob.Q && ~isempty(cur)
        routes{end+1} = cur;   %#ok<AGROW>
        cur = [];
        load = 0;
    end
    cur(end+1) = c;            %#ok<AGROW>
    load = load + q;
end
if ~isempty(cur)
    routes{end+1} = cur;       %#ok<AGROW>
end
end
