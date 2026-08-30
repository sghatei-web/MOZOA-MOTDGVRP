function chrom = clarkeWrightChrom(prob)
%CLARKEWRIGHTCHROM  Clarke & Wright (1964) savings heuristic -> giant tour.
%
%   chrom = clarkeWrightChrom(prob)
%
%   Classic parallel savings construction respecting capacity Q, then the
%   resulting routes are concatenated into a single customer permutation to
%   seed the population (Section 5.1, Clarke & Wright heuristic).

n = prob.n;
D = prob.dist;
depot = 1;
cust = 2:(n+1);

% Savings s(i,j) = d(0,i) + d(0,j) - d(i,j)
S = [];
for a = 1:numel(cust)
    for b = a+1:numel(cust)
        i = cust(a); j = cust(b);
        s = D(depot,i) + D(depot,j) - D(i,j);
        S(end+1,:) = [s, i, j]; %#ok<AGROW>
    end
end
S = sortrows(S, -1);              % descending savings

% Each customer starts in its own route
route = num2cell(cust);          % cell of routes
routeOf = containers.Map('KeyType','double','ValueType','double');
for r = 1:numel(route), routeOf(route{r}) = r; end
load = arrayfun(@(c) prob.demand(c), cust);
load = containers.Map(num2cell(cust), num2cell(load));

for s = 1:size(S,1)
    i = S(s,2); j = S(s,3);
    if ~isKey(routeOf,i) || ~isKey(routeOf,j), continue; end
    ri = routeOf(i); rj = routeOf(j);
    if ri == rj || isempty(route{ri}) || isempty(route{rj}), continue; end
    Ri = route{ri}; Rj = route{rj};
    % i must be at an end of Ri, j at an end of Rj (interior merge only)
    if Ri(end) ~= i && Ri(1) ~= i, continue; end
    if Rj(1) ~= j && Rj(end) ~= j, continue; end
    if load(i) + load(j) > prob.Q, continue; end %#ok<*NASGU>
    % merge so that i and j become adjacent
    if Ri(end) ~= i, Ri = fliplr(Ri); end
    if Rj(1)  ~= j, Rj = fliplr(Rj); end
    merged = [Ri, Rj];
    newLoad = load(i) + load(j);
    route{ri} = merged; route{rj} = [];
    for c = merged, routeOf(c) = ri; end
    for c = merged, load(c) = newLoad; end
end

chrom = [];
for r = 1:numel(route)
    if ~isempty(route{r}), chrom = [chrom, route{r}]; end %#ok<AGROW>
end
% Safety: ensure it is a valid permutation of all customers
if numel(chrom) ~= n || numel(unique(chrom)) ~= n
    chrom = cust;                % fall back to identity
end
end
