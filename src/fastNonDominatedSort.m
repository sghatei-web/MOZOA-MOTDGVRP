function [fronts, rank] = fastNonDominatedSort(F)
%FASTNONDOMINATEDSORT  Deb et al. (2002) NSGA-II non-dominated sorting.
%
%   [fronts, rank] = fastNonDominatedSort(F)
%
%   F     : P x M matrix of objective vectors (rows = individuals).
%   fronts: cell array; fronts{k} = indices of individuals in front k.
%   rank  : P x 1 vector with the front number of each individual.

P = size(F,1);
S = cell(P,1);          % set of solutions each one dominates
nDom = zeros(P,1);      % domination counter
rank = zeros(P,1);
fronts = {[]};

for p = 1:P
    S{p} = [];
    for q = 1:P
        if p == q, continue; end
        if dominates(F(p,:), F(q,:))
            S{p}(end+1) = q;
        elseif dominates(F(q,:), F(p,:))
            nDom(p) = nDom(p) + 1;
        end
    end
    if nDom(p) == 0
        rank(p) = 1;
        fronts{1}(end+1) = p; %#ok<AGROW>
    end
end

k = 1;
while ~isempty(fronts{k})
    next = [];
    for p = fronts{k}
        for q = S{p}
            nDom(q) = nDom(q) - 1;
            if nDom(q) == 0
                rank(q) = k + 1;
                next(end+1) = q; %#ok<AGROW>
            end
        end
    end
    k = k + 1;
    fronts{k} = next;
end
fronts(end) = [];       % remove trailing empty front
end
