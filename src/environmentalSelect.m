function keep = environmentalSelect(F, N)
%ENVIRONMENTALSELECT  NSGA-II survivor selection.
%
%   keep = environmentalSelect(F, N)
%
%   Given objective matrix F (P x M), return indices of the N survivors:
%   fill whole fronts in order; for the last partially-fitting front, choose
%   the most widely spread individuals by crowding distance.

[fronts, ~] = fastNonDominatedSort(F);
keep = [];
for f = 1:numel(fronts)
    fr = fronts{f};
    if numel(keep) + numel(fr) <= N
        keep = [keep, fr]; %#ok<AGROW>
    else
        cd = crowdingDistance(F(fr,:));
        [~, ord] = sort(cd, 'descend');
        need = N - numel(keep);
        keep = [keep, fr(ord(1:need))]; %#ok<AGROW>
        break;
    end
end
end
