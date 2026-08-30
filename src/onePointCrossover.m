function child = onePointCrossover(p1, p2)
%ONEPOINTCROSSOVER  Order-based one-point crossover for permutations.
%
%   child = onePointCrossover(p1, p2)
%
%   The "1-point" crossover described in Section 5.1 of the paper, adapted to
%   permutation chromosomes: take the prefix of p1 up to a random cut point,
%   then append the remaining customers in the order they appear in p2. This
%   guarantees a valid permutation (each customer exactly once).

n = numel(p1);
cut = randi([1, n-1]);
head = p1(1:cut);
tail = p2(~ismember(p2, head));
child = [head, tail];
end
