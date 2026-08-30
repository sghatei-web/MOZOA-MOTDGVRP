function c = shiftMutation(chrom)
%SHIFTMUTATION  Relocate one customer to a new random position ("shift").
%   One of the two mutation operators used in the paper's local search.
c = chrom;
n = numel(c);
if n < 2, return; end
i = randi(n);
gene = c(i);
c(i) = [];
j = randi(numel(c)+1);
c = [c(1:j-1), gene, c(j:end)];
end
