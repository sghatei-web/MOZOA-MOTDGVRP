function c = twoOptMutation(chrom)
%TWOOPTMUTATION  Reverse a random segment of the chromosome (2-opt move).
%   The second mutation operator used in the paper's local-search mutation.
c = chrom;
n = numel(c);
if n < 2, return; end
i = randi(n); j = randi(n);
if i > j, [i,j] = deal(j,i); end
c(i:j) = c(j:-1:i);
end
