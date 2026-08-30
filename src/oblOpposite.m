function opp = oblOpposite(chrom)
%OBLOPPOSITE  Clockwise (CW) opposite tour, Ergezer & Simon (2011), Eq. (40).
%
%   opp = oblOpposite(chrom)
%
%   Builds the opposite of a permutation by the CW construction used for
%   Opposition-Based Learning in MLNSGA-2 (Section 5.2.1). For a tour with
%   n customers numbered along positions 1..n, the opposite interleaves the
%   first and second halves:
%       T_CW = [1, 1+n/2, 2, 2+n/2, ..., n/2, n]
%   When n is odd, a dummy position is appended, the rule is applied, and the
%   dummy is removed afterwards.

c = chrom(:)';
n = numel(c);

oddPad = false;
if mod(n,2) == 1
    n = n + 1;
    oddPad = true;
    padPos = n;                  % dummy position index
end

half = n/2;
order = zeros(1,n);
idx = 1;
for k = 1:half
    order(idx) = k;        idx = idx + 1;     % position k (first half)
    order(idx) = k + half; idx = idx + 1;     % position k+n/2 (second half)
end

if oddPad
    order(order == padPos) = [];              % drop the dummy position
    n = n - 1;
end

% Map positions back to customers of the original chromosome
opp = c(order);
end
