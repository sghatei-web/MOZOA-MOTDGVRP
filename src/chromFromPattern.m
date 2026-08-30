function chrom = chromFromPattern(patterns, prob)
%CHROMFROMPATTERN  Build an individual embedding a frequent sub-sequence.
%
%   chrom = chromFromPattern(patterns, prob)
%
%   Used in MLNSGA-2 to generate the half of the next population that is
%   derived from frequent contiguous sequences (Section 5.2.2). One pattern
%   is chosen at random (longer patterns favoured because the list is sorted
%   by length), placed at a random position, and the remaining customers are
%   filled in random order.

n = prob.n;
allCust = 2:(n+1);

if isempty(patterns)
    chrom = allCust(randperm(n));
    return;
end

% Favour longer patterns: pick from the top of the (length-sorted) list
topK = max(1, round(numel(patterns)*0.5));
p = patterns{randi(topK)};

% Remaining customers not in the pattern, shuffled
rest = allCust(~ismember(allCust, p));
rest = rest(randperm(numel(rest)));

% Insert the pattern block at a random split of the remaining sequence
cut = randi(numel(rest)+1) - 1;
chrom = [rest(1:cut), p, rest(cut+1:end)];

% Guard against accidental duplication / loss
if numel(unique(chrom)) ~= n
    chrom = allCust(randperm(n));
end
end
