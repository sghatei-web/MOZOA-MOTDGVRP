function patterns = mineContiguousPatterns(elite, lambda, maxLen)
%MINECONTIGUOUSPATTERNS  Frequent contiguous sub-sequences (Apriori-style).
%
%   patterns = mineContiguousPatterns(elite, lambda, maxLen)
%
%   Implements the pattern-discovery step of MLNSGA-2 (Section 5.2.2). Given
%   a set of elite (non-dominated) chromosomes, it extracts contiguous
%   sub-sequences (sub-routes) whose frequency across the elite set is at
%   least lambda. The Apriori principle is used: a length-(k+1) candidate can
%   only be frequent if its length-k prefix is frequent, so we grow patterns
%   incrementally and prune.
%
%   INPUTS
%     elite  : cell array of chromosomes (each a permutation row vector).
%     lambda : minimum number of occurrences (support threshold).
%     maxLen : maximum pattern length to mine (default 4).
%
%   OUTPUT
%     patterns : cell array of frequent contiguous sequences (row vectors),
%                sorted by descending length then descending support.

if nargin < 3 || isempty(maxLen), maxLen = 4; end
if isempty(elite), patterns = {}; return; end

% ---- length-2 contiguous pairs and their support --------------------------
countMap = containers.Map('KeyType','char','ValueType','double');
seqMap   = containers.Map('KeyType','char','ValueType','any');

% Seed with length-2
[countMap, seqMap, freqKeys] = countContiguous(elite, 2, countMap, seqMap, lambda);

allPatterns = {};
support     = [];
for k = 1:numel(freqKeys)
    seq = seqMap(freqKeys{k});
    allPatterns{end+1} = seq;           %#ok<AGROW>
    support(end+1) = countMap(freqKeys{k}); %#ok<AGROW>
end

% ---- grow longer patterns by extending frequent ones ----------------------
prevFreq = freqKeys;
for len = 3:maxLen
    if isempty(prevFreq), break; end
    cMap = containers.Map('KeyType','char','ValueType','double');
    sMap = containers.Map('KeyType','char','ValueType','any');
    [cMap, sMap, fk] = countContiguous(elite, len, cMap, sMap, lambda, prevFreq);
    for k = 1:numel(fk)
        seq = sMap(fk{k});
        allPatterns{end+1} = seq;       %#ok<AGROW>
        support(end+1) = cMap(fk{k});   %#ok<AGROW>
    end
    prevFreq = fk;
end

% ---- sort by length desc, then support desc -------------------------------
if isempty(allPatterns)
    patterns = {};
    return;
end
lens = cellfun(@numel, allPatterns);
[~, ord] = sortrows([-lens(:), -support(:)]);
patterns = allPatterns(ord);
end

% -------------------------------------------------------------------------
function [countMap, seqMap, freqKeys] = countContiguous(elite, len, countMap, seqMap, lambda, prevFreq)
% Count contiguous sub-sequences of given length across elite chromosomes.
% If prevFreq is supplied, only extend sequences whose length-(len-1) prefix
% is frequent (Apriori pruning).
if nargin < 6, prevFreq = []; end
usePrune = ~isempty(prevFreq);
prevSet = containers.Map('KeyType','char','ValueType','logical');
for k = 1:numel(prevFreq), prevSet(prevFreq{k}) = true; end

for e = 1:numel(elite)
    c = elite{e};
    for s = 1:numel(c)-len+1
        seq = c(s:s+len-1);
        if usePrune
            prefixKey = key(seq(1:end-1));
            if ~isKey(prevSet, prefixKey), continue; end
        end
        kk = key(seq);
        if isKey(countMap, kk)
            countMap(kk) = countMap(kk) + 1;
        else
            countMap(kk) = 1;
            seqMap(kk) = seq;
        end
    end
end

freqKeys = {};
if nargin >= 5 && ~isempty(lambda)
    allKeys = countMap.keys;
    for k = 1:numel(allKeys)
        if countMap(allKeys{k}) >= lambda
            freqKeys{end+1} = allKeys{k}; %#ok<AGROW>
        end
    end
end
end

function k = key(seq)
k = strtrim(sprintf('%d_', seq));
end
