function [coords, demand, Q, depotIdx, info, distMat] = readGoldenSingle(filename)
%READGOLDENSINGLE  Reads a Golden et al. single-instance CVRP file in the
%   format actually used by the individually-distributed "kelly01.txt"
%   .. "kelly20.txt" files (as opposed to the combined "kelly-pbs.txt",
%   which uses a different layout -- see readGolden.m for that one).
%
%   [coords, demand, Q, depotIdx, info, distMat] = readGoldenSingle(filename)
%
%   This format has NO section keywords at all -- not even the
%   "[k] with N Customers" / "Vehicle capacity = ..." labels of
%   kelly-pbs.txt. It looks like:
%
%       240 550 650 0 5646.46
%          0.0000      0.0000
%         30.0000      0.0000     10
%         29.6306      4.6930     30
%         28.5317      9.2705     30
%         ...
%
%   Line 1: <numCustomers> <capacity> <routeLengthLimit> <flag> <bestKnown>
%           (flag and bestKnown are parsed into info but not used further;
%           routeLengthLimit is parsed into info.RouteLengthLimit for
%           reference only -- see the note on this in loadGoldenInstance.m)
%   Line 2: <depotX> <depotY> [optionally a third value, the depot's
%           demand, which is always 0 when present -- some files include
%           this column for the depot line and some do not, both are
%           handled here]
%   Lines 3..end: <customerX> <customerY> <demand>  (ONE row per customer,
%           note: unlike kelly-pbs.txt, there is NO leading point-index
%           column here -- rows are implicitly numbered in file order)
%
%   OUTPUTS mirror readCVRP.m/readGolden.m:
%     coords   : (D,2) array, row 1 = depot, rows 2..D = customers in file order
%     demand   : (D,1) array, demand(1) = 0 (the depot)
%     Q        : vehicle capacity
%     depotIdx : always 1
%     info     : struct with fields Name, Dimension, EdgeWeightType='EUC_2D',
%                RouteLengthLimit, BestKnown
%     distMat  : [] (coordinates are given; Euclidean distance is computed
%                by the caller, as for EUC_2D CVRPLIB files)

fid = fopen(filename, 'r');
if fid == -1
    error('readGoldenSingle:fileNotFound', 'File "%s" was not found or could not be opened.', filename);
end
raw = fread(fid, '*char')';
fclose(fid);
lines = strsplit(strtrim(raw), {'\r\n','\n'});
lines = lines(~cellfun(@(s) isempty(strtrim(s)), lines));  % drop blank lines

if numel(lines) < 2
    error('readGoldenSingle:tooShort', 'File "%s" does not have enough lines to be a valid instance.', filename);
end

% ---- line 1: numCustomers capacity routeLimit flag bestKnown ----
headerVals = sscanf(lines{1}, '%f');
if numel(headerVals) < 3
    error('readGoldenSingle:badHeader', ...
        'First line of "%s" does not look like "<numCustomers> <capacity> <routeLimit> ...": "%s"', ...
        filename, lines{1});
end
numCustomers = round(headerVals(1));
Q = headerVals(2);
routeLimit = headerVals(3);
if routeLimit >= 999999
    routeLimit = Inf;  % this benchmark's convention for "no route-length limit"
end
bestKnown = NaN;
if numel(headerVals) >= 5
    bestKnown = headerVals(5);
end

% ---- line 2: depot coordinates (2 or 3 numbers; 3rd, if present, is
%      always a demand of 0 and is simply discarded) ----
depotVals = sscanf(lines{2}, '%f');
if numel(depotVals) < 2
    error('readGoldenSingle:badDepotLine', ...
        'Second line of "%s" does not look like depot coordinates: "%s"', filename, lines{2});
end
depotXY = depotVals(1:2)';

% ---- remaining lines: customer x y demand (exactly 3 numbers each) ----
custLines = lines(3:end);
pts = zeros(numel(custLines), 3);
nParsed = 0;
for i = 1:numel(custLines)
    vals = sscanf(custLines{i}, '%f');
    if numel(vals) == 3
        nParsed = nParsed + 1;
        pts(nParsed, :) = vals(:)';
    end
end
pts = pts(1:nParsed, :);

if nParsed ~= numCustomers
    warning('readGoldenSingle:countMismatch', ...
        ['Header of "%s" declares %d customers, but %d customer data rows ' ...
         'were parsed. Proceeding with the %d rows actually found.'], ...
        filename, numCustomers, nParsed, nParsed);
end

coords = [depotXY; pts(:,1:2)];
demand = [0; pts(:,3)];
depotIdx = 1;

[~, fbase, ~] = fileparts(filename);
info = struct();
info.Name = sprintf('Golden-%s-n%d', fbase, size(coords,1) - 1);
info.Dimension = size(coords,1);
info.EdgeWeightType = 'EUC_2D';
info.RouteLengthLimit = routeLimit;
info.BestKnown = bestKnown;

distMat = [];
end
