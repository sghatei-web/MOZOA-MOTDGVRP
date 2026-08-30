function [coords, demand, Q, depotIdx, info, distMat] = readGolden(filename, whichInstance)
%READGOLDEN  Reads a Golden et al. large-scale CVRP instance file
%   (the "kelly01.txt".."kelly20.txt" / "kelly-pbs.txt" format downloaded
%   from https://neo.lcc.uma.es/vrp/vrp-instances/capacitated-vrp-instances/,
%   NOT the standard CVRPLIB "NAME:/DIMENSION:/NODE_COORD_SECTION" format
%   that readCVRP.m parses).
%
%   [coords, demand, Q, depotIdx, info, distMat] = readGolden(filename, whichInstance)
%
%   This format looks like:
%       [6] with 280 Customers
%        Vehicle capacity = 900      Route-length Limit = 1500
%        Point   --- x ---   --- y ---   ---q---
%           1     30.0000       .0000     10.0000
%           2     29.2478      6.6756     30.0000
%           ...
%   with NO "NAME:"/"DIMENSION:"/"NODE_COORD_SECTION" header keywords at
%   all, so readCVRP.m (which looks for those keywords) returns nothing
%   useful on this file -- this is why loading a renamed "kelly01.vrp"
%   file with loadCVRPInstance.m fails or silently misparses.
%
%   A single file such as "kelly-pbs.txt" may contain SEVERAL instances
%   back-to-back (one block per instance, each starting with its own
%   "[k] with N Customers" header); the single-instance files
%   "kelly01.txt".."kelly20.txt" contain exactly one block each. This
%   function handles both: if the file has more than one block, use
%   whichInstance to select which one (1-based, in file order); if you
%   pass a single-instance file, whichInstance is ignored (or leave it
%   empty / omit it).
%
%   OUTPUTS mirror readCVRP.m exactly, so this can be used as a drop-in
%   inside a loadGoldenInstance-style wrapper (see loadGoldenInstance.m):
%     coords   : (D,2) array of x,y coordinates (depot included)
%     demand   : (D,1) array of demands (0-based row 1 = whichever node
%                the file lists first, which is the depot by convention)
%     Q        : vehicle capacity
%     depotIdx : always 1 (Golden format always lists the depot first)
%     info     : struct with fields Name, Dimension, EdgeWeightType='EUC_2D'
%     distMat  : [] (always empty -- Golden format gives coordinates, not
%                an explicit distance matrix, so distances are computed
%                as Euclidean by the caller, exactly as for EUC_2D
%                CVRPLIB files)
%
%   NOTE ON ROUTE-LENGTH LIMIT: the Golden benchmark's "Route-length
%   Limit" field (a maximum total distance per route) is NOT a concept
%   that exists in this project's MOTDGVRP model (Section 3 of the
%   paper has no route-distance-limit constraint), so it is parsed into
%   info.RouteLengthLimit for reference/logging only and is not enforced
%   anywhere in the solvers. If you want to enforce it, that would need
%   to be added to evaluate.m/decodeRoutes.m separately; this reader
%   only exposes the value.

if nargin < 2, whichInstance = []; end

fid = fopen(filename, 'r');
if fid == -1
    error('readGolden:fileNotFound', 'File "%s" was not found or could not be opened.', filename);
end

raw = fread(fid, '*char')';
fclose(fid);
lines = strsplit(raw, {'\r\n','\n'});

% ---- find every block header line: "[k] with N Customers" or
%      "Data of New LSVRP No. [k] with N Customers" ----
headerIdx = [];
headerNums = [];
for i = 1:numel(lines)
    tok = regexp(lines{i}, '\[(\d+)\]\s*with\s*(\d+)\s*Customers', 'tokens', 'once');
    if ~isempty(tok)
        headerIdx(end+1) = i; %#ok<AGROW>
        headerNums(end+1) = str2double(tok{2}); %#ok<AGROW>
    end
end

if isempty(headerIdx)
    error('readGolden:noInstanceFound', ...
        'No "[k] with N Customers" header found in "%s". Is this really a Golden-format file?', filename);
end

nBlocks = numel(headerIdx);
if nBlocks > 1
    if isempty(whichInstance)
        error('readGolden:multipleInstances', ...
            ['File "%s" contains %d instances. Pass whichInstance (1..%d) ' ...
             'to select which one to read.'], filename, nBlocks, nBlocks);
    end
    if whichInstance < 1 || whichInstance > nBlocks
        error('readGolden:badIndex', 'whichInstance must be between 1 and %d for this file.', nBlocks);
    end
else
    whichInstance = 1;
end

blockStart = headerIdx(whichInstance);
if whichInstance < nBlocks
    blockEnd = headerIdx(whichInstance + 1) - 1;
else
    blockEnd = numel(lines);
end
block = lines(blockStart:blockEnd);

% ---- parse "Vehicle capacity = X   Route-length Limit = Y" ----
Q = NaN; routeLimit = NaN;
for i = 1:numel(block)
    capTok = regexp(block{i}, 'Vehicle\s*capacity\s*=\s*([\d.]+)', 'tokens', 'once');
    if ~isempty(capTok)
        Q = str2double(capTok{1});
    end
    limTok = regexp(block{i}, 'Route-length\s*Limit\s*=\s*([\w.]+)', 'tokens', 'once');
    if ~isempty(limTok)
        if strcmpi(limTok{1}, 'Infinity')
            routeLimit = Inf;
        else
            routeLimit = str2double(limTok{1});
        end
    end
    if ~isnan(Q), break; end  % capacity line always precedes the data rows
end
if isnan(Q)
    error('readGolden:noCapacity', 'Could not find "Vehicle capacity = ..." in the selected block of "%s".', filename);
end

% ---- parse the "Point x y q" data rows ----
pts = [];
for i = 1:numel(block)
    s = strtrim(block{i});
    if isempty(s), continue; end
    % a data row is: <int> <float> <float> <float>  (4 numeric tokens);
    % skip the header/capacity/column-label lines, which don't match this
    vals = sscanf(s, '%f');
    if numel(vals) == 4 && vals(1) == round(vals(1))
        pts(end+1, :) = vals(:)'; %#ok<AGROW>
    end
end

if isempty(pts)
    error('readGolden:noPoints', 'No coordinate/demand rows found in the selected block of "%s".', filename);
end

% sort by point index just in case the file ever lists them out of order
[~, ord] = sort(pts(:,1));
pts = pts(ord, :);

coords = pts(:, 2:3);
demand = pts(:, 4);
depotIdx = 1;  % Golden format always lists the depot as point 1

info = struct();
info.Name = sprintf('Golden-%d-n%d', headerNums(whichInstance), size(pts,1));
info.Dimension = size(pts,1);
info.EdgeWeightType = 'EUC_2D';
info.RouteLengthLimit = routeLimit;  % informational only, not enforced

distMat = [];
end
