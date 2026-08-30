function [coords, demand, Q, depotIdx, info, distMat] = readCVRP(filename)
%READCVRP  Read a CVRPLIB/TSPLIB .vrp instance (E-n51-k5.vrp, etc.).
%
%   [coords, demand, Q, depotIdx, info, distMat] = readCVRP(filename)
%
%   Supports two edge-weight conventions used by Christofides Set E:
%     * EUC_2D  : NODE_COORD_SECTION with x,y coordinates (distances are
%                 Euclidean and computed by buildInstance).
%     * EXPLICIT: EDGE_WEIGHT_SECTION with an explicit distance matrix in
%                 LOWER_ROW / UPPER_ROW / FULL_MATRIX format (no coordinates).
%
%   OUTPUTS
%     coords   : D x 2 coordinates if EUC_2D, else [] (no coordinates given).
%     demand   : D x 1 demands (depot demand forced to 0).
%     Q        : vehicle capacity.
%     depotIdx : 1-based depot node index from DEPOT_SECTION (default 1).
%     info     : struct with Name, Type, Dimension, EdgeWeightType,
%                EdgeWeightFormat, Capacity.
%     distMat  : D x D distance matrix if EXPLICIT, else [] (caller computes
%                Euclidean distances from coords).
%
%   The reader matches coordinates and demands to node ids, so it is robust
%   to ordering and to the depot not being listed first.

fid = fopen(filename,'r');
if fid == -1
    error('readCVRP:fileNotFound','Cannot open "%s".', filename);
end

info = struct('Name','','Type','','Dimension',0, ...
              'EdgeWeightType','EUC_2D','EdgeWeightFormat','','Capacity',0);
section = '';
nodeIds = []; xy = [];
demIds  = []; dem = [];
ewVals  = [];
depotIdx = 1;

line = fgetl(fid);
while ischar(line)
    s = strtrim(line);
    if isempty(s), line = fgetl(fid); continue; end
    up = upper(s);

    % --- section headers ---
    if contains(up,'NODE_COORD_SECTION'), section = 'coord';  line = fgetl(fid); continue; end
    if contains(up,'DEMAND_SECTION'),     section = 'demand'; line = fgetl(fid); continue; end
    if contains(up,'EDGE_WEIGHT_SECTION'),section = 'ew';     line = fgetl(fid); continue; end
    if contains(up,'DEPOT_SECTION'),      section = 'depot';  line = fgetl(fid); continue; end
    if contains(up,'DISPLAY_DATA_SECTION'),section = 'skip';  line = fgetl(fid); continue; end
    if strcmp(up,'EOF'), break; end

    % --- key : value spec lines (only when not inside a data section) ---
    cp = strfind(s, ':');
    if ~isempty(cp) && (isempty(section) || any(strcmp(section,{'skip'})))
        key = strtrim(s(1:cp(1)-1));
        val = strtrim(s(cp(1)+1:end));
        switch upper(key)
            case 'NAME',               info.Name = val;
            case 'TYPE',               info.Type = val;
            case 'DIMENSION',          info.Dimension = str2double(val);
            case 'EDGE_WEIGHT_TYPE',   info.EdgeWeightType = val;
            case 'EDGE_WEIGHT_FORMAT', info.EdgeWeightFormat = val;
            case 'CAPACITY',           info.Capacity = str2double(val);
        end
        line = fgetl(fid); continue;
    end

    % --- data rows ---
    nums = sscanf(s, '%f')';
    switch section
        case 'coord'
            if numel(nums) >= 3
                nodeIds(end+1) = nums(1);   %#ok<AGROW>
                xy(end+1,:)    = nums(2:3);  %#ok<AGROW>
            end
        case 'demand'
            if numel(nums) >= 2
                demIds(end+1) = nums(1);    %#ok<AGROW>
                dem(end+1)    = nums(2);    %#ok<AGROW>
            end
        case 'ew'
            if ~isempty(nums)
                ewVals = [ewVals, nums]; %#ok<AGROW>
            end
        case 'depot'
            if ~isempty(nums) && nums(1) > 0
                depotIdx = nums(1);
            end
    end
    line = fgetl(fid);
end
fclose(fid);

D = info.Dimension;
if D == 0
    % infer dimension if header omitted it
    D = max(numel(nodeIds), numel(demIds));
end

% --- demands aligned by id (ids assumed 1..D when coords absent) ---
demand = zeros(D,1);
if ~isempty(nodeIds)
    [nodeIds, order] = sort(nodeIds);
    xy = xy(order,:);
    for k = 1:numel(demIds)
        pos = find(nodeIds == demIds(k), 1);
        if ~isempty(pos), demand(pos) = dem(k); end
    end
    depotPos = find(nodeIds == depotIdx, 1);
    coords = xy;
else
    % EXPLICIT: no coordinates; demand ids are 1..D in order
    for k = 1:numel(demIds)
        id = demIds(k);
        if id >= 1 && id <= D, demand(id) = dem(k); end
    end
    depotPos = depotIdx;
    coords = [];
end
if isempty(depotPos), depotPos = 1; end
depotIdx = depotPos;

% --- build explicit distance matrix if provided ---
distMat = [];
if strcmpi(strtrim(info.EdgeWeightType),'EXPLICIT') && ~isempty(ewVals)
    distMat = expandMatrix(ewVals, D, upper(strtrim(info.EdgeWeightFormat)));
end

Q = info.Capacity;
end

% -------------------------------------------------------------------------
function M = expandMatrix(vals, D, fmt)
%EXPANDMATRIX  Reconstruct a full symmetric D x D matrix from a flat list.
M = zeros(D,D);
idx = 1;
switch fmt
    case {'LOWER_ROW','LOWERROW'}
        % strictly lower triangle, row by row (no diagonal)
        for r = 2:D
            for c = 1:r-1
                M(r,c) = vals(idx); M(c,r) = vals(idx); idx = idx + 1;
            end
        end
    case {'UPPER_ROW','UPPERROW'}
        % strictly upper triangle, row by row (no diagonal)
        for r = 1:D-1
            for c = r+1:D
                M(r,c) = vals(idx); M(c,r) = vals(idx); idx = idx + 1;
            end
        end
    case {'LOWER_DIAG_ROW','LOWERDIAGROW'}
        for r = 1:D
            for c = 1:r
                M(r,c) = vals(idx); M(c,r) = vals(idx); idx = idx + 1;
            end
        end
    case {'UPPER_DIAG_ROW','UPPERDIAGROW'}
        for r = 1:D
            for c = r:D
                M(r,c) = vals(idx); M(c,r) = vals(idx); idx = idx + 1;
            end
        end
    case {'FULL_MATRIX','FULLMATRIX'}
        for r = 1:D
            for c = 1:D
                M(r,c) = vals(idx); idx = idx + 1;
            end
        end
    otherwise
        % default to LOWER_ROW
        for r = 2:D
            for c = 1:r-1
                M(r,c) = vals(idx); M(c,r) = vals(idx); idx = idx + 1;
            end
        end
end
end
