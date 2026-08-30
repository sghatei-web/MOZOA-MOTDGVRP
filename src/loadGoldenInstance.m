function prob = loadGoldenInstance(filename, whichInstance, catMode, seed)
%LOADGOLDENINSTANCE  MOTDGVRP instance from a Golden et al. large-scale
%   CVRP benchmark file, auto-detecting which of the TWO distinct
%   Golden-benchmark text formats this project has encountered is in use:
%
%     (a) the "combined" format used by kelly-pbs.txt: one or more blocks,
%         each starting with a "[k] with N Customers" / "Data of New
%         LSVRP No. [k] with N Customers" header line, followed by
%         "Vehicle capacity = ... Route-length Limit = ..." and then rows
%         of "<pointIndex> <x> <y> <demand>". Parsed by readGolden.m.
%
%     (b) the "single-instance" format used by the individually
%         distributed kelly01.txt..kelly20.txt files: NO section keywords
%         at all -- just a header line "<numCustomers> <capacity>
%         <routeLimit> <flag> <bestKnown>", then one depot coordinate
%         line, then one "<x> <y> <demand>" row per customer (no leading
%         point-index column). Parsed by readGoldenSingle.m.
%
%   Neither format is the standard CVRPLIB "NAME:/DIMENSION:/
%   NODE_COORD_SECTION" format that readCVRP.m/loadCVRPInstance.m expect,
%   which is why those functions cannot read Golden-benchmark files
%   directly (renaming to .vrp does not help -- the content differs).
%
%   prob = loadGoldenInstance(filename, whichInstance, catMode, seed)
%
%   Reads the file with whichever of readGolden.m/readGoldenSingle.m
%   matches its format, then attaches the time-dependent green layer (3
%   slots, speed Table 3, GPH Table 4, load surcharge Eq. 39) via
%   buildInstance -- exactly the same pipeline loadCVRPInstance.m uses
%   for the Christofides Set E instances, so every solver (solveZOA8op,
%   solveMOTDGVRP, solveSPEA2, solveMOEAD) and run30_mozoa_paper.m work on
%   a Golden instance with NO code changes.
%
%   INPUTS
%     filename      path to the .txt file (e.g. 'instances_large/kelly01.txt')
%     whichInstance ONLY relevant for combined-format files (kelly-pbs.txt):
%                   which instance block to read if the file contains
%                   several back-to-back (1-based, in file order). Ignored
%                   for single-instance files; leave [] or omit for those.
%     catMode       'uniform' (paper default) or 'conflict' (default: 'conflict')
%     seed          RNG seed for the road-category layout (default: 7)
%
%   NOTE ON VEHICLE COUNT: Golden-format files specify only a vehicle
%   CAPACITY, not a fixed number of vehicles NV (unlike E-n**-k* files,
%   whose name encodes NV). NV is therefore estimated here as
%   ceil(sum(demand)/Q) -- the minimum number of vehicles that could
%   possibly cover total demand -- exactly the same fallback
%   loadCVRPInstance.m itself uses whenever a name has no "-k<NV>" suffix.
%   If you want a specific NV (e.g. matching a published best-known
%   solution's vehicle count), set prob.NV manually after this call.
%
%   USAGE
%     prob = loadGoldenInstance('instances_large/kelly01.txt');       % single-file format
%     prob = loadGoldenInstance('instances_large/kelly-pbs.txt', 6);  % combined format, block [6]

if nargin < 2, whichInstance = []; end
if nargin < 3 || isempty(catMode), catMode = 'conflict'; end
if nargin < 4 || isempty(seed),    seed = 7;             end

if isGoldenCombinedFormat(filename)
    [coords, demand, Q, depotIdx, info, distMat] = readGolden(filename, whichInstance);
else
    [coords, demand, Q, depotIdx, info, distMat] = readGoldenSingle(filename);
end

N = size(coords,1);
ord = [depotIdx, setdiff(1:N, depotIdx)];
coords = coords(ord,:);
demand = demand(ord);
if demand(1) ~= 0
    warning('loadGoldenInstance:nonzeroDepotDemand', ...
        ['Point 1 in "%s" has a non-zero raw demand (%.4g), but is being ' ...
         'treated as the depot (demand forced to 0) under this benchmark''s ' ...
         'convention that the depot is always listed/read first. If this ' ...
         'instance is unusual and point 1 is NOT actually the depot, this ' ...
         'will silently drop that demand -- inspect the raw file if this ' ...
         'instance''s results look wrong.'], filename, demand(1));
end
demand(1) = 0;  % depot demand is 0 by definition

NV = max(2, ceil(sum(demand)/Q));

prob = buildInstance(coords, demand, Q, NV, seed, catMode, distMat);
prob.name = info.Name;
prob.routeLengthLimit = info.RouteLengthLimit;  % informational only, not enforced by the solvers
end

% =========================================================================
function tf = isGoldenCombinedFormat(filename)
% Distinguish the two Golden-benchmark text formats by checking whether
% the file contains the combined format's telltale header phrase
% "with N Customers" anywhere; the single-instance format never contains
% this phrase (its header line is purely numeric).
fid = fopen(filename, 'r');
if fid == -1
    error('loadGoldenInstance:fileNotFound', 'File "%s" was not found or could not be opened.', filename);
end
raw = fread(fid, '*char')';
fclose(fid);
tf = ~isempty(regexp(raw, '\[\d+\]\s*with\s*\d+\s*Customers', 'once'));
end
