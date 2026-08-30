function prob = loadCVRPInstance(filename, catMode, seed)
%LOADCVRPINSTANCE  MOTDGVRP instance from a real CVRPLIB .vrp file.
%
%   prob = loadCVRPInstance(filename, catMode, seed)
%
%   Reads an actual Christofides Set E (or Augerat) .vrp file with readCVRP,
%   reorders nodes so the depot is row 1 (the convention used throughout this
%   package), and keeps the REAL demands and vehicle capacity. The
%   time-dependent green layer (3 slots, speed Table 3, GPH Table 4, load
%   surcharge Eq. 39) is then attached by buildInstance.
%
%   This reproduces the exact routing instances used by Nyako et al. (who
%   take E-n** from Christofides Set E and, for E-n51/76/101, multiply the
%   requests and capacity by 100). To match that scaling, pass a *-scaled
%   file or set prob.demand/prob.Q *100 after loading (see scaleInstance).
%
%   catMode : 'uniform' (paper default) or 'conflict'.
%   seed    : RNG seed for the category layout.

if nargin < 2 || isempty(catMode), catMode = 'uniform'; end
if nargin < 3 || isempty(seed),    seed = 7;            end

% Resolve bundled Set E instances when only a filename is supplied.
if ~isfile(filename)
    srcDir = fileparts(mfilename('fullpath'));
    bundledFile = fullfile(srcDir, '..', 'benchmarks', 'set_e', filename);
    if isfile(bundledFile)
        filename = bundledFile;
    end
end

[coords, demand, Q, depotIdx, info, distMat] = readCVRP(filename);

% Reorder so the depot is node 1 (the convention used throughout the package)
if ~isempty(coords)
    N = size(coords,1);
else
    N = size(distMat,1);
end
ord = [depotIdx, setdiff(1:N, depotIdx)];

if ~isempty(coords)
    coords = coords(ord,:);
end
if ~isempty(distMat)
    distMat = distMat(ord, ord);     % permute rows AND columns consistently
end
demand = demand(ord);
demand(1) = 0;                       % depot demand is 0 by definition

% Number of vehicles: parse from the name (…-k<NV>) if present, else estimate
NV = parseK(info.Name);
if isempty(NV) || NV < 1
    NV = max(2, ceil(sum(demand)/Q));
end

prob = buildInstance(coords, demand, Q, NV, seed, catMode, distMat);
prob.name = info.Name;
end

% -------------------------------------------------------------------------
function k = parseK(name)
% Extract the trailing -k<NV> from a CVRPLIB instance name.
k = [];
tok = regexp(name, '-k(\d+)', 'tokens', 'once');
if ~isempty(tok), k = str2double(tok{1}); end
end
