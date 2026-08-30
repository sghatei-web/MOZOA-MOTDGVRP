function prob = buildInstance(coords, demand, Q, NV, seed, catMode, distMat)
%BUILDINSTANCE  Construct a MOTDGVRP problem instance.
%
%   prob = buildInstance(coords, demand, Q, NV, seed, catMode, distMat)
%
%   catMode (optional): 'uniform' (default, 20% per category as in the paper)
%   or 'conflict' (length-biased categories that produce a richer Pareto
%   front for testing the multi-objective solvers).
%
%   distMat (optional): an N x N precomputed distance matrix. When supplied
%   (e.g. for EXPLICIT TSPLIB instances that give distances directly and have
%   no coordinates), it is used as-is and coords may be empty. When omitted,
%   distances are computed as Euclidean from coords.
%
%   Implements the time-dependent green VRP data model of
%   Nyako, Tayachi & Ben Abdelaziz (2025), "Machine learning multi-objective
%   optimization for time-dependent green vehicle routing problem",
%   Energy Economics 148, 108628.
%
%   INPUTS
%     coords : (n+1) x 2 matrix of node coordinates. Row 1 is the depot (node 0).
%     demand : (n+1) x 1 vector of demands. demand(1) (the depot) must be 0.
%     Q      : vehicle capacity (homogeneous fleet).
%     NV     : number of vehicles available.
%     seed   : RNG seed used when assigning a road category to each arc
%              (categories are uniformly distributed, 20% each, as in the paper).
%
%   OUTPUT  prob struct with fields used by every other routine.
%
%   The three time slots (Section 4.1 of the paper) are:
%       slot 1 : 08:00 -> 10:00   (z0 -> z1)
%       slot 2 : 10:00 -> 14:00   (z1 -> z2)
%       slot 3 : 14:00 -> beyond  (z2 -> inf)
%   All times are kept in HOURS measured from 08:00 (so z0 = 0).

if nargin < 5 || isempty(seed), seed = 42; end
if nargin < 6 || isempty(catMode), catMode = 'uniform'; end
if nargin < 7, distMat = []; end

demand = double(demand(:));

% Determine N (number of nodes) from coords or distMat
if ~isempty(coords)
    coords = double(coords);
    N = size(coords,1);
elseif ~isempty(distMat)
    N = size(distMat,1);
    coords = [];                       % no coordinates for EXPLICIT instances
else
    error('buildInstance:noGeometry', ...
          'Provide either coords or distMat.');
end
n = N - 1;                             % number of customers (excluding depot)

prob.n       = n;
prob.coords  = coords;
prob.demand  = demand;
prob.Q       = Q;
prob.NV      = NV;
prob.serviceTime = 10/60;             % 10 minutes, expressed in hours

% ---- Distance matrix: use supplied one, else Euclidean from coords --------
if ~isempty(distMat)
    D = double(distMat);
else
    D = zeros(N,N);
    for i = 1:N
        for j = 1:N
            D(i,j) = hypot(coords(i,1)-coords(j,1), coords(i,2)-coords(j,2));
        end
    end
end
prob.dist = D;

% ---- Time-slot boundaries (hours from 08:00) ------------------------------
% z(0)=0 (08:00), z(1)=2 (10:00), z(2)=6 (14:00)
prob.slotEnd = [2, 6, inf];           % upper bound of slot 1,2,3
prob.numSlots = 3;

% ---- Speed levels per category & slot  (Table 3, km/h) --------------------
%   rows = category 1..5, cols = period 1..3
prob.speed = [ 60 40 60 ;
               80 60 80 ;
               80 40 80 ;
               80 40 60 ;
               60 40 60 ];

% ---- Gallons-per-hour per category & slot (Table 4) -----------------------
prob.GPH = [ 3.692 3.125 6.417 ;
             4.923 4.688 8.556 ;
             4.923 3.125 8.556 ;
             4.923 3.125 6.417 ;
             3.692 3.125 8.556 ];

% ---- Load-dependent fuel surcharge (Eq. 39) -------------------------------
%   2% extra fuel per additional 100 load units  ->  p/L = 0.02/100
prob.p = 0.02;     % percentage increase
prob.L = 100;      % per L load units

% ---- Assign a road category (1..5) to every arc ---------------------------
% Two modes:
%   'uniform' : each category with 20% probability (the paper's default
%               simulation assumption).
%   'conflict': bias categories by arc length so that short arcs tend to be
%               slow/fuel-thirsty and long arcs efficient. This decouples
%               fuel/time from pure distance and yields a richer (multi-point)
%               Pareto front, which is useful for exercising the
%               multi-objective machinery.
if ~exist('catMode','var') || isempty(catMode), catMode = 'uniform'; end
rng(seed);
cat = ones(N,N);
dmax = max(D(:));
for i = 1:N
    for j = i+1:N
        if strcmpi(catMode,'conflict')
            r = D(i,j)/max(dmax,eps);
            if r < 1/3
                pool = [1 5 3];
            elseif r < 2/3
                pool = [4 3 1];
            else
                pool = [2 2 4];
            end
            c = pool(randi(numel(pool)));
        else
            c = randi(5);
        end
        cat(i,j) = c; cat(j,i) = c;
    end
end
prob.arcCat = cat;

end
