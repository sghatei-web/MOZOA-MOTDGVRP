function prob = scaleInstance(prob, factor)
%SCALEINSTANCE  Multiply demands and capacity by a factor (paper's x100).
%
%   prob = scaleInstance(prob, factor)
%
%   Nyako et al. multiply the requests and capacity of E-n51-k5, E-n76-k7,
%   and E-n101-k8 by 100. Apply this after loadCVRPInstance:
%       prob = loadCVRPInstance('E-n51-k5.vrp');
%       prob = scaleInstance(prob, 100);
%
%   Distances, speeds, and GPH are unaffected; only load-related quantities
%   (which drive the fuel surcharge Eq. 39 and the route split) are scaled.

if nargin < 2 || isempty(factor), factor = 100; end
prob.demand = prob.demand * factor;
prob.Q      = prob.Q      * factor;
end
