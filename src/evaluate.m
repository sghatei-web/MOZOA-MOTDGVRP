function [obj, routes, penalty] = evaluate(chrom, prob)
%EVALUATE  Three objectives of the MOTDGVRP for one chromosome.
%
%   [obj, routes, penalty] = evaluate(chrom, prob)
%
%   obj = [totalDistance, totalTime, totalFuel]   (Eqs. 1-3 of the paper).
%
%   Routing logic per vehicle:
%     - vehicle leaves the depot at z(0) = 08:00 (t = 0) carrying the total
%       load of its route (Eq. 7-8 load conservation: load decreases as
%       customers are served);
%     - travel time and fuel on each arc are time-dependent (arcTimeFuel);
%     - service time st_i is added at each customer (Eq. 13);
%     - the load carried on arc (i,j) is the demand still to be delivered to
%       the customers from j onward, which drives the fuel surcharge Eq. 39.
%
%   penalty counts vehicles used beyond prob.NV (infeasible fleet size); the
%   driver/algorithm may use it to discourage infeasible splits.

routes = decodeRoutes(chrom, prob);

totalDist = 0; totalTime = 0; totalFuel = 0;
depot = 1;

for r = 1:numel(routes)
    route = routes{r};
    seq = [depot, route, depot];          % depot ... depot

    % Load carried when leaving depot = sum of demands on this route
    remainingLoad = sum(prob.demand(route));
    t = 0;                                % depart depot at 08:00

    for k = 1:numel(seq)-1
        i = seq(k); j = seq(k+1);

        [tij, FFij, dij] = arcTimeFuel(prob, i, j, t, remainingLoad);

        totalDist = totalDist + dij;
        totalTime = totalTime + tij;
        totalFuel = totalFuel + FFij;

        t = t + tij;
        % After arriving at j, serve it (if it is a customer) and drop its load
        if j ~= depot
            t = t + prob.serviceTime;
            remainingLoad = remainingLoad - prob.demand(j);
        end
    end
end

obj = [totalDist, totalTime, totalFuel];
penalty = max(0, numel(routes) - prob.NV);
end
