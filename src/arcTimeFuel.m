function [tij, FFij, dij] = arcTimeFuel(prob, i, j, departTime, load)
%ARCTIMEFUEL  Time-dependent travel time and fuel for one arc (i->j).
%
%   [tij, FFij, dij] = arcTimeFuel(prob, i, j, departTime, load)
%
%   Implements the FIFO time-dependent procedure of Section 4 of
%   Nyako et al. (2025), Eqs. (23)-(39). Node indices i,j are 1-based
%   (depot = 1). departTime is in hours from 08:00. load is the weight
%   carried on the arc (y^k_ij), used by the load surcharge Eq. (39).
%
%   The arc has a fixed total distance dij; the vehicle drives it at the
%   slot-dependent speed, possibly spanning two or three slots. We integrate
%   distance over time at each slot's speed until the whole arc is covered.

cat   = prob.arcCat(i,j);
speed = prob.speed(cat,:);      % 1x3 speeds (km/h) for this arc's category
gph   = prob.GPH(cat,:);        % 1x3 gallons-per-hour for this arc
zEnd  = prob.slotEnd;           % slot upper bounds in hours
dij   = prob.dist(i,j);

% Determine starting slot from departure time
s = currentSlot(departTime, zEnd);

remDist = dij;                  % distance still to cover
tNow    = departTime;
tij     = 0;                    % accumulated travel time on this arc
FEij    = 0;                    % accumulated base fuel (before load surcharge)

while remDist > 1e-12
    v  = speed(s);              % km/h in this slot
    gp = gph(s);                % gallons/hour in this slot

    % Time left before this slot ends
    if isinf(zEnd(s))
        tSlotLeft = inf;
    else
        tSlotLeft = zEnd(s) - tNow;
    end

    % Time needed to finish the arc at this speed
    tNeed = remDist / v;

    if tNeed <= tSlotLeft + 1e-12
        % Arc finishes within the current slot
        tij  = tij  + tNeed;
        FEij = FEij + tNeed * gp;
        remDist = 0;
        tNow = tNow + tNeed;
    else
        % Only part of the arc is covered in this slot; roll to next slot
        tij  = tij  + tSlotLeft;
        FEij = FEij + tSlotLeft * gp;
        remDist = remDist - v * tSlotLeft;
        tNow = zEnd(s);
        s = min(s + 1, prob.numSlots);   % advance slot (clamp at last)
    end
end

% Load-dependent fuel surcharge, Eq. (39):  FF = FE * (1 + (p/L)*load)
FFij = FEij * (1 + (prob.p/prob.L) * load);
end

% -------------------------------------------------------------------------
function s = currentSlot(t, zEnd)
% Return slot index (1..3) for time t given slot upper bounds zEnd.
s = find(t < zEnd, 1, 'first');
if isempty(s), s = numel(zEnd); end
end
