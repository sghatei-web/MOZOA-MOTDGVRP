function v = igdPlus(F, R)
%IGDPLUS  Inverted Generational Distance plus (IGD+) of front F w.r.t. a
%   reference set R (minimisation convention throughout this project).
%
%   v = igdPlus(F, R)
%
%   F : K x m matrix of objective vectors produced by the algorithm being
%       scored (a Pareto front / archive).
%   R : L x m matrix of reference-set objective vectors. Since the true
%       Pareto front of MOTDGVRP is unknown, this project builds R as the
%       non-dominated subset of the POOLED objective vectors across every
%       algorithm and every run on a given instance (see
%       buildReferenceSet.m), i.e. the best-known approximation to the
%       true front, exactly mirroring the shared-reference-point
%       convention already used for hypervolume (Section 5 of the paper).
%
%   IGD+ (Ishibuchi et al., 2015) modifies the classical IGD distance so
%   that, for each reference point r in R, the "distance" to the
%   candidate front F only counts objectives where F is WORSE than r
%   (i.e. it does not penalise a candidate solution for being better than
%   the reference point on some objectives, which the original
%   Euclidean-distance IGD does). For minimisation:
%
%       d+(r, f) = sqrt( sum_k [ max(f_k - r_k, 0) ]^2 )
%       IGD+(F, R) = (1/|R|) * sum_{r in R} min_{f in F} d+(r, f)
%
%   Objectives are normalised to [0,1] using the combined range of F and
%   R (per-objective min/max) before computing distances, so that
%   distance/time/fuel -- which live on very different numeric scales --
%   contribute comparably; this mirrors how hypervolume is always
%   reported after a shared reference-point normalisation elsewhere in
%   this project. Lower IGD+ is better (closer to the reference set).
%
%   USAGE
%     v = igdPlus(mozoaFront, referenceSet);

if isempty(F) || isempty(R)
    v = Inf;
    return;
end

% ---- normalise both F and R to [0,1] per objective using their
%      combined range, so no single objective dominates the distance ----
combined = [F; R];
lo = min(combined, [], 1);
hi = max(combined, [], 1);
rng_ = hi - lo;
rng_(rng_ == 0) = 1;  % guard against a degenerate (constant) objective

Fn = (F - lo) ./ rng_;
Rn = (R - lo) ./ rng_;

L = size(Rn, 1);
dplus = zeros(L, 1);
for i = 1:L
    diffs = max(Fn - Rn(i, :), 0);         % only count "worse than r" gaps
    dists = sqrt(sum(diffs.^2, 2));
    dplus(i) = min(dists);
end
v = mean(dplus);
end
