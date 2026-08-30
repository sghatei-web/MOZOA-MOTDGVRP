function c = cMetric(A, B)
%CMETRIC  Coverage metric C(A,B), Zitzler et al. (2003), Eq. (41).
%
%   c = cMetric(A, B)
%
%   A, B : objective matrices (rows = solutions, minimisation).
%   c    : fraction of B dominated by at least one solution in A, in [0,1].
%          Note C(A,B) is generally NOT 1 - C(B,A); compute both directions.

nB = size(B,1);
if nB == 0, c = 0; return; end
covered = 0;
for b = 1:nB
    for a = 1:size(A,1)
        if dominates(A(a,:), B(b,:)) || isequal(A(a,:), B(b,:))
            covered = covered + 1;
            break;
        end
    end
end
c = covered / nB;
end
