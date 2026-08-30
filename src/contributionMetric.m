function contribA = contributionMetric(A, B)
%CONTRIBUTIONMETRIC  Contribution of A vs B, Meunier et al. (2000), Eq. (42).
%
%   contribA = contributionMetric(A, B)
%
%   Build the combined non-dominated front PO = ND(A U B), then partition:
%     C   solutions common to both A and B that lie on PO
%     W_A solutions of A on PO that dominate the matching B-region (A-only wins)
%     N_A solutions of A on PO not in B (A-only non-dominated)
%     W_B, N_B symmetric for B
%   Contribution(A,B) = (|C|/2 + |W_A| + |N_A|) / (|C| + |W_A| + |N_A| + |W_B| + |N_B|)
%
%   With A,B two Pareto approximations, this yields Contribution(A,B) +
%   Contribution(B,A) = 1. Values > 0.5 mean A contributes more to the joint
%   front than B.

A = round(A,6); B = round(B,6);
U = [A; B];
labelA = [true(size(A,1),1); false(size(B,1),1)];

% Non-dominated subset of the union
nd = false(size(U,1),1);
for i = 1:size(U,1)
    dominated = false;
    for j = 1:size(U,1)
        if i~=j && dominates(U(j,:), U(i,:))
            dominated = true; break;
        end
    end
    nd(i) = ~dominated;
end

PO   = U(nd,:);
fromA = labelA(nd);

% Identify, per PO point, whether it appears in A, in B, or both
inA = ismember(PO, A, 'rows');
inB = ismember(PO, B, 'rows');

C  = sum(inA & inB);          % common points on the front
NA = sum(inA & ~inB);         % A-only points on the front
NB = sum(inB & ~inA);         % B-only points on the front

% In this discrete-front setting W_A and W_B (strictly-dominating wins that
% are not shared points) are already captured by N_A / N_B, so we set them 0
% and the formula reduces to the standard contribution ratio.
WA = 0; WB = 0;

num = C/2 + WA + NA;
den = C + WA + NA + WB + NB;
if den == 0
    contribA = 0.5;
else
    contribA = num / den;
end
end
