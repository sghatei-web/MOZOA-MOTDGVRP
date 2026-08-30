function idx = tournamentSelect(rank, cd, k)
%TOURNAMENTSELECT  Binary tournament on (rank, crowding distance).
%   idx = tournamentSelect(rank, cd, k) returns k selected indices.
%   Lower rank wins; ties broken by larger crowding distance.
P = numel(rank);
idx = zeros(1,k);
for t = 1:k
    a = randi(P); b = randi(P);
    if rank(a) < rank(b)
        idx(t) = a;
    elseif rank(b) < rank(a)
        idx(t) = b;
    elseif cd(a) >= cd(b)
        idx(t) = a;
    else
        idx(t) = b;
    end
end
end
