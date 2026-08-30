function d = dominates(a, b)
%DOMINATES  True if objective vector a dominates b (minimisation).
%   a dominates b iff a <= b in every component and a < b in at least one.
d = all(a <= b) && any(a < b);
end
