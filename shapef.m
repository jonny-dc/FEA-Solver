function N = shapef(xi,eta)

% N = shapef(xi,eta)

% Given a coordinate in natural coordinates (xi,eta), output the value of
% the shape function matrix at that point:
% N = [N1 N2 N3 N4]

% INPUTS:
% - xi  : double
%         Value of xi at the requested point
% - eta : double
%         Value of eta at the requested point

% OUTPUTS:
% - N : double
%       Shape function matrix at requested point

N1 = (1-xi) .* (1-eta);  % Each individual shape function
N2 = (1+xi) .* (1-eta);
N3 = (1+xi) .* (1+eta);
N4 = (1-xi) .* (1+eta);

N = 0.25 .* [N1 N2 N3 N4];

end
