function dN = shapefgrad(xi,eta)

% dN = shapefgrad(xi,eta)

% Given a point in natural coordinates (xi,eta), return the gradient matrix
% of the shape functions at that point:
% dN = [[dN1/dxi  dN2/dxi  dN3/dxi  dN4/dxi ]
%       [dN1/deta dN2/deta dN3/deta dN4/deta]]

% INPUTS:
% xi  : double
%       Value of xi at the requested point
% eta : double
%       Value of eta at the requested point

% OUTPUTS:
% dN : double
%      Value of dN at the requested point

dN = zeros(2,4);

dN(1,1) = -(1-eta);
dN(1,2) = 1 - eta;
dN(1,3) = 1 + eta;
dN(1,4) = -(1+eta);
dN(2,1) = -(1-xi);
dN(2,2) = -(1+xi);
dN(2,3) = 1 + xi;
dN(2,4) = 1 - xi;

dN = 0.25 .* dN;

end
