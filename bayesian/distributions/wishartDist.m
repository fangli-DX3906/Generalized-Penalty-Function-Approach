function A = wishartDist(phi, nu)
% Command:  s=wish(h,n)
% Purpose:  Draws an m x m matrix from a wishart distribution
%                with scale matrix phi and degrees of freedom nu.
%                This procedure uses Bartlett's decomposition.
% Inputs:   phi   -- m x m scale matrix.
%               nu   -- scalar degrees of freedom.
% Outputs:  A   -- m x m matrix draw from the wishart distribution.
% Note: Parameterized so that mean is nu*phi
phi = inv(phi);
A = chol(phi)'*randn(size(phi,1),nu);
A = A*A';

end

