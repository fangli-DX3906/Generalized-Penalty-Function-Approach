function A= inverseWishartDist(phi, nu)

% Purpose:  Draws an m x m matrix from a inverse wishart distribution

wish = wishartDist(phi, nu);
A = inv(wish);

end

