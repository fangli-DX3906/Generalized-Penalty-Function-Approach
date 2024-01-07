function A = multiNormalDist(psi, omega, size)

% size indicate the size of the random variables.
% For a VAR(p) with n variables and an intercept,
% size should be (n, 1+n*p).

num_of_elements = size(1)*size(2);
norm_rv = vec(psi) + chol(omega)'*randn(num_of_elements,1);
A = reshape(norm_rv, size(2), size(1));

end

