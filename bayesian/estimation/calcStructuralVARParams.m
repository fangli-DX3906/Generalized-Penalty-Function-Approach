function [A_0, A_plus] = calcStructuralVARParams(reduced_var_param, covariance_mat, delta_1)

% for a n-variable VAR(p) model with an intercept,
% reduced_var_param has the size of (n*p+1)*n
% covariance_mat has the size of n*n
% delta_1 is uniformly drawn from [0, 1]

B = reduced_var_param;
sigma = covariance_mat;
C = chol(sigma);   % C is a upper triangular Cholesky decomp, ie C' in the pdf

Q = calcGPFARotation(covariance_mat, delta_1);

A_0 = inv(C)*Q;
A_plus = B*inv(C)*Q;

end

