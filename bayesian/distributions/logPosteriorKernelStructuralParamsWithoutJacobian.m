function quasi_posterior_kernel = logPosteriorKernelStructuralParamsWithoutJacobian(nu, phi, psi, omega, A_0, A_plus)

quasi_posterior_kernel =(nu+5)* log(abs(det(A_0)));      % according to Marco JMP: posterior nu=T-n;
% quasi_posterior_kernel =(nu+3)* log(abs(det(A_0)));     
xxx = vec(A_0)' * kron(eye(2), phi) * vec(A_0);
yyy = vec(A_plus - psi*A_0)' * inv(kron(eye(2), omega)) * vec(A_plus - psi*A_0);
quasi_posterior_kernel = quasi_posterior_kernel -xxx/2 -yyy/2;

end

