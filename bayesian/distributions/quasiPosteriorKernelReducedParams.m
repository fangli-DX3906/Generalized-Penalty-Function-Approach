function quasi_posterior_kernel = quasiPosteriorKernelReducedParams(nu, phi, psi, omega, B, sigma)
    n = size(sigma,1);
    p = (size(B, 1)-1)/n;
    quasi_posterior_kernel = (det(psi))^(nu/2) * (det(omega))^(-(p+0.5)) * (det(sigma))^(-(nu+5)/2);
    quasi_posterior_kernel = quasi_posterior_kernel * exp(-trace(phi*inv(sigma))/2);
    xxx = vec(B-psi)' * inv(kron(sigma, omega)) * vec(B-psi);
    quasi_posterior_kernel = quasi_posterior_kernel * exp(-xxx/2);
end

