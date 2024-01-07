function J2 = calcJ2(n, H, p, A_0, covariance_mat, companion_mat)

sigma = covariance_mat;
IRF = calcReducedIRF(companion_mat, sigma, p, H);

J2 = zeros(n^2*(H+1), n^2*(H+1));
J2(1:n^2, 1:n^2) = eye(n^2);

for i=2:H+1
    J2((i-1)*n^2+1:i*n^2, 1:n^2) = kron(eye(n), reshape(IRF(:, i), n, []));
    J2((i-1)*n^2+1:i*n^2, (i-1)*n^2+1:i*n^2) = kron(inv(A_0), eye(n));
end

end

