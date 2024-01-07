function IRF = calcStructuralIRF(companion_mat, covariance_mat, p, H, rotation)

A = companion_mat;
sigma = covariance_mat;

n   = size(sigma,1);
J   = [eye(n,n) zeros(n,n*(p-1))];
IRF = reshape(J*A^0*J'*chol(sigma)'* rotation,n^2,1);

for i = 1:H
    IRF = ([IRF reshape(J*A^i*J'*chol(sigma)'*rotation,n^2,1)]);
end

end

