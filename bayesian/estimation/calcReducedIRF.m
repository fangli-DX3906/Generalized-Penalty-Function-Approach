function IRF = calcReducedIRF(comp_mat,sigma, p, irf_length)
% IRFVAR.M
% Lutz Kilian
% University of Michigan
% April 1997

A = comp_mat;
n   = size(sigma,1);
J   = [eye(n,n) zeros(n,n*(p-1))];
IRF = reshape(J*A^0*J'*chol(sigma)',n^2,1);

for i = 1:irf_length
	IRF = ([IRF reshape(J*A^i*J'*chol(sigma)',n^2,1)]);
end

end

