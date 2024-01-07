function [IRF] = irfvar_PF(A,SIGMA, p, h, GAMMA)

n   = size(SIGMA,1);
J   = [eye(n,n) zeros(n,n*(p-1))];
IRF = reshape(J*A^0*J'*chol(SIGMA)'* GAMMA,n^2,1);
% IRF = reshape(J*A^0*J'*GAMMA,n^2,1);

for i = 1:h
    IRF = ([IRF reshape(J*A^i*J'*chol(SIGMA)'*GAMMA,n^2,1)]);
%      IRF = ([IRF reshape(J*A^i*J'*GAMMA,n^2,1)]);
end

