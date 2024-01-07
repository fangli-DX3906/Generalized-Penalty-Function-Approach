% IRFVAR.M
% Lutz Kilian
% University of Michigan
% April 1997

function [IRF] = irfvar(A,SIGMA,p,h)

n   = size(SIGMA,1);
J   = [eye(n,n) zeros(n,n*(p-1))];
IRF = reshape(J*A^0*J'*chol(SIGMA)',n^2,1);

for i = 1:h
	IRF = ([IRF reshape(J*A^i*J'*chol(SIGMA)',n^2,1)]);
end
