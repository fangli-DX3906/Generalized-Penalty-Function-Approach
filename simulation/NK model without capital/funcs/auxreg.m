function [A, SIGMAaux] = auxreg(s, y, p, iscon)

[t,qy]  = size(y);
y = y';
Y = y(:, p+1:end);
for i=1:p
	Y = [Y; y(:,p+1-i:t-i)];
end

if iscon
    X = [ones(1, size(Y, 2)); Y];
    divisor = t - (p+1)*qy - 1;
else
    X = Y;
    divisor = t - (p+1)*qy;
end

S = s';
S = S(:, p+1:end);
A = (S*X')/(X*X'); 
U = S - A*X;

SIGMAaux = U*U'/divisor;

end

