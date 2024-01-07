function [A, SIGMAaux] = olsvaraux(s, y, p)

global q

[t,qy]  = size(y);
qs = size(s,2);
y = y';
Y = y(:, p+1:end);
for i=1:p
	Y = [Y; y(:,p+1-i:t-i)];
end

rr = size(Y, 2);
X = [Y];

S = s';
S = S(:, p+1:end);
A = (S*X')/(X*X'); 
U     = S - A*X;
SIGMAaux = U*U'/(t-qs-p*(qs*qy)-1);

end
