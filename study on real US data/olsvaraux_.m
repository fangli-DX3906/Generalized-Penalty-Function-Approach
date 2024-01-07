function [bet] = olsvaraux_(s, y, p)

global q

[t,q]  = size(y);
Y      = lagmatrix(y,0:p);
x       = [ones(t, 1), Y];
rnan  = ~any(isnan([x s]),2);
X      = x(rnan,:);
S       = s(rnan,:);
bet   = X'*X\X'*S;

end
