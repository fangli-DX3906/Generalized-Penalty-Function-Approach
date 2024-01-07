function [B, sigma, A, X, sigma2] = calcReducedVARParams(y, p)
% Program estimates a level VAR with intercept in companion format by LS (by Luts Kilian)
% y is the data matrix with size t*n
% p is the lag order

[T, n] = size(y);
y = y';
Y = y(:,p:T);	

for i = 1:p-1
 	Y = [Y; y(:,p-i:T-i)];
end

X = [ones(1,T-p); Y(:,1:T-p)];
Y = Y(:,2:T-p+1);

A = (Y*X')/(X*X');   
U = Y - A*X;

B = A(1:n, :)';              % size: (n*p+1, n)

% check with Prof Brianti, which covariance matris should be the correct one.
sigma = U*U'/(T-p-p*n-1);	
sigma2 = U*U'/T;	   % for now, not adjust for DF
sigma = sigma(1:n, 1:n);
sigma2 = sigma2(1:n, 1:n);

A = A(:, 2:end);
X = X';

end

