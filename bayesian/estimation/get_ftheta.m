function ftheta = get_ftheta(bet,bethat,X,p,sigma,sigmahat,U,T,SW,nuW)

% This function calculate the posterior probability of impulse response theta,
% f(theta). See pages 1 and 2 of Corrigendum to "Inference on impulse response 
% functions in structural VAR models" by Inoue and Kilian (2013, Journal of Econometrics)

% Inputs:
% bet: simulated value of the slopes of the reduced-form VAR
% bethat: OLS value of the slopes of the reduced-form VAR
% T: length of the data (time dimension)
% X: regressors in the reduced-form VAR
% p: number of lags in the reduced-form VAR
% sigma: simulated value of the variance of the reduced-form residuals
% sigmahat: OLS value of the variance of the reduced-form residuals
% U: orthogonal matrix U from the QR decomposition such that chol(sigma)*U is a candidate impact matrix
% SW: first argument of the inverse wishart
% nuW: second argument of the inverse wishart

[n,n2] = size(U);
if n ~= n2
    error('U must be a square matrix.')
end
A = chol(sigma)';
k = length(bet);

% Normalize the sign of U and A
W      = eye(n);
W(1,1) = -1;
detU = det(U);
if (detU - 1)^2 < eps
    Utilde = U;
    Atilde = A;
elseif (detU + 1)^2 < eps
    Utilde = W*U;
    Atilde = A*W;
else
    error('Determinant of U must be equal to one or minus one.')
end

% Derive S - Eq. 1
S = eye(n) - 2*(eye(n) + Utilde)^1;
obj1 = 1;
for in = 2:n
    obj1 = obj1*gamma(in/2)/(pi^(in/2));
end
obj2 = 2^((n-1)*(n-2)/2);
obj3 = (det(eye(n) + S))^(n-1);

% Derive fs - Eq. 2
fs = obj1*obj2/obj3;

% Derive f(B,sigma|data) - for diffuse prior, see among others Furlanetto
% et al. 2017 EJ
det(sigma);
fBsigma = det(sigma)^(-(T+n+1)/2) * exp(-1/2*(bet-bethat)' * kron(sigma^-1,X'*X) * (bet-bethat)) * ...
    exp(-1/2*trace(sigma^(-1)*sigmahat));

% Derive f(theta|data) - Eq. 4
ftheta = 2^(n*(n+1)/2) * (abs(det(Atilde))^(-n*p+1)) * obj3^(-1) * fBsigma * fs;

% % Posterio probability density function of bet given sigma
% fbet = (2*pi)^(-k/2) * det(sigma)^(-1/2) * exp(-1/2*(bet-bethat)'* kron(sigma^-1,X'*X) *(bet-bethat));
% 
% % Posterior probability density function of sigma
% fsigma = det(SW)^(nuW/2) / ( 2^(nuW*n/2) * gamma(nuW/2)) * det(sigma)^(-(nuW+n+1)/2) * exp(-1/2*trace(SW*sigma^(-1)));








end