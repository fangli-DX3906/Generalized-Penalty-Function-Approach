function [AIC,BIC,HQ] = selectOptimLagOrder(y, max_lags)

n = size(y, 2);
T = size(y, 1);

%These criteria minimize the following objective functions:
%AIC(p) = T ln detOmega_hat + 2(n^2*p)
%BIC(p) = T ln detOmega_hat + (n^2*p)lnT
%HQ(p)  = T ln detOmega_hat + 2(n^2*p)lnlnT
%where n is the number of variables, T is the length of the dataset, and p
%is the number of lags

if nargin < 2
    max_lags = 8;
end

for nlags = 1:max_lags
      phim = (n^2*nlags + n);
      [~, ~, ~, ~, sigma] = calcReducedVARParams(y, nlags);
%       AIC(nlags)        = T*log(det(sigma)) + 2*(n^2*nlags);
%       BIC(nlags)        = T*log(det(sigma)) + (n^2*nlags)*log(T);
%       HQ(nlags)         = T*log(det(sigma)) + 2*(n^2*nlags)*log(log(T));
      AIC(nlags)        =  log(det(sigma)) + 2*phim/ T;
      BIC(nlags)        =  log(det(sigma)) + phim*log(T)/T;
      HQ(nlags)         =  log(det(sigma)) + 2*phim*log(log(T))/T;
end

[AIC_minvalue, AIC]    = min(AIC);
[BIC_minvalue, BIC]    = min(BIC);
[HQ_minvalue, HQ]      = min(HQ);

end

