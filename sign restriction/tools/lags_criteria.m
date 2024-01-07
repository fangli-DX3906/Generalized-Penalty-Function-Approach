function [AIC,BIC,HQ] = lags_criteria(sigma,N,T,max_lags)

%These criteria minimize the following objective functions:
%AIC(p) = T ln detOmega_hat + 2(n^2*p)
%BIC(p) = T ln detOmega_hat + (n^2*p)lnT
%HQ(p)  = T ln detOmega_hat + 2(n^2*p)lnlnT
%where n is the number of variables, T is the length of the dataset, and p
%is the number of lags

if nargin < 4
    max_lags = 8;
end

for nlags = 1:max_lags
      AIC(nlags)        = T*log(det(sigma)) + 2*(N^2*nlags);
      BIC(nlags)        = T*log(det(sigma)) + (N^2*nlags)*log(T);
      HQ(nlags)         = T*log(det(sigma)) + 2*(N^2*nlags)*log(log(T));
end

[AIC_minvalue, AIC]    = min(AIC);
[BIC_minvalue, BIC]    = min(BIC);
[HQ_minvalue, HQ]      = min(HQ);

end