function  vc = vdvar(A,SIGMA,p,H,D)

% Calculate variance decomposition following the code by Kiliand and
% Lutkepol, Structural Vector AutoRegressive Anlysis - Table 4.1

% Input
% A is (Kp,Kp) and represents autoregressive matrix in companion form
% SIGMA is the variance-covariance matrix of reduced-form residual
% p is the number of lags
% H is the horizon up to whicht he variance is calculated
% D is an orthogonal matrix 

[K,~] = size(SIGMA);
vc    = zeros(K,K,H);

if nargin < 5
    D = eye(K);
end

for ih = 1:H
    h=ih;
    J=[eye(K,K) zeros(K,K*(p-1))];
    TH1=J*A^0*J'; TH=TH1*chol(SIGMA)'*D; TH=TH'; TH2=(TH.*TH); TH3=TH2;
    for i=2:h
        TH=J*A^(i-1)*J'*chol(SIGMA)'*D; TH=TH'; TH2=(TH.*TH); TH3=TH3+TH2;
    end
    TH4=sum(TH3);
    for j=1:K
        vc(j,:,ih)=TH3(j,:)./TH4;
    end
end
vc = permute(100*vc,[2 1 3]);
vc = reshape(vc,[],H);


end