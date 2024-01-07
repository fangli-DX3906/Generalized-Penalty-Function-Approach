function [uci, lci] = bootstrapCI(pemodel, numPath, alpha, irHoriz, precision)
% TODO: involking sovler.m func

rng(1);

if isempty(precision)
    precision = 10^(-5);
end
    
resid = pemodel.auxMats.Residuals;
info = pemodel.PFidentifyInfo;
nvar = pemodel.basicModel.NumSeries;
A_temp = pemodel.auxMats.cholVC;
onesMat = ones(size(resid));
CI = [];

% getting the bootstrap samples
[~, sampleIdx] = bootstrp(numPath, [], resid);

% MC simulation
for i = 1:numPath
    
    idx = sampleIdx(:,i);
    residual = resid(idx, :);
    
    mval = mean(residual);
    residual_centered = residual - onesMat .* mval;
    vc_eps = cov(residual_centered);

    vc_eps = cov(residual_centered);
    A = chol(vc_eps)';
    delta = 0;
    pemodel.auxMats.cholVC = A;
    disp(num2str(i))
    [GAMMA, ~] = solver(pemodel, info.var,  info.horiz, info.sign, delta, info.shockpos, info.isCC, precision, 0);
    A_mat = A * GAMMA;
    CI =[CI, {calcOrthoIRFs(pemodel, A_mat, 51, NaN)}];
    
end

uci = zeros(irHoriz,nvar,nvar);
lci = zeros(irHoriz,nvar,nvar);

for shock=1:nvar
    for var=1:nvar
        this = [];
        for n=1:numPath
            irf_this = CI{n};
            this = [this; irf_this(shock,:,var)];
        end
        upper = quantile(this, 1-(1-alpha)/2);
        lower = quantile(this, (1-alpha)/2);
        uci(:, var, shock) = upper;
        lci(:, var, shock) = lower;
    end
end

pemodel.auxMats.cholVC = A_temp;  % recover the original cholVC

end