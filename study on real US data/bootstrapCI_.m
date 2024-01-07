function [uci, lci] = bootstrapCI_(pemodel, numPath, alpha, irHoriz, precision)
% TODO: involking sovler.m func

rng(1);

if isempty(precision)
    precision = 10^(-5);
end
    
resid = pemodel.auxMats.Residuals;
info = pemodel.PFidentifyInfo;
nvar = pemodel.basicModel.NumSeries;
CI = [];

% getting the bootstrap sample
[~, sampleIdx] = bootstrp(numPath, [], resid);

% MC simulation
for i = 1:numPath
    
    idx = sampleIdx(:,i);
    residual = resid(idx, :);
%     mval = mean(residual);
%     residual_centered = residual - onesMat .* mval;
%     vc_eps = cov(residual_centered);
    vc_eps = cov(residual);
    A = chol(vc_eps)';
    
    objgam = 1;
    threshold = 0.1;   
    delta = 0; 
    add = 0.01;
    
    disp('Iterations Starts...')
    disp('***********************************************')
    while objgam >= precision
          while objgam >= threshold
                warning on
                [GAMMA, GAM] = BBBBB(pemodel, A, info.var,  info.horiz, info.sign, delta, info.shockpos, info.isCC);
                objgam = GAM(:,info.shockpos(1))' * GAM(:, info.shockpos(2));
                delta  = delta + add;
          end
          threshold = threshold/10;
          add = add/10;
          disp(['In this epoch, the target inner product is ', num2str(objgam)])
          disp(['In this epoch, delta is ', num2str(delta)])
          disp('***********************************************')
    end
    if abs(objgam) == objgam
          fprintf('\n')
          disp([num2str(i), '  Identification achieved!'])      
          fprintf('\n')
    else
          fprintf('\n')
          warning(['The target inner product is ' , num2str(objgam) ])
          fprintf('\n')
    end

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

end

