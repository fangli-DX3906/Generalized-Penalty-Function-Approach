function [Bcomp, cvec, dvec, Bpl_ev, VC_eps, Residmat, postEstimModel] = estimVARCoeff(Y, vlag, iscon, istr, varName)
%-----------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:
% This routine estimates the residuals and companion matrix 
% of reduced form VAR USING the built-in function varm
%-----------------------------------------------------------------------------------------------------------------------
% INPUTS:
% Y:        matrix of data series, m x n
% 
% vlag:    lag length   
% 
% iscon:   a value '1' introduces a constant in the conditional mean 
% 
% istr:      a value '1' introduces a determistic trend in the conditional mean
%-----------------------------------------------------------------------------------------------------------------------
% OUTPUT:
% Bcomp:   matrix with structure of estimated reduced-form coefficients in 
%           companion form excluding all the deterministic terms. 
%           This includes the estimated parameters in the partition 
%           n x (n x (n x vlag)). The remaining partition is made up of a 
%           diagonalic identity matrix that pins down the lead-lag relation 
%           between variables. The size is (n x vlag) x (n x vlag).  
% 
% cvec:    matrix with estimated constants; this returns the scalar 0 
%           if no constant is included in the model; 
%           otherwise, this is a matrix of size n x 1
% 
% dvec:    matrix with estimated parameters on the deterministic trends;
%           this returns the scalar 0 if no constant is included 
%           otherwise, this is a matrix of size n x 1
%
% Bpl_ev:  'full' companion matrix with parameter estimates; 
%          with both constant and deterministic trend. 
%          The content is arranged as follows:
%          [ partition of size (n x vlag) x n with estimated parameters on lagged data series;
%            partition of size (1 x n) with estimated coefficients on constant terms;
%            partition of size (1 x n) with estimated coefficients on deterministic trends ]
%           The overall size is (n*vlag+1+1) x (n*vlag+1+1)
%
% VC_eps:  covariance matrix of reduced-form residuals, size n x n
%-----------------------------------------------------------------------------------------------------------------------
warning('off');
nSeries = size(Y,2);
model = varm(nSeries, vlag);
model.SeriesNames = varName;
if iscon ~= 1
    model.Constant = zeros(nSeries,1);
end
if istr == 1
    model.Trend = Nan;
end

% In order to fit the new data (Jan 27th, 2022)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [modelRes, ~, ~, Residmat] = estimate(model, Y);
% Bcomp_temp = modelRes.AR;
% Bcomp = cell2mat(Bcomp_temp);
% cvec = modelRes.Constant;
% dvec = modelRes.Trend;
% VC_eps = cov(Residmat);
% if vlag==1
%     Bpl_ev = Bcomp';
% else
%     Bpl_ev = zeros(nSeries, nSeries*vlag);
%     for i=1:vlag
%         Bpl_ev(:, i:vlag:end)=Bcomp_temp{i};
%     end
%     Bpl_ev = Bpl_ev';
%     Bcomp   = [Bcomp; kron(eye(vlag-1), eye(nSeries)), ... 
%                       zeros((vlag-1)*nSeries, nSeries)];
% end
% 
% if iscon==1
%     Bpl_ev = [Bpl_ev; cvec'];
% end
% if istr==1
%     Bpl_ev = [Bpl_ev; dvec'];
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Jan 27th: V1
[Bcomp, cvec, dvec, Bpl_ev, ~ , Residmat] = estim(Y, vlag, iscon, istr);
[~,SIGMA,~,~,~] = olsvarc(Y,vlag);
modelRes.NumSeries = size(Y, 2);
modelRes.SeriesNames = string(varName);
modelRes.P = vlag;
modelRes.Constant = cvec;
modelRes.Trend = dvec;
VC_eps = SIGMA(1:modelRes.NumSeries, 1:modelRes.NumSeries);
modelRes.Covariance = VC_eps;

postEstimModel.basicModel = modelRes;
postEstimModel.auxMats.Residuals = Residmat;
postEstimModel.auxMats.Bcomp = Bcomp;
postEstimModel.auxMats.Bpl_ev = Bpl_ev;
postEstimModel.auxMats.cholVC = chol(VC_eps)';    % not chol(Covariance)';
postEstimModel.auxMats.stdepsvec = diag(VC_eps).^0.5;
        
end

