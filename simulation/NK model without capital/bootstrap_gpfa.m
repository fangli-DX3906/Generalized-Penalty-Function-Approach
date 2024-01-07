%**************************************************************
% MAIN_PROG - bootstrap the CI for IRF using GPFA
%
%**************************************************************
%cd '/Users/fangli/Dropbox/Marco_Fang_project/SimpleNKGPFA/NK model without capital'
cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/SimpleNKGPFA/NK model without capital'
addpath './funcs'
addpath './data'
% addpath './data/CRS'
% addpath './data/HF'
% addpath './data/Calvo'
% addpath './data/Nonsep'

clc
clear
close all

load dlt.mat
load economies2000.mat
load shocksim2000.mat
load IRFs_model.mat
load IRFs_population.mat

%% parameters
%%
% baseline, CRS
totalVars = {'Output', 'Consumption', 'Hours', 'Wage', 'SDF', 'Multiplier', 'Inflation', 'TFP',...
                    'Goverment spending', 'Monetary policy shock', 'Discount factor', 'Interest rate'};
% HF
% totalVars = {'Output', 'Consumption', 'Hours', 'Wage', 'SDF', 'Multiplier', 'Inflation', 'Multiplier1', 'TFP',...
%                      'Goverment spending', 'Monetary policy shock', 'Discount factor', 'Interest rate', 'Comp Consumption'};
% Calvo
% totalVars = {'Output', 'Consumption', 'Hours', 'Wage', 'X1', 'X2', 'Marginal Cost', 'Inflation', 's',...
%                    'SDF', 'TFP', 'Goverment spending', 'Monetary policy shock', 'Discount factor',...
%                    'Interest rate', 'Price Dispersion'};
% Nonsep
% totalVars = {'Output', 'Consumption', 'Hours', 'Wage', 'SDF', 'Multiplier', 'Inflation', 'Multiplier1', 'TFP',...
%                      'Goverment spending', 'Monetary policy shock', 'Discount factor', 'Interest rate'};
                 
% baseline, CRS, HF, Nonsep
aux_names = {'Consumption', 'Hours', 'Wage', 'Interest rate', 'SDF', 'Multiplier'};
% Calvo
% aux_names = {'Consumption', 'Hours', 'Wage', 'Interest rate', 'SDF', 'Marginal Cost'};
               
%---------------------------------------------------%
shock_names = {'demand shock', 'supply shock'};
var_names  = {'Output', 'Inflation'};

% maximum lag order for selecting the VAR lag
% maxOrder = 8;   % original setting

% maximum lag order for selecting the aux VAR lag
% maxAuxOrder = 12;  % original setting

% length for IRF (including the zero period)
horz = 41;
h = horz-1;  

% intercept indicator for VAR
iscon = 1;

% intercept indicator for aux VAR
auxcon = 0;

% time trend for VAR
istr = 0;

% preprocessing
nsim = size(Ysim, 1);
dataset_names =  [var_names, aux_names];
var_idx = findingVars(var_names, totalVars);
aux_idx = findingVars(aux_names, totalVars);
nShocks = size(shock_names,2);
nVar = size(var_names,2);
nAux = size(aux_names,2);
nTotal = nVar + nAux;
nSystem = size(totalVars,2);
selectedVar = [var_idx, aux_idx, var_idx+nSystem, aux_idx+nSystem];
IRFs_model = IRFs_model(selectedVar,:);

%% gfpa
%%
optimList = [];
for isim=1:nsim
    ydata = Ysim{isim};
    ydata = ydata(:, var_idx);
    try
        oo = orderSelect(ydata, maxOrder);
    catch
        oo = -1;
    end
    optimList = [optimList, oo];
end

% B) estimation
q = size(var_idx,2);
t = size(Ysim{1},1);
IRF_gpfa = zeros(nsim,(q^2)*horz);
IRF_aux = zeros(nsim, nAux*nShocks*horz);

for isim=1:nsim  
    isim
    Y0 = Ysim{isim};
    Yy = Y0(:, var_idx);
    Ss =  Y0(:, aux_idx);
    p = optimList(isim);
    

    % PF
    try
        [~, ~, ~, ~, ~, ~, mVAR] = estimVARCoeff(Yy, p, iscon, istr, var_names);  % lag order selection
        [A,SIGMA] = olsvarc(Yy, p);
        if ~ any(abs(eig(A))>=1)
            [A] = asybc(A,SIGMA,t, p);
        end
%         initDelta = 0;
%         [GAMMA_, delta, ~] = solver(mVAR, {{'Output', 'Inflation'}, {'Output', 'Inflation'}}, [[1,1]; [1,1]], [[1,1]; [1,-1]], initDelta, [1,2], 1, [], 1);
%         delta
%         deltaList(isim) = delta;
        delt = 1- mVAR.auxMats.stdepsvec(2)/(mVAR.auxMats.stdepsvec(1)+mVAR.auxMats.stdepsvec(2));
        [GAMMA_, ~] = identifyPFA(mVAR, {{'Output', 'Inflation'}, {'Output', 'Inflation'}}, [[1,1]; [1,1]], [[1,1]; [1,-1]], delt, [1,2], 1);
        IRFr = irfvar_PF(A,  SIGMA(1:q,1:q), p, h, GAMMA_);
    catch
        disp(['In isim = ', num2str(isim), ' something is wrong with the identification!'])
        continue
    end
    
    % correlation
    shocks_gpfa = (inv(chol(SIGMA(1:q,1:q))'*GAMMA_)*(mVAR.auxMats.Residuals)')';
    for i=1:nShocks
        tempshock = Shocksim{isim};
        temp = corrcoef(shocks_gpfa(:, i), tempshock(p+1:end,i));
        corrl(i) = temp(1,2);
    end
    corrs(isim, :) = corrl;
    
    % aux reg
    [As, ps] = olsvaraux(Ss, Yy, auxcon, maxAuxOrder);
    IRF_temp = recoverIRF(nShocks, nVar, nAux, IRFr, As, ps, h, auxcon);
    IRF_aux(isim,:) = IRF_temp;
    IRF_gpfa(isim,:) = vec(IRFr)';   
end

% construct IRF
CI1_aux = prctile(IRF_aux,[2.5 97.5]);
% CI2_aux = prctile(IRF_aux, 50);
CILO1_aux = [];
CIUP1_aux = [];
% CImid_aux = [];
for i = 1:size(aux_idx,2)*size(shock_names,2)
    CILO1_aux = [CILO1_aux; CI1_aux(1, (i-1)*horz+1:i*horz)];
    CIUP1_aux = [CIUP1_aux; CI1_aux(2, (i-1)*horz+1:i*horz)];
%     CImid_aux = [CImid_aux; CI2_aux(1, (i-1)*horz+1:i*horz)];
end

CI1_gpfa   = prctile(IRF_gpfa,[2.5 97.5]);
% CI2_gpfa   = prctile(IRF_gpfa, 50);
CILO1_gpfa = reshape(CI1_gpfa(1,:),nVar^2,h+1);
CIUP1_gpfa = reshape(CI1_gpfa(2,:),nVar^2,h+1);
% CImid_gpfa = reshape(CI2_gpfa,nVar^2,h+1);

% assembling IRF together
varD = 1: nVar;
varS = nVar+1:2*nVar;
auxD = 1:nAux;
auxS = nAux+1:2*nAux;
CILO1 = [CILO1_gpfa(varD, :); CILO1_aux(auxD,:); CILO1_gpfa(varS, :); CILO1_aux(auxS,:)];
CIUP1 = [CIUP1_gpfa(varD, :); CIUP1_aux(auxD,:); CIUP1_gpfa(varS, :); CIUP1_aux(auxS,:)];
% CIMID = [CImid_gpfa(varD, :); CImid_aux(auxD,:); CImid_gpfa(varS,  :); CImid_aux(auxS,:)];

% save './data/HF/dlt.mat' deltaList
% save './data/HF/correlation2000.mat' corrs
% save './data/HF/CI2000.mat' CILO1 CIUP1

%% plotting
%%
IRFmat = {IRFs_population, IRFs_model};
CImat = {CILO1, CIUP1};
legendString = {'95\% Confidence Interval $ \ \ \ $', 'GFPA IRF $\ \ \ $', 'Model Implied IRF $ \ \ \ $', 'Demand shock $ \ \ \ $';...
                         '95\% Confidence Interval $ \ \ \ $', 'GFPA IRF $\ \ \ $', 'Model Implied IRF $ \ \ \ $', 'Supply shock $ \ \ \ $'};
plots(40, {'b','r'}, dataset_names, shock_names, IRFmat, 1, CImat, legendString);

%% dealing with correlation
%%
CI_corr = prctile(corrs,[2.5 97.5])
mid_corr = prctile(corrs, 50)