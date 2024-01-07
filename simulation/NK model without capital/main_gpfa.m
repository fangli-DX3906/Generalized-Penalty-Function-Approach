%**************************************************************
% MAIN_PROG - Solves the IRF using GPFA for the population economy
%
%**************************************************************
% cd '/Users/fangli/Dropbox/Marco_Fang_project/SimpleNKGPFA/NK model without capital'
cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/SimpleNKGPFA/NK model without capital'
addpath './funcs'
addpath './data'
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/bayesian GPFA/estimation');
% addpath './data/CRS'
% addpath './data/HF'
% addpath './data/Calvo'
% addpath './data/Nonsep'

clc
clear
close all

load IRFs_model.mat
load population.mat
load trueShockSeries.mat

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
% maxOrder = 4;   % original setting
maxOrder = 24;

% maximum lag order for selecting the aux VAR lag
maxAuxOrder = 24;  % original setting

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

%% gpfa
%%
% A) VAR order selection
ydata = Y;
ydata = ydata(:, var_idx);
% oo = orderSelect(ydata, maxOrder);

% B) VAR estimation
q = size(var_idx,2);
t = size(Y,1);
Yy = Y(:, var_idx);
Ss =  Y(:, aux_idx);
% p = oo;
% p = 12;
p = 4;

try
    [~, ~, ~, ~, ~, ~, mVAR] = estimVARCoeff(Yy, p, iscon, istr, var_names);
    [A,SIGMA] = olsvarc(Yy, p);	
%     [IRFxxx] = irfvar(A,SIGMA,p,h)
    if ~ any(abs(eig(A))>=1)
        [A] = asybc(A,SIGMA,t, p);
    end
    initDelta = 0;    
    [GAMMA_, delta, ~] = solver(mVAR, {{'Output', 'Inflation'}, {'Output', 'Inflation'}}, [[1,1]; [1,1]], [[1,1]; [1,-1]], initDelta, [1,2], 1, [], 1);
    IRFr = irfvar_PF(A,  SIGMA(1:q,1:q), p, h, GAMMA_);
catch
    disp(['In isim = ', num2str(isim), ' covariance matrix is not PD'])
end
xxx =  calcGPFARotation(SIGMA(1:q,1:q), 1-delta);   % exactly the same rotation

% calculate the shock correlation
shocks_gpfa = (inv(chol(SIGMA(1:q,1:q))'*GAMMA_)*(mVAR.auxMats.Residuals)')';
for i=1:nShocks
    temp = corrcoef(shocks_gpfa(:, i), trueShocks(p+1:end,i));
    corrl(i) = temp(1,2);
end
corrl

% C) auxillary VAR
[As, plist] = olsvaraux(Ss, Yy, auxcon, maxAuxOrder);
IRF_new = recoverIRF(nShocks, nVar, nAux, IRFr, As, plist, h, auxcon);
IRF_aux = IRF_new;
IRF_gpfa = vec(IRFr)';

% reshape
IRFs_aux = [];
for i = 1:size(aux_idx,2)*size(shock_names,2)
    IRFs_aux = [IRFs_aux; IRF_aux(1, (i-1)*horz+1:i*horz)];
end

% assembling IRF together
IRFs_gpfa = reshape(IRF_gpfa, nVar^2, h+1);
varD = 1: nVar;
varS = nVar+1:2*nVar;
auxD = 1:nAux;
auxS = nAux+1:2*nAux;
IRFs_population = [IRFs_gpfa(varD, :); IRFs_aux(auxD,:); IRFs_gpfa(varS, :); IRFs_aux(auxS,:)];

% save './data/IRFs_population.mat' IRFs_population 
% save './data/HF/IRFs_population.mat' IRFs_population 
% save './data/Calvo/IRFs_population.mat' IRFs_population 
% save './data/Nonsep/IRFs_population.mat' IRFs_population 
% save '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/IRFs_GPFA_population.mat' IRFs_population
sdsds

%% plotting
%%
IRFmat = {IRFs_population, IRFs_model};
legendString = {'GFPA IRF $\ \ \ $', 'Model Implied IRF $ \ \ \ $','Demand shock $ \ \ \ $';...
                        'GFPA IRF $\ \ \ $', 'Model Implied IRF $ \ \ \ $','Supply shock $ \ \ \ $'};
plots(40, {'b','r'}, dataset_names, shock_names, IRFmat, 0, 0, legendString);