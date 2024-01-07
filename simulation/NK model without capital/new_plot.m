%% adds path
%%
% new plot just to make comparison with the sign restriction approach
% date: 2022/10/7

cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/SimpleNKGPFA/NK model without capital'
addpath './funcs'
addpath './data'

clc
clear
close all

%% import data
%%
load IRFs_model.mat
IRFs_model = IRFs_model([1 7 13 19], :);
load economies2000.mat
load dlt.mat


%% parameters
%%
totalVars = {'Output', 'Consumption', 'Hours', 'Wage', 'SDF', 'Multiplier', 'Inflation', 'TFP',...
                    'Goverment spending', 'Monetary policy shock', 'Discount factor', 'Interest rate'};
aux_names = {'Consumption', 'Hours', 'Wage', 'Interest rate', 'SDF', 'Multiplier'};
shock_names = {'demand shock', 'supply shock'};
var_names  = {'Output', 'Inflation'};

% to check if this is better (yes it is)
maxOrder = 24;
maxAuxOrder =  24;

H = 41;
h = H-1;
iscon = 1;
auxcon = 0;
istr = 0;
nsim = size(Ysim, 1);
dataset_names =  [var_names, aux_names];
var_idx = findingVars(var_names, totalVars);
aux_idx = findingVars(aux_names, totalVars);
n_shocks = size(shock_names,2);
n_var = size(var_names,2);
n_aux = size(aux_names,2);

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
IRF_GPFA = zeros(q^2, H, nsim);
% IRF_gpfa = zeros(nsim,(q^2)*H);
% IRF_aux = zeros(nsim, n_aux*n_shocks*H);

for isim=1:nsim  
    isim
    Y0 = Ysim{isim};
    Yy = Y0(:, var_idx);
    Ss =  Y0(:, aux_idx);
    p = optimList(isim);

    [~, ~, ~, ~, ~, ~, mVAR] = estimVARCoeff(Yy, p, iscon, istr, var_names);  % lag order selection
    [A,SIGMA] = olsvarc(Yy, p);
    if ~ any(abs(eig(A))>=1)
        [A] = asybc(A,SIGMA,t, p);
    end
    delt = 1- mVAR.auxMats.stdepsvec(2)/(mVAR.auxMats.stdepsvec(1)+mVAR.auxMats.stdepsvec(2));
    [GAMMA_, ~] = identifyPFA(mVAR, {{'Output', 'Inflation'}, {'Output', 'Inflation'}}, [[1,1]; [1,1]], [[1,1]; [1,-1]], delt, [1,2], 1);
    IRFr = irfvar_PF(A,  SIGMA(1:q,1:q), p, h, GAMMA_);
    IRF_GPFA(:,:, isim) = IRFr;
    
    % aux reg
%     [As, ps] = olsvaraux(Ss, Yy, auxcon, maxAuxOrder);
%     IRF_temp = recoverIRF(n_shocks, n_var, n_aux, IRFr, As, ps, h, auxcon);
%     IRF_aux(isim,:) = IRF_temp;
%     IRF_gpfa(isim,:) = vec(IRFr)';   
end

save '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/IRF_GPFA_lag_24.mat' IRF_GPFA

%% construct IRF
%%
% CI1_aux = prctile(IRF_aux,[2.5 97.5]);
% CI2_aux = prctile(IRF_aux, 50);
% CILO1_aux = [];
% CIUP1_aux = [];
% CImid_aux = [];
% for i = 1:size(aux_idx,2)*size(shock_names,2)
%     CILO1_aux = [CILO1_aux; CI1_aux(1, (i-1)*H+1:i*H)];
%     CIUP1_aux = [CIUP1_aux; CI1_aux(2, (i-1)*H+1:i*H)];
%     CImid_aux = [CImid_aux; CI2_aux(1, (i-1)*H+1:i*H)];
% end
% 
% CI1_gpfa   = prctile(IRF_gpfa,[2.5 97.5]);
% CI2_gpfa   = prctile(IRF_gpfa, 50);
% CILO1_gpfa = reshape(CI1_gpfa(1,:), n_var^2, []);
% CIUP1_gpfa = reshape(CI1_gpfa(2,:), n_var^2, []);
% CImid_gpfa = reshape(CI2_gpfa, n_var^2, []);
% 
% % assembling IRF together
% varD = 1: n_var;
% varS = n_var+1:2*n_var;
% auxD = 1:n_aux;
% auxS = n_var+1:2*n_var;
% CILO1 = [CILO1_gpfa(varD, :); CILO1_aux(auxD,:); CILO1_gpfa(varS, :); CILO1_aux(auxS,:)];
% CIUP1 = [CIUP1_gpfa(varD, :); CIUP1_aux(auxD,:); CIUP1_gpfa(varS, :); CIUP1_aux(auxS,:)];
% CIMID = [CImid_gpfa(varD, :); CImid_aux(auxD,:); CImid_gpfa(varS,  :); CImid_aux(auxS,:)];
% 
% save '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/SimpleNKGPFA/NK model without capital/data/CI95_new_20221007.mat' CILO1 CIUP1 CIMID