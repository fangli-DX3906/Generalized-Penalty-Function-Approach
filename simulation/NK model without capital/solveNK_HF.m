%**************************************************************
% MAIN_PROG - Solves the NK model with capital
%
% Code by Marco Brianti, University of Alberta
% January, 1 2022
%**************************************************************
%cd '/Users/fangli/Dropbox/Marco_Fang_project/SimpleNKGPFA/NK model without capital'
cd 'C:\Users\brianti\Dropbox\Share_Marco_Fang\SimpleNKGPFA\NK model without capital'
addpath './funcs'

close all
clear
clc

disp('**********************************************************************')
disp('NK model without physical capital, external habit formation')
disp('Shocks: TFP, Government spending, monetary policy, inter-temporal preference shock')
disp('**********************************************************************')

%% parameters value
%%
% y c n w m lam pii cc xi a t chi betv r
subplottitles = {'Output, $y$',...     
                        'Consumption, $c$',...
                        'Hours, $n$',...
                        'Wage, $w$',...
                        'Stochastic discount factor, $m$',...
                        'Multiplier, $\lambda$',...
                        'Inflation, $\pi$',...
                        'Multiplier1, $\xi$',...
                        'Technology, $a$',...
                        'Goverment spending, $t$',...
                        'Monetary policy, $\chi$',...
                        'Discount factor, $\beta$',...
                        'Interest rate, $r$',...
                        'Comp Consumption, $CC$'};
shocktitles = {'Technology shock', 'Goverment spending shock', 'Monetary policy shock',...
    'Inter-temporal preference shock'};
plotting  = [1, 7, 2, 3, 4, 13, 9, 10, 11, 12];  % just for plotting
nss = size(shocktitles, 2);
horz = 41;


%% solving the model
%%
% Load Parameters
[flexp, setp] = paramset_HF;

% Compute state-space representation of model
[fyn, fxn, fypn, fxpn] = model_HF(flexp,setp);
[gx,hx]                = gx_hx_alt(fyn,fxn,fypn,fxpn);

%% model implied IRFs
%%
% Construct eta matrix
nx                  = length(hx);
ny                  = size(gx,1);
eta                 = zeros(nx,nss);
eta(1:nss,1:nss)    = eye(nss);

for is = 1:nss
    shockmat     = zeros(1,nss);
    shockmat(is) = 1;
    irfs         = ir(gx,hx,eta*shockmat',horz);
    if is==4
        demandShock = irfs';
%         demandShock = demandShock(total_idx, :);
    elseif is==1
        supplyShock = irfs';
%         supplyShock = supplyShock(total_idx, :);
    end
    figure
    for j=1:length(plotting)
        subplot(3,4,j)
        plot([0:horz-1],irfs(:,plotting(j)))
        title(subplottitles{plotting(j)},'interprete','latex','fontsize',16)
    end
    sgtitle(shocktitles{is},'interprete','latex','fontsize',24)
end

IRFs_model=[demandShock; supplyShock];
% save './data/HF/IRFs_model.mat' IRFs_model