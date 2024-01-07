%**************************************************************
% MAIN_PROG - Solves the NK model with capital
%
% Code by Marco Brianti, University of Alberta
% January, 1 2022
%**************************************************************
cd '/Users/fangli/Dropbox/Marco_Fang_project/SimpleNKGPFA/NK model without capital'
addpath './funcs'

close all
clear
clc

disp('**********************************************************************')
disp('NK model without physical capital, Calvo pricing')
disp('Shocks: TFP, Government spending, monetary policy, inter-temporal preference shock')
disp('**********************************************************************')

%% parameters value
%%
% y c n w x1 x2 mc pii s m a t chi b int vp
subplottitles = {'Output, $y$', 'Consumption, $c$', 'Hours, $n$', 'Wage, $w$',...
                        '$x1$', '$x2$', 'Marginal Cost, $mc$', 'Inflation, $\pi$',...
                        '$s$', 'SDF', 'TFP, $a$', 'Goverment spending, $t$','Monetary policy, $\chi$',...
                        'Inter-temporal, $b$', 'Interest rate, $r$', 'price dispersion, $v^p$'};
shocktitles = {'Technology shock', 'Goverment spending shock', 'Monetary policy shock',...
    'Inter-temporal preference shock'};
plotting  = [1, 8, 2, 3, 4, 15, 11, 12, 13, 14];  % just for plotting
nss = size(shocktitles, 2);
horz = 41;


%% solving the model
%%
% Load Parameters
[flexp, setp] = paramset_Calvo;

% Compute state-space representation of model
[fyn, fxn, fypn, fxpn] = model_Calvo(flexp,setp);
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
% save './data/Calvo/IRFs_model.mat' IRFs_model