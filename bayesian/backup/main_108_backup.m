%****************************
% MAIN_PROG - GPFA Inference
% date: Dec 29, 2022
%****************************
cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/bayesian GPFA'
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/bayesian GPFA/data')
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/bayesian GPFA/distributions')
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/bayesian GPFA/estimation')
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/code/tools')

% in Dropbox, all the needed mfiles are in the tools folder
% addpath('/Users/fangli/Dropbox/Marco_Fang_project/code/tools')
% addpath('/Users/fangli/Dropbox/Marco_Fang_project/code/data')

clc
clear
close all

%% import data and set global params
%%
load IRFs_model.mat
IRFs_model = IRFs_model([1 7 13 19], :);

% load population.mat
% data = Y;
load economies2000.mat
data = Ysim{626}; 

total_names = {'Output', 'Consumption', 'Hours', 'Wage', 'SDF', 'Multiplier', 'Inflation', 'TFP',...
                     'Goverment spending', 'Monetary policy shock', 'Discount factor', 'Interest rate'};
var_names  = {'Output', 'Inflation'};
aux_names = {'Consumption', 'Hours', 'Wage', 'Interest rate', 'SDF', 'Multiplier'};
shock_names = {'demand shock', 'supply shock'};

var_idx = findVarsPos(var_names, total_names);
y = data(:, var_idx);

aux_idx = findVarsPos(aux_names, total_names);
s =  data(:, aux_idx);

max_lag= 12;
irf_length = 40;

[p, p_1, p_2] = selectOptimLagOrder(y, max_lag);
T = size(data,1);
n = length(var_names);

precision = 10^(-5);
num_burn = 2000; 
num_sim = 1000; 
num_tot = num_sim + num_burn;
pct = 0.95;

% params used in sign restriction
signs = [1 1; 1, -1];
h1 = 1; 
h2 = 1;
M = 2; 

%% start with bayesian approach
%%
% get the reduced form params
[B, sigma, comp_mat, X] = calcReducedVARParams(y, p);

% set up the posterior params
phi = T*sigma; 
% nu = T-n;
if p > 1
    nu = (T-n)*(p-1);    
else
    nu = T-n;   
end
omega = (X'*X)^(-1);
psi = B;

rng(3906);
IRFs = zeros(n^2, irf_length+1, num_sim);
VDs = zeros(n^2, irf_length+1, num_sim);
IRFs_sign = zeros(n^2, irf_length+1, num_sim);
kernels = zeros(num_sim, 1);
kernels_ = zeros(num_sim, 1);
for i=1:num_tot
    
    % draw from posterior, just to burn the first 2000 draws to make the distribution stationary
    % actually there is only feedback from sigma to B, but no feedback from B to sigma
    % so actually these steps (line 71-72) seems to have no effect, but for now, just keep them
    % they are not time-consuming anyway.
    dsigma = inverseWishartDist(phi, nu); 
    dB = multiNormalDist(psi, kron(dsigma,omega), [n, n*p+1]);
    
    if i>num_burn
        disp([num2str(i - num_burn), ' simulation'])
        
        % check whether the solution exist, according to the derivation
        % if correlative is negative, then the delta1 should be in the
        % range: [-dsigma(1,2)/dsigma(2,2), -dsigma(1,1)/dsigma(1,2)]
        dtheta1 = 100;
        if dsigma(1,2)<0
            while dtheta1> -dsigma(1,1)/dsigma(1,2) || dtheta1 < -dsigma(1,2)/dsigma(2,2)
                ddelta1 = uniformDist;
                dtheta1 = ddelta1^(-1)-1;
            end
        else
            ddelta1 = uniformDist;
        end
        
        % calculate the rotation Q
        dQ = calcGPFARotation(dsigma, ddelta1);
        
        % recover the companion matrix from the draw
        dcomp_mat = makeCompMat(dB, n, p);
        
        % calculate A_0, A_+
        [dA_0, dA_plus] = calcStructuralVARParams(dB, dsigma, ddelta1);
        
        % calculate the jacobian of (A_0, A_plus) wrt (B, \Sigma, \delta_1)
        dJ =  calcJacobianStructuralwrtReduced(dB, dsigma, ddelta1, precision);  % n^2+n(np+1)+1 by n^2+n(np+1)
        
        % calulcate the posterior kernel
        a = logPosteriorKernelStructuralParamsWithoutJacobian(nu, phi, psi, omega, dA_0, dA_plus);
        pos_kernel = a - log(abs(det(dJ'*dJ)))/2;
        kernels_(i-num_burn) = pos_kernel;    

        % calculate the jacobian
        dJ1 = calcJ1([], n, irf_length);
        dJ2 = calcJ2(n, irf_length, p, dA_0, dsigma, dcomp_mat);
        dJ3 = calcJ3(n, irf_length, p, dA_0, dA_plus, dsigma, dcomp_mat);
        dJ4 = calcJ4(n, irf_length, p, dA_0, dA_plus);
        
        % calculate the jacobian of IRF wrt (A_0, A_plus)
        dJJ = calcJacobianIRFwrtStructural(dJ1, dJ2, dJ3, dJ4);   % n^2(H+1) by n^2(min(H, p)+1)
        pos_kernel = pos_kernel - log(abs(det(dJJ'*dJJ)))/2;
        
        IRFs_sign(:,:, i-num_burn) = calcSignRestriction(dcomp_mat, dsigma, n, p, M, irf_length, signs, h1, h2);
        dirf = calcStructuralIRF(dcomp_mat, dsigma, p, irf_length, dQ);
        IRFs(:,:, i-num_burn) = dirf; 
        VDs(:,:, i-num_burn) = calcStructuralVD(n, irf_length, dirf);
        kernels(i-num_burn) = pos_kernel;    
    end
end

%% sign restriction
%%
% fix the code, make it more readable. make the two methods using the same loop
% IRFxxx = bayesianSignRestriction(y, s, max_lag, 12, num_burn, num_sim, M, irf_length+1, signs, h1, h2);
% IRFs_sign = quantile(IRFxxx, [(1-pct)/2 0.5 1-(1-pct)/2], 3);
% IRFs_sign = IRFs_sign([1 2 9 10], :, :);

%% plot
%%
num_pct = pct*num_sim;
[V, I] = sort(kernels, 'descend');

IRFs = IRFs(:,:, I);
VDs = VDs(:,:,I);
IRFs_GPFA(:,:,1) = IRFs(:,:,1);   % mode
IRFs_GPFA(:,:,2:3) = quantile(IRFs(:, :, 1:num_pct), [0 1], 3);   % 95% credible set
IRFs_GPFA(:,:,4) = quantile(IRFs, 0.5, 3);  % median

VDs_GPFA(:,:,1) = VDs(:,:,1);
VDs_GPFA(:,:,2) = quantile(VDs, 0.5, 3);  % median

IRFs_sign = quantile(IRFs_sign, [(1-pct)/2 0.5 1-(1-pct)/2], 3);

% plot IRFs
sizelabel = 24;
sizeticks = 18;
H = irf_length;
for ifig = 1:size(shock_names,2)
    fig = figure(ifig); 
    fig.Color = 'w'; 
    fig.Position = [-1681 40 1431 954];
    
    if ifig == 1
        Colrs = [0 0 1];
    else
        Colrs = [1 0 0];
    end
    
    for in = 1: size(var_names,2)
        subplot(ceil(size(var_names,2)/2), 2,in)
        plot_xtick = [0 5:5:H];   
        set(gca,'XTick',plot_xtick)
        grid on

        ic = in + size(var_names,2)*(ifig-1);
        hold on
        x2 = [0:H, fliplr(0:H)];
        
        inBetween = [IRFs_GPFA(ic,1:H+1,2), fliplr(IRFs_GPFA(ic,1:H+1,3))];
        hh1 = fill(x2, inBetween,Colrs,'LineStyle','none');
        set(hh1,'facealpha', .25)
        
        inBetween = [IRFs_sign(ic,1:H+1,1), fliplr(IRFs_sign(ic,1:H+1,3))];
        hh2 = fill(x2, inBetween,Colrs,'LineStyle','none');
        set(hh2,'facealpha', .15)

        plot(0:H, IRFs_GPFA(ic,1:H+1,1),'-','linewidth',3,'Color',Colrs)   % modal IRF
%         plot(0:H, IRFs_GPFA(ic,1:H+1,4),'-.','linewidth',3,'Color',Colrs)
        plot(0:H, IRFs_model(ic,1:H+1),'-o','linewidth',3,'Color',Colrs)   % model-implied IRF
        plot(0:H, IRFs_sign(ic,1:H+1,2),'-.','linewidth',3,'Color',Colrs)    % sign restriction IRF

        xt = get(gca,'XTick');
        set(gca, 'FontSize', sizeticks)
        set(gca,'TickLabelInterpreter','latex')
        axis tight
        hold off
        clear new_yticks num_ticks a pct
        if ic == 1 + size(var_names,2)*(ifig-1)
            xlabel('Quarter','fontsize',sizeticks,'interpreter','latex')
            ylabel('Percent','fontsize',sizeticks,'interpreter','latex')
        end
        title(var_names{in},'fontsize',sizelabel,'interpreter','latex') % 24
    end
    sgtitle(shock_names{ifig}, 'fontsize', sizelabel+7,'interpreter','latex')
    lgd = legend(' 95\% GPFA credible set $ \ \ \ $', ' 95\% SR credible set $ \ \ \ $', ' modal GPFA IRF $ \ \ \ $', ' model-implied IRF $ \ \ \ $', ' sign restriction IRF $ \ \ \ $');
    set(lgd,'Position', [0.126041666666667 -0.00975710738586667 0.777788036318183 0.0564541461102402],...
         'Orientation','horizontal','FontSize',20,'interpreter','latex');
    pause(1)
    legend boxoff
end

% plot VDs
% close all
% sizelabel = 24;
% sizeticks = 18;
% H = irf_length;
% for ifig = 1:size(shock_names,2)
%     fig = figure(ifig); 
%     fig.Color = 'w'; 
%     fig.Position = [-1681 40 1431 954];
%     
%     if ifig == 1
%         Colrs = [0 0 1];
%     else
%         Colrs = [1 0 0];
%     end
%     
%     for in = 1: size(var_names,2)
%         subplot(ceil(size(var_names,2)/2), 2,in)
%         plot_xtick = [0 5:5:H];   
%         set(gca,'XTick',plot_xtick)
%         grid on
% 
%         ic = in + size(var_names,2)*(ifig-1);
%         hold on
%         x2 = [0:H, fliplr(0:H)];
%         
%         plot(0:H, VDs_GPFA(ic,1:H+1,1),'-','linewidth',3,'Color',Colrs)
%         plot(0:H, VDs_GPFA(ic,1:H+1,2),'-.','linewidth',3,'Color',Colrs)
% 
%         xt = get(gca,'XTick');
%         set(gca, 'FontSize', sizeticks)
%         set(gca,'TickLabelInterpreter','latex')
%         axis tight
%         hold off
%         clear new_yticks num_ticks a pct
%         if ic == 1 + size(var_names,2)*(ifig-1)
%             xlabel('Quarter','fontsize',sizeticks,'interpreter','latex')
%             ylabel('Percent','fontsize',sizeticks,'interpreter','latex')
%         end
%         title(var_names{in},'fontsize',sizelabel,'interpreter','latex') % 24
%     end
%     sgtitle(shock_names{ifig}, 'fontsize', sizelabel+7,'interpreter','latex')
%     lgd = legend(' modal VD $ \ \ \ $', ' median VD $ \ \ \ $');
%     set(lgd,'Position', [0.126041666666667 -0.00975710738586667 0.777788036318183 0.0564541461102402],...
%          'Orientation','horizontal','FontSize',20,'interpreter','latex');
%     pause(1)
%     legend boxoff
% end
