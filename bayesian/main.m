%****************************
% MAIN_PROG - GPFA Inference
% date: Jan 27, 2023
%****************************
cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/bayesian GPFA'
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/bayesian GPFA/data')
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/bayesian GPFA/distributions')
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/bayesian GPFA/estimation')
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/code/tools')
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/bayesian GPFA/Export_Fig')

clc
clear
close all

%% import data and set global params
%%
load IRFs_model.mat
IRFs_model = IRFs_model([1 7 13 19], :);

load population.mat
data = Y;
% load economies2000.mat
% data = Ysim{101}; 

total_names = {'Output', 'Consumption', 'Hours', 'Wage', 'SDF', 'Multiplier', 'Inflation', 'TFP',...
                       'Goverment spending', 'Monetary policy shock', 'Discount factor', 'Interest rate'};
var_names  = {'Output', 'Inflation'};
aux_names = {'Consumption', 'Hours', 'Wage', 'Interest rate', 'SDF', 'Multiplier'};
shock_names = {'demand shock', 'supply shock'};

var_idx = findVarsPos(var_names, total_names);
y = data(:, var_idx);

aux_idx = findVarsPos(aux_names, total_names);
s =  data(:, aux_idx);

max_lag= 4; 
irf_length = 20;

[p, p_1, p_2] = selectOptimLagOrder(y, max_lag);
T = size(data,1);
n = length(var_names);
p = 8;

precision = 10^(-5);
num_burn = 1000; 
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
[B, ~, comp_mat, X, sigma] = calcReducedVARParams(y, p);

% set up the posterior params
phi = T*sigma; 
% nu = T+n;
nu = T-n;
% if p > 1
%     nu = (T-n)*(p-1);    
% else
%     nu = T-n;   
% end
omega = (X'*X)^(-1);
psi = B;

rng(3906);
IRFs_gpfa = zeros(n^2, irf_length+1, num_sim);
% VDs = zeros(n^2, irf_length+1, num_sim);
IRFs_sign = zeros(n^2, irf_length+1, num_sim);
kernels_gpfa = zeros(num_sim, 1);
% kernels_ = zeros(num_sim, 1);
kernels_sign = zeros(num_sim, 1);
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
        dJ_plus =  calcJacobianStructuralwrtReduced(dB, dsigma, ddelta1, precision);  % n^2+n(np+1)+1 by n^2+n(np+1)
        dJ_minus =  calcJacobianStructuralwrtReduced(dB, dsigma, ddelta1, -precision);  % n^2+n(np+1)+1 by n^2+n(np+1)
        dJ = (dJ_plus + dJ_minus)/2;
        
        % calulcate the posterior kernel
        a = logPosteriorKernelStructuralParamsWithoutJacobian(nu, phi, psi, omega, dA_0, dA_plus);
        pos_kernel = a - log(abs(det(dJ'*dJ)))/2;
%         kernels_(i-num_burn) = pos_kernel;    

        % calculate the jacobian
        dJ1 = calcJ1([], n, irf_length);
        dJ2 = calcJ2(n, irf_length, p, dA_0, dsigma, dcomp_mat);
        dJ3 = calcJ3(n, irf_length, p, dA_0, dA_plus, dsigma, dcomp_mat);
        dJ4 = calcJ4(n, irf_length, p, dA_0, dA_plus);
        
        % calculate the jacobian of IRF wrt (A_0, A_plus)
        dJJ = calcJacobianIRFwrtStructural(dJ1, dJ2, dJ3, dJ4);   % n^2(H+1) by n^2(min(H, p)+1)
        pos_kernel = pos_kernel - log(abs(det(dJJ'*dJJ)))/2;
        
        [irfsign, Qsign] = calcSignRestriction(dcomp_mat, dsigma, n, p, M, irf_length, signs, h1, h2);
        kernels_sign(i-num_burn) = get_ftheta(vec(dB), vec(psi), X, p, dsigma, sigma, Qsign, T, phi, nu);
        IRFs_sign(:,:, i-num_burn) = irfsign;
        
        dirf = calcStructuralIRF(dcomp_mat, dsigma, p, irf_length, dQ);
        IRFs_gpfa(:,:, i-num_burn) = dirf; 
%         VDs(:,:, i-num_burn) = calcStructuralVD(n, irf_length, dirf);
        kernels_gpfa(i-num_burn) = pos_kernel;    
    end
end

%% plot_1: median and 95% CI
%%
lower_bound = (1-pct)/2;
upper_bound = 1-(1-pct)/2;

IRFs_GPFA1(:,:,1:3) = quantile(IRFs_gpfa, [lower_bound 0.5 upper_bound], 3);   % 95% GPFA CI
IRFs_SIGN1(:,:,1:3) = quantile(IRFs_sign, [lower_bound 0.5 upper_bound], 3);   % 95% Sign CI

%% plot_2: modal and 95% credible sets
%%
num_pct = pct*num_sim;
[V_gpfa, I_gpfa] = sort(kernels_gpfa, 'descend');
[V_sign, I_sign] = sort(kernels_sign, 'descend');

IRFs_gpfa_sorted = IRFs_gpfa(:,:, I_gpfa);
IRFs_sign_sorted = IRFs_sign(:,:, I_sign);

IRFs_GPFA2(:,:,1) = IRFs_gpfa_sorted(:,:,1);   % modal IRF
IRFs_GPFA2(:,:,2:3) = quantile(IRFs_gpfa_sorted(:, :, 1:num_pct), [0 1], 3);   % 95% GPFA credible set
IRFs_SIGN2(:,:,1) = IRFs_sign_sorted(:,:,1);   % modal IRF
IRFs_SIGN2(:,:,2:3) = quantile(IRFs_sign_sorted(:, :, 1:num_pct), [0 1], 3);   % 95% GPFA credible set

plotIRFMedian;
close all
plotIRFModal;
close all

%% plot VDs
%%
% num_pct = pct*num_sim;
% [V, I] = sort(kernels, 'descend');
% VDs = VDs(:,:,I);
% VDs_GPFA(:,:,1) = VDs(:,:,1);
% VDs_GPFA(:,:,2) = quantile(VDs, 0.5, 3);  % median

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