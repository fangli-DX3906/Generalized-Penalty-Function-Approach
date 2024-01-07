%****************************
% MAIN_PROG - GPFA Inference
% date: Dec 22, 2022
%****************************
cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/bayesian GPFA'
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/bayesian GPFA/data')
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/bayesian GPFA/distributions')
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/bayesian GPFA/estimation')

clc
clear
close all

%% import data and set global params
%%
load IRFs_model.mat
% load population.mat
% data = Y;
load economies2000.mat
data = Ysim{1}; 

total_names = {'Output', 'Consumption', 'Hours', 'Wage', 'SDF', 'Multiplier', 'Inflation', 'TFP',...
                     'Goverment spending', 'Monetary policy shock', 'Discount factor', 'Interest rate'};
var_names  = {'Output', 'Inflation'};
aux_names = {'Consumption', 'Hours', 'Wage', 'Interest rate', 'SDF', 'Multiplier'};
shock_names = {'demand shock', 'supply shock'};

var_idx = findVarsPos(var_names, total_names);
y = data(:, var_idx);

max_order = 12;
irf_length = 40;

p = 4;    % for now, set p=4
T = size(data,1);
n = length(var_names);

precision = 10^(-4);
num_burn = 2000; 
num_sim = 1000; 
num_tot = num_sim + num_burn;
pct = 0.95;

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
kernels = zeros(num_sim, 1);
for i=1:num_tot
    
    % draw from posterior     
    dsigma = inverseWishartDist(phi, nu); 
    dB = multiNormalDist(psi, kron(dsigma,omega), [n, n*p+1]);
    
    if i>num_burn
        disp([num2str(i - num_burn), ' simulation'])
        
        % check whether the solution exist       
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
        dJ =  calcJacobianStructuralwrtReduced(dB, dsigma, ddelta1, precision);
        
        % calulcate the posterior kernel
        a = logPosteriorKernelStructuralParamsWithoutJacobian(nu, phi, psi, omega, dA_0, dA_plus);
        pos_kernel = a - log(abs(det(dJ'*dJ)))/2;

        % calculate the jacobian
        dJ1 = calcJ1([], n, irf_length);
        dJ2 = calcJ2(n, irf_length, p, dA_0, dsigma, dcomp_mat);
        dJ3 = calcJ3(n, irf_length, p, dA_0, dA_plus, dsigma, dcomp_mat);
        dJ4 = calcJ4(n, irf_length, p, dA_0, dA_plus);
        
        % calculate the jacobian of IRF wrt (A_0, A_plus)
        dJJ = calcJacobianIRFwrtStructural(dJ1, dJ2, dJ3, dJ4);
        pos_kernel = pos_kernel - log(abs(det(dJJ'*dJJ)))/2;
        
        dirf = calcStructuralIRF(dcomp_mat, dsigma, p, irf_length, dQ);
        IRFs(:,:, i-num_burn) = dirf;
        VDs(:,:, i-num_burn) = calcStructuralVD(n, irf_length, dirf);
        kernels(i-num_burn) = pos_kernel;    
    end
end

%% plot
%%
num_pct = pct*num_sim;
[V, I] = sort(kernels, 'descend');

IRFs = IRFs(:,:, I);
VDs = VDs(:,:,I);
irf(:,:,1) = IRFs(:,:,1);   % mode
irf(:,:,2:3) = quantile(IRFs(:, :, 1:num_pct), [0 1], 3);   % 95% credible set
irf(:,:,4) = quantile(IRFs, 0.5, 3);  % median

vd(:,:,1) = VDs(:,:,1);
vd(:,:,2) = quantile(VDs, 0.5, 3);  % median

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
        
        inBetween = [irf(ic,1:H+1,2), fliplr(irf(ic,1:H+1,3))];
        hh1 = fill(x2, inBetween,Colrs,'LineStyle','none');
        set(hh1,'facealpha',.05)

        plot(0:H, irf(ic,1:H+1,1),'-','linewidth',3,'Color',Colrs)
        plot(0:H, irf(ic,1:H+1,4),'-.','linewidth',3,'Color',Colrs)

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
    lgd = legend(' 95\% confidence interval $ \ \ \ $', ' modal IRF $ \ \ \ $', ' median IRF $ \ \ \ $');
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
%         plot(0:H, vd(ic,1:H+1,1),'-','linewidth',3,'Color',Colrs)
%         plot(0:H, vd(ic,1:H+1,2),'-.','linewidth',3,'Color',Colrs)
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
%           'Orientation','horizontal','FontSize',20,'interpreter','latex');
%     pause(1)
%     legend boxoff
% end
