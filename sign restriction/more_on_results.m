%% set up
% This m-file plots the graphs that are presented on the first meeting 
%%
clear
close all
clc

% path
cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/code'
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1');
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/code/tools');  
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/Export_Fig');
addpath('/Users/fangli/Dropbox/Marco_Fang_project/SimpleNKGPFA/NK model without capital/data')

% import data
load IRFs_model.mat
IRFs_model = IRFs_model([1 7 13 19], :);
load IRF_simulated_bayesian.mat
load IRF_GPFA.mat
load IRFs_population.mat
% load IRF_GPFA_lag_24.mat
% load IRFs_population_lag_24.mat
IRFs_population = IRFs_population([1 2 9 10], :);


% some helper variables
totalVars = {'Output', 'Consumption', 'Hours', 'Wage', 'SDF', ... 
                   'Multiplier', 'Inflation', 'TFP','Goverment spending', ...
                   'Monetary policy shock', 'Discount factor', 'Interest rate'};
aux_names = {'Consumption', 'Hours', 'Wage', 'Interest rate', 'SDF', 'Multiplier'};   
var_names  = {'Output', 'Inflation'};
shock_names = {'demand shock', 'supply shock'};
yfnames = [var_names, aux_names];
var_idx = findingVars(var_names, totalVars);
aux_idx = findingVars(aux_names, totalVars);
neconomy = 2000;
maxOrder = 8;
maxAuxOrder =  12;
n_shock = length(shock_names);
n_var = length(var_names);
n_aux = length(aux_names);
n = n_var + n_aux;
H = 40;
median_pos = 10;
max_pos = 19;
min_pos = 1;
ci_95_max = 16;
ci_95_min = 4;
rows = [1 2 1+n 2+n];

%% camparison cross economies
%%
maxs = zeros(n_shock*n_var, H+1, neconomy);
mins = zeros(n_shock*n_var, H+1, neconomy);
medians = zeros(n_shock*n_var, H+1, neconomy);
averages = zeros(n_shock*n_var, H+1, neconomy);
ci_95_up = zeros(n_shock*n_var, H+1, neconomy);
ci_95_low = zeros(n_shock*n_var, H+1, neconomy);
CI_GPFA = quantile(IRF_small, [0 0.005 0.025 0.5 0.975 0.995 1], 3);
CI_GPFA(:,:,end+1) = mean(IRF_small, 3);

for i=1:neconomy
    IRF_temp = IRF_simulated{i};
    maxs(:,:,i) = IRF_temp(rows, :, max_pos);
    mins(:,:,i) = IRF_temp(rows, :, min_pos);
    medians(:,:,i) = IRF_temp(rows, :, median_pos);
    ci_95_up(:,:,i) = IRF_temp(rows, :, 16);
    ci_95_low(:,:,i) = IRF_temp(rows, :, 4);
end
CI_GPFA = quantile(IRF_small, [0 0.005 0.025 0.5 0.975 0.995 1], 3);
CI_GPFA(:,:,end+1) = mean(IRF_small, 3);

% keep
% CI95_max = quantile(ci_95_up, 0.5, 3);
% CI95_min = quantile(ci_95_low, 0.5, 3);
% CI95_max(:,:,end+1) = CI95_min;
% CI95_max(:,:,end+1) = IRFs_model;
% CI95_max(:,:,end+1) = CI_GPFA(:,:,3);
% CI95_max(:,:,end+1) = CI_GPFA(:,:,5);
% % CI_max(:,:,end+1) = mean(maxs, 3);
% % CI_max(:,:,end+1) = IRFs_model;
% titles = {'\textbf{Bayesian sign restriction on simulated economies: demand shock (max)}', '\textbf{Bayesian sign restriction on simulated economies: supply shock (max)}'};
% legend_string = {'GPFA 95\% Confidence Interval $ \ \ \ $', 'Sign 95\% Confidence Interval $ \ \ \ $', 'Model implied IRF $\ \ \ $';...
%                          'GPFA 95\% Confidence Interval $ \ \ \ $', 'Sign 95\% Confidence Interval $ \ \ \ $', 'Model implied IRF $\ \ \ $'};
% plot_cross_economy_CI(CI95_max, var_names, shock_names, {'b', 'r'}, H, titles, legend_string, false, '95 interval');
% pause(2)

%  h = 1:3;
CI95_min = quantile(ci_95_up, 0.5, 3);
CI95_max = quantile(ci_95_low, 0.5, 3);
CI_median = quantile(medians, [0 0.025 0.5 0.975 1], 3);
CI95(:,:,1) = CI95_min;                  % median of lower 95 bounds
CI95(:,:,end+1) = CI95_max;          % median of upper 95 bound
CI95(:,:,end+1) = CI_GPFA(:,:,3);   % upper 95 bounds of GPFA
CI95(:,:,end+1) = CI_GPFA(:,:,5);   % lower 95 bounds of GPFA
CI95(:,:,end+1) = IRFs_model;         % theoretical IRF
CI95(:,:,end+1) = CI_median(:,:,3);   % median of medians 
CI95(:,:,end+1) = IRFs_population;   % point estimates of GPFA
titles = {'\textbf{Demand shocks}', '\textbf{Supply shocks}'};
legend_string = {'sign restriction: 95\% CI  $ \ \ \ $', 'GPFA: 95\% CI $ \ \ \ $', 'theoretical IRF $\ \ \ $',  'sign restriction: median IRF $\ \ \ $',  'GPFA: point est IRF $\ \ \ $';...
                          'sign restriction: 95\% CI $ \ \ \ $', 'GPFA: 95\% CI $ \ \ \ $', 'theoretical IRF $\ \ \ $',  'sign restriction: median IRF $\ \ \ $',  'GPFA: point est IRF $\ \ \ $'};
h_to_plot = 40;
plot_cross_economy_CI(CI95, var_names, shock_names, {'b', 'r'}, h_to_plot, titles, legend_string, false, '95 interval');
pause(2)

% CI_median = quantile(medians, [0 0.025 0.5 0.975 1], 3);
% CI_median(:,:,end+1) = CI_GPFA(:,:,4);
% CI_median(:,:,end+1) = IRFs_model;
% titles = {' ',' '};
% legend_string = {'95\% Confidence Interval of median $ \ \ \ $', 'Sign Restriction Median IRF $\ \ \ $',  'GPFA Median IRF $\ \ \ $', 'Model implied IRF $\ \ \ $';...
%                          '95\% Confidence Interval of median $ \ \ \ $', 'Sign Restriction Median IRF $\ \ \ $', 'GPFA Median IRF $ \ \ \ $', 'Model implied IRF $\ \ \ $'};
% % legend_string = {' ' ; ' '};
% plot_cross_economy(CI_median, var_names, shock_names, {'b', 'r'}, H, titles, legend_string, true, 'median');
% pause(2)

% CI_max = quantile(maxs, [0 0.025 0.5 0.975 1], 3);
% CI_max(:,:,end+1) = mean(maxs, 3);
% CI_max(:,:,end+1) = IRFs_model;
% titles = {'\textbf{Bayesian sign restriction on simulated economies: demand shock (max)}', '\textbf{Bayesian sign restriction on simulated economies: supply shock (max)}'};
% legend_string = {'95\% Confidence Interval $ \ \ \ $', 'Median IRF $\ \ \ $', 'Mean IRF $ \ \ \ $';...
%                          '95\% Confidence Interval $ \ \ \ $', 'Median IRF $\ \ \ $', 'Mean IRF $ \ \ \ $'};
% plot_cross_economy(CI_max, var_names, shock_names, {'b', 'r'}, H, titles, legend_string, false, 'max');
% pause(2)
% 
% CI_min = quantile(mins, [0 0.025 0.5 0.975 1], 3);
% CI_min(:,:,end+1) = mean(mins, 3);
% titles = {'\textbf{Bayesian sign restriction on simulated economies: demand shock (min)}', '\textbf{Bayesian sign restriction on simulated economies: supply shock (min)}'};
% legend_string = {'95\% Confidence Interval $ \ \ \ $', 'Median IRF $\ \ \ $', 'Mean IRF $ \ \ \ $';...
%                          '95\% Confidence Interval $ \ \ \ $', 'Median IRF $\ \ \ $', 'Mean IRF $ \ \ \ $'};
% plot_cross_economy(CI_min, var_names, shock_names, {'b', 'r'}, H, titles, legend_string, false, 'min');
% pause(2)

%% GPFA IRF
%%
CI_GPFA = quantile(IRF_small, [0 0.005 0.025 0.5 0.975 0.995 1], 3);
CI_GPFA(:,:,end+1) = mean(IRF_small, 3);
save_fig  = false;
sizelabel = 20;
sizeticks = 18;
min_pos = 1;
max_pos = 7;
ci95_min_pos = 3;
ci95_max_pos = 5;
ci99_min_pos = 2;
ci99_max_pos = 6;
median_pos = 4;
mean_pos = 8;
color_scheme = {'b', 'r'};
fnames = {'\textbf{Generalized penalty function approach on simulated economies: demand shock}', '\textbf{Generalized penalty function approach on simulated economies: supply shock}'};
legend_string = {'99\% Confidence interval $ \ \ \ $', '95\% Confidence interval $ \ \ \ $', 'Median IRF $\ \ \ $', 'Mean IRF $ \ \ \ $', 'Model implied IRF $ \ \ \ $';...
                         '99\% Confidence interval $ \ \ \ $', '95\% Confidence interval $ \ \ \ $', 'Median IRF $\ \ \ $', 'Mean IRF $ \ \ \ $', 'Model implied IRF $ \ \ \ $'};

for ifig = 1:n_shock
    fig             = figure(ifig); 
    fig.Color       = 'w'; 
    fig.Position    = [-1681 40 1431 954]; 
    
    for in = 1: n_var
        subplot(1, n_var, in)
        plot_xtick = [0 5:5:H];   % x-axis ticks
        set(gca,'XTick',plot_xtick)
        grid on
        ic = in + n_var*(ifig-1);
        hold on
        x2 = [0:H, fliplr(0:H)];
        
        % The bigger CI
        inBetween = [CI_GPFA(ic,:, ci99_min_pos), fliplr(CI_GPFA(ic,:, ci99_max_pos))];
        hh1 = fill(x2, inBetween,color_scheme{ifig},'LineStyle','none');
        set(hh1,'facealpha',.05)

        % The smaller CI
        inBetween = [CI_GPFA(ic,:, ci95_min_pos), fliplr(CI_GPFA(ic,:, ci95_max_pos))];
        hh1 = fill(x2, inBetween,color_scheme{ifig},'LineStyle','none');
        set(hh1,'facealpha',.1)
        
        plot(0:H, CI_GPFA(ic,:, median_pos),'-','linewidth',2,'Color',color_scheme{ifig})
        plot(0:H, CI_GPFA(ic,:, mean_pos), '-o','linewidth',2,'Color',color_scheme{ifig})
        plot(0:H, IRFs_model(ic,:), '-o','linewidth',2,'Color','k')
        plot(0:H, zeros(1,H+1),'-','Color','k')
        
        xt = get(gca,'XTick');
        set(gca, 'FontSize', sizeticks)
        set(gca,'TickLabelInterpreter','latex')
        axis tight
        hold off
        clear new_yticks num_ticks a pct
        
        if ic == 1 + n_var*(ifig-1)
            xlabel('Quarter','fontsize',sizeticks,'interpreter','latex')
            ylabel('Percent','fontsize',sizeticks,'interpreter','latex')
        end
        title(var_names{in},'fontsize',sizelabel,'interpreter','latex')
    end
    sgtitle(fnames{ifig}, 'fontsize', sizelabel+7,'interpreter','latex')
    
    lgd = legend(legend_string{ifig,:});
    set(lgd,'Position',...
        [0.126041666666667 -0.00975710738586667 0.777788036318183 0.0564541461102402],...
        'Orientation','horizontal','FontSize',18,'interpreter','latex');
    pause(1)
    legend boxoff
    if save_fig
        cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1'
        export_fig(shock_names{ifig})
        cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/code'
    end
end
