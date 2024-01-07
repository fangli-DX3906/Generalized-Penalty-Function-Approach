%% Sign-restriction identified VAR with Bayesian estimation using natural conjugate Gaussian
%%
% picking the IRF according to Kilian and Inoue (2020)

%% warm ups
%%
clear
close all
clc

cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1'
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/code');
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/code/tools');   
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/SimpleNKGPFA/NK model without capital/data');

% Data read
load IRF_bayesian_3906.mat
load economies2000.mat
load IRFs_model.mat
% IRFs_model = IRFs_model([1 7 2 3 4 12 5 6 13 19 14 15 16 24 17 18], :);
IRFs_model = IRFs_model([1 7 13 19], :);
load IRFs_GPFA_population.mat
% IRFs_population = IRFs_population([1 7 2 3 4 12 5 6 13 19 14 15 16 24 17 18], :);
IRFs_population = IRFs_population([1 2 9 10], :);

% Endogenous variables in the system
totalVars = {'Output', 'Consumption', 'Hours', 'Wage', 'SDF', 'Multiplier',... 
                   'Inflation', 'TFP','Goverment spending', ...
                   'Monetary policy shock', 'Discount factor', 'Interest rate'};   % order for IRFs_model
aux_names = {'Consumption', 'Hours', 'Wage', 'Interest rate', 'SDF', 'Multiplier'};    % order for IRFs_population
var_names  = {'Output', 'Inflation'};
shock_names = {'demand shock', 'supply shock'};
yfnames = [var_names, aux_names];
var_idx = findingVars(var_names, totalVars);
aux_idx = findingVars(aux_names, totalVars);
neconomy = size(Ysim, 1);
maxOrder = 8;
maxAuxOrder =  12;
n_shock = length(shock_names);
N = length(var_idx);
T = size(Ysim{1},1);

%% Reduced-form VAR with frequentist approach and checking
%%
tic
for isim=1:neconomy
    ydata = Ysim{isim};
    y = ydata(:, var_idx);
    p = orderSelect(y, maxOrder);
    [A,~,Uhat,const,X,SIGMA] = olsvarc(y,p);	 % Frequentist approach: VAR with intercept
    [AIC,BIC,HQ]             = lags_criteria(SIGMA(1:N,1:N),N,T);

    % Check if the matrix algebra is correct
    Bhat        = [const(1:N,1), A(1:N,:)]'; % re-add constant and remove identity matrix from the companion form
    sigmautilde = SIGMA(1:N,1:N);
    Y           = y(1+p:end,:);
    X           = X';
    checkmat    = Y - X*Bhat - Uhat(1:N,:)';   % check if the matrix combo is correct, should be teensy-weensy if correct!!

    % Check if the vectorization algebra is correct
    ybar     = reshape(y(1+p:end,:),[],1);
    ubar     = reshape(Uhat(1:N,:)',[],1);
    bethat   = reshape(Bhat,[],1);
    checkvec = ybar - kron(eye(N),X)*bethat - ubar;

    if sum(sum(checkmat.^2)) > 10^(-10) && sum(checkvec.^2) > 10^(-10)
        error('Matrix algebra or vectorization algebra is not correct')
        sound(sin(2*pi*25*(1:4000)/100));
    else
        disp(['Simulation', num2str(isim), ' : Everything seems good.'])
    end 
end
toc

%% Bayesian approach
%%
rng(3906);
nburn = 2000; 
nsim = 2500; 
M = 2; % for each simulation, required number of accepted rotations
H = 41; % horizon IRFs
signs = [1 1; 1, -1];
h1 = 1; % first horizon to implement the restriction
h2 = 1; % last horizon to implement the restriction. Cumulate from h1 to h2
IRF_bk = cell(neconomy, 4);

% p = parpool(4);
% parfor n = 1:2000
for n = 1:2000
    disp(['---> Estimating ', num2str(n), ' economy <---' ]);
    ydata = Ysim{n};
    y = ydata(:, var_idx);
    s = ydata(:, aux_idx);
    % first QL sets, second QL opt; third and fourth AL; fifth and sixth angular
    [IRF_bk{n, 1},  IRF_bk{n, 2},  IRF_bk{n, 3},  IRF_bk{n, 4}]= bayesianSignRestrictionKilian(y, s, maxOrder, maxAuxOrder, nburn, nsim, M, H, signs, h1, h2);
end
% delete(p);

% save '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/IRF_bayesian_kilian.mat' IRF_bk

%% Plot
%%
pick = 1128;
% tic
% Y0 = Ysim{pick};
% Yy = Y0(:, var_idx);
% Ss =  Y0(:, aux_idx);
% [A, B, C, D, E, F] = bayesianSignRestrictionKilian(Yy, Ss, maxOrder, maxAuxOrder, nburn, nsim, M, H, signs, h1, h2);
% toc
% save '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/econ1128_KI.mat' A B C D E F
load econ1128_KI.mat
CI_q = quantile(A, [0, 1], 3);
CI_a = quantile(C, [0, 1], 3);
CI_c = quantile(E, [0, 1], 3);
IRF_conventional = IRF_bayesian{pick};
IRF_conventional = IRF_conventional([1 2 9 10],:,:);
legend_string = {'95\% K\&I CI $ \ \ \ $',  '95\% BS CI $\ \ \ $', 'BS median IRF $\ \ \ $', 'Opt K\&I IRF $ \ \ \ $',  'GPFA population IRF $ \ \ \ $', 'Model implied IRF $ \ \ \ $';...
                          '95\%  K\&I CI $ \ \ \ $',  '95\% BS CI $\ \ \ $', 'BS median IRF $\ \ \ $', 'Opt K\&I IRF $ \ \ \ $', 'GPFA population IRF $ \ \ \ $', 'Model implied IRF $ \ \ \ $'};
opt = F;
band = CI_c;

% Figure impulse response functions: quadratic
sizelabel = 24;
sizeticks = 18;
for ifig = 1:n_shock
    fig             = figure(ifig); % name figure
    fig.Color       = 'w'; % white backgroung
    %fig.WindowState = 'maximized'; % maximize figure size
    fig.Position    = [-1681 40 1431 954]; % 
    if ifig == 1
        Colrs = [0 0 1];
    else Colrs = [1 0 0];
    end
    for in = 1: 2
        subplot(1,2,in)
        plot_xtick = [0 5:5:H-1];   % x-axis ticks
        set(gca,'XTick',plot_xtick)
        grid on
        ic = in + 2*(ifig-1);
        hold on
        x2 = [0:H-1, fliplr(0:H-1)];

        % credible set
        inBetween = [band(ic,1:H,1), fliplr(band(ic,1:H,2))];
        hh1 = fill(x2, inBetween,Colrs,'LineStyle','none');
        set(hh1,'facealpha',.1)
        % conventional confidence bands for bayesian
        inBetween = [IRF_conventional(ic,1:H,4), fliplr(IRF_conventional(ic,1:H,16))];
        hh3 = fill(x2, inBetween,Colrs,'LineStyle','none');
        set(hh3,'facealpha',.25)
        % convengtional median
        plot(0:H-1,IRF_conventional(ic,1:H, 10),'--d','linewidth',2,'Color',Colrs)
        % GPFA IRF
%         plot(0:H-1,IRF_GPFA(ic,1:H),'--+','linewidth',2.5,'Color',Colrs)
        % optimal IRF
        plot(0:H-1,opt(ic,1:H),'-','linewidth',1.5,'Color',Colrs)
        % GPFA population IRF
        plot(0:H-1,IRFs_population(ic,1:H), '--+','linewidth',1.5,'Color','k')
        % model implied IRF
        plot(0:H-1,IRFs_model(ic,1:H), '--o','linewidth',1.5,'Color','k')
        plot(0:H-1,zeros(1,H),'-','Color','k')
        xt = get(gca,'XTick');
        set(gca, 'FontSize', sizeticks)
        %title([varname],'fontsize',20,'interpreter','latex') % 24
        set(gca,'TickLabelInterpreter','latex')
        axis tight
        hold off
        %grid on
        clear new_yticks num_ticks a pct
        if ic == 1 + size(yfnames,2)*(ifig-1)
            xlabel('Quarter','fontsize',sizeticks,'interpreter','latex')
            ylabel('Percent','fontsize',sizeticks,'interpreter','latex')
        end
        title(yfnames{in},'fontsize',sizelabel,'interpreter','latex') % 24
    end
    
    lgd = legend(legend_string{ifig,:});
    set(lgd,'Position',...
        [0.126041666666667 -0.00975710738586667 0.777788036318183 0.0564541461102402],...
        'Orientation','horizontal','FontSize',18,'interpreter','latex');
    pause(1)
    legend boxoff
end
