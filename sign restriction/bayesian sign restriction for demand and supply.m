%% Sign-restriction identified VAR with Bayesian estimation using natural conjugate Gaussian
%%
% VAR coefficient piror: prior gaussian dist
% VAR covariance piror: inverse whishart dist 
% Diffuse prior for redu-form VAR. See Furlanetto et al 2017 for references
% The model can be partially idenfified
% Code by Brianti, M. (see Furlanetto et al. 2017, Appendix A)
% Today 09/20/2022

%% warm ups
%%
clear
close all
clc

cd '/Users/fangli/dissertation/chapter1/code/sign restriction matlab'
addpath('/Users/fangli/dissertation/chapter1/code/sign restriction matlab/tools');   
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation Data/CHAPTER_1/SimpleNK GPFA data')
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation Data/CHAPTER_1/Bayesian GPFA data')
addpath('/Users/fangli/dissertation/chapter1/code/bayesian GPFA matlab/distributions')
addpath('/Users/fangli/dissertation/chapter1/code/bayesian GPFA matlab/estimation')

% Data read
load economies2000.mat
load IRFs_model.mat
IRFs_model = IRFs_model([1 7 2 3 4 12 5 6 13 19 14 15 16 24 17 18], :);

% Endogenous variables in the system
totalVars = {'Output', 'Consumption', 'Hours', 'Wage', 'SDF', 'Multiplier',... 
                   'Inflation', 'TFP','Goverment spending', ...
                   'Monetary policy shock', 'Discount factor', 'Interest rate'};
aux_names = {'Consumption', 'Hours', 'Wage', 'Interest rate', 'SDF', 'Multiplier'};   
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
%     p = orderSelect(y, maxOrder);
    p=4;
    [A,~,Uhat,const,X,SIGMA] = olsvarc(y,p);	 % Frequentist approach: VAR with intercept
     [a, b, c, d] = calcReducedVARParams(y, p);
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
IRF_bayesian = cell(neconomy,1);

% p = parpool(4);
% parfor n = 1:2000
for n = 1:2000
    disp(['---> Estimating ', num2str(n), ' economy <---' ]);
    ydata = Ysim{n};
    y = ydata(:, var_idx);
    s = ydata(:, aux_idx);
    temp = bayesianSignRestriction(y, s, maxOrder, maxAuxOrder, nburn, nsim, M, H, signs, h1, h2);
    IRF = quantile(temp, [0 0.005 0.01 0.025 0.05 0.10 0.16 0.2 0.32 0.5 0.68 0.8 0.84 0.90 0.95 0.975 0.99 0.995 1], 3);
    IRF(:,:,end+1) = mean(temp,3);
    IRF_bayesian{n} = IRF;
end
% delete(p);

% save '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/IRF_bayesian.mat' IRF_bayesian

%% Plot
%%
tic
Y0 = Ysim{1024};
Yy = Y0(:, var_idx);
Ss =  Y0(:, aux_idx);
IRFxxx = bayesianSignRestriction(Yy, Ss, maxOrder, maxAuxOrder, nburn, nsim, M, H, signs, h1, h2);
toc

% Select the credible intervals and the median
% IR   = quantile(IRFxxx,[0.025 0.5 0.975],3);
% FEVD = quantile(100*VD,[0.10 0.16 0.5 0.84 0.90],3);
% YHAT = median(yhat,4);

% Calculate the variance decomposition for the median as Furlanetto et al.
% IRsq              = reshape(IR(:,1:H+1,3).^2,N,N,H+1);
% IRsqsumH          = cumsum(IRsq,3);
% totFEV            = squeeze(sum(IRsqsumH,2));
% for ih = 1:size(totFEV,2)
%     vd(:,:,ih) = IRsqsumH(:,:,ih)./totFEV(:,ih);
% end
% FEVDmed           = reshape(100*vd,N*N,H+1);
%FEVD              = FEVD - FEVD(:,:,3) + FEVDmed;

% Figure impulse response functions
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
    for in = 1: size(yfnames,2)
        subplot(2,ceil(size(yfnames,2)/2),in)
        plot_xtick = [0 5:5:H-1];   % x-axis ticks
        set(gca,'XTick',plot_xtick)
        grid on
        %varname = varnames{i_var};
        %shockname = names{i_shock};
        %name = names{i_shock};
        ic = in + size(yfnames,2)*(ifig-1);
        hold on
        x2 = [0:H-1, fliplr(0:H-1)];
        % The bigger CI
        inBetween = [IRFxxx(ic,1:H,4), fliplr(IRFxxx(ic,1:H,16))];
        hh1 = fill(x2, inBetween,Colrs,'LineStyle','none');
        set(hh1,'facealpha',.05)
%         % The smaller CI
%         inBetween = [IR(ic,1:H,6), fliplr(IR(ic,1:H,8))];
%         hh2 = fill(x2, inBetween,Colrs,'LineStyle','none');
%         set(hh2,'facealpha',.2)
        plot(0:H-1,IRFxxx(ic,1:H,10),'-','linewidth',3,'Color',Colrs)
        plot(0:H-1,IRFs_model(ic,1:H), '--o','linewidth',1,'Color','k')
        plot(0:H-1,zeros(1,H),'-','Color','k')
        xt = get(gca,'XTick');
        set(gca, 'FontSize', sizeticks)
        %title([varname],'fontsize',20,'interpreter','latex') % 24
        set(gca,'TickLabelInterpreter','latex')
        axis tight
        hold off
        %grid on
        clear new_yticks num_ticks a pct
        if ic == 1 + size(yfnames,2)*(ifig-1);
            xlabel('Quarter','fontsize',sizeticks,'interpreter','latex')
            ylabel('Percent','fontsize',sizeticks,'interpreter','latex')
        end
        title(yfnames{in},'fontsize',sizelabel,'interpreter','latex') % 24
    end

%     % Print figure
%     tt         = now;
%     DateString = datestr(tt);
%     newStr = strrep(DateString,'-','_');
%     newStr = strrep(newStr,' ','_');
%     newStr = strrep(newStr,':','_');
%     %print([pwd, '/figures/','figure',num2str(ifig),DateString],'-depsc','-r0')
%     print([pwd, '\figures\','figure_IRFs_',ssnames{ifig},'_shock_',newStr],'-dpng','-r0')
end