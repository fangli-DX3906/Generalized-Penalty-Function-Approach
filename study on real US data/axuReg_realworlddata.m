% Simple Cholesky identification with Efron bootstrap with Kilian correction
% Code by Brianti, M. (parts of the code are from Kilian's website)
% 12/15/2021
% testing auxiliary regression

%% Data read
%%
clear
close all
clc
cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Empirics/thisOne'

%% params value 
%%
system_names = {'GDP', 'CPI', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate', 'TFP', 'Employee'}; 
START = 1960;  % 1960-2021.9.30
END = 2020;
iscon = 1; 
istr = 0; 
% ynames = {'GDP', 'CPI', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate', 'TFP'}; 
% ynames = {'GDP', 'Inflation', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate', 'TFP'}; 
yfnames = {'GDP', 'Inflation', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate', 'TFP'}; % names for the figure
shocknames = {'demand shock', 'supply shock'}; 

%% import the data
%%
[dataset, txt, ~] = xlsread('2021data.xlsx');
tf = isreal(dataset);
if tf == 0
      fprintf('\n')
      warning('Dataset has complex variables in it.')
      dataset = real(dataset);
      fprintf('\n')
end
time_start  = dataset(1,1);
time_end = dataset(end,1);
Time = dataset(:,1);

% Assess names to each variable as an array
for i = 1:size(dataset,2)
      eval([txt{i} ' = dataset(:,i);']);
end

% Standard Transformation
Consumption = log(RealConsumptionNonDurables + RealConsumptionService);
Investment = log(RealConsumptionDurables + RealInvestment);
GDP = log(RealGDP);
CPI = log(CPI);
Def = log(GDPDeflator);
TotalHours = log(TotalHours);
Employee = log(Employee); 
TFP = cumsum(TFP,'omitnan')/400;

% Build System (data set)
for i = 1:length(system_names)
      Y(:,i) = eval(system_names{i});
end

% Truncate data in case of NaN
[XXX, ~, ~] = truncate_data([Y Time]);
Y = XXX(:,1:end-1);
TimeXXX = XXX(:,end);
locSTART = find(TimeXXX == START);
locEND = find(TimeXXX == END);
TimeXXX = TimeXXX(1:locEND);
Y = Y(1:locEND, :);

if isempty(locSTART) == 1
      Y = Y(1:end,:);
      TimeXXX = TimeXXX(1:end);
      fprintf('\n')
      warning(['First available observation is in ',num2str(TimeXXX(1))])
      fprintf('\n')
else
      Y = Y(locSTART:end,:);
      TimeXXX = TimeXXX(locSTART:end);
end

%% different cases
%%
% Y1 = Y(:, 1:end-2);
% system_names1 = {'GDP', 'CPI', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate'};

Inflation = Y(2:end, 2) - Y(1:end-1, 2);  % Y(:,2) is CPI
% Inflation = Y(2:end,2) ./Y(1:end-1, 2) -1;
Y1 = Y(2:end, 1:end-2);
Y1(:, 2) = Inflation*1000;
system_names1 = {'GDP', 'Inflation', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate'}; 

% timeIdx = 1:size(Y, 1);
% XX = [ones(size(Y, 1), 1) timeIdx'];
% Y2 = zeros(size(Y, 1), size(Y, 2)-1);
% for i = 1: size(Y2, 2)
%     [~, ~, Y2(:, i)] = regress(Y(:, i), XX);
% end
% system_names2 = {'GDP', 'CPI', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate', 'TFP'}; 
% 
% Y3 = Y;
% Y3(:, 1) = Y3(:, 1) - Y3(:, end);
% Y3(:, 5) = Y3(:, 5) - Y3(:, end);
% Y3(:, 3) = Y3(:, 3) - Y3(:, end);
% Y3(:, 4) = Y3(:, 4) - Y3(:, end);
% Y3 = Y3(:, 1:end-1);
% system_names3 = {'GDP', 'CPI', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate', 'TFP'}; 

%% estimation
%%
%%%%%%%%%%%%%%%%%%%
y_data = Y1;
systemNames = system_names1;
%%%%%%%%%%%%%%%%%%%

% estimation
maxOrder = 4;
oo = orderSelect(y_data, maxOrder);
p = oo;
[t,q] = size(y_data);
h = 40; 
[A,SIGMA,Uhat,V,X] = olsvarc(y_data, p);	
% A is the companion matrix, same as Bcomp
% SIGMA: VC_eps
% V: cvec
% Uhat: Residmat

% adjust for stationarity if needed
if ~ any(abs(eig(A))>=1)
    disp('*******')
    [A] = asybc(A,SIGMA,t,p);
end
SIGMA = SIGMA(1:q,1:q);

%% IRF
%%
% VAR Impulse response analysis (OIR)
% [IRF]= irfvar(A,SIGMA,p,h);

% switch for PF indetification
PFsign = 1;

% PF
if PFsign
    [~, ~, ~, ~, sigma, ~, mVAR] = estimVARCoeff(y_data, p, iscon, istr, systemNames); 
    initDelta = 0;    
    [GAMMA_, delta, ~] = solver(mVAR, {{'GDP', 'Inflation'}, {'GDP', 'Inflation'}}, [[1,1]; [1,1]], [[1,1]; [1,-1]], initDelta, [1,2], 1, [], 1);  
    [IRF] = irfvar_PF(A, SIGMA, p, h, GAMMA_);
end
 
% VAR bootstrap
% randn('seed',1234);
rng(12);
nrep   = 10;  % 200
[t,q] = size(y_data);
y_data = y_data';
YY = y_data(:,p:t);
for i = 1:p-1
    YY = [YY; y_data(:,p-i:t-i)];
end
U  = Uhat; 

% making data
Ysim = cell(nrep,1);
for r=1:nrep
    r
    Ur = zeros(q*p,t-p);
    Yr = zeros(q*p,t-p+1);
    
    pos     = fix(rand(1,1)*(t-p+1))+1;
    Yr(:,1) = YY(:,pos);

    % iid resampling
    index         = fix(rand(1,t-p)*(t-p))+1;
    Ur(:,2:t-p+1) = U(:,index);

    for i = 2:t-p+1
        Yr(:,i) = V + A*Yr(:,i-1) + Ur(:,i);
    end

    yr = [Yr(1:q,:)];
    for i = 2:p
        yr = [Yr((i-1)*q+1:i*q,1) yr];
    end
    
    Ysim{r} = yr';
    clear Yr Ur
end

% choosing lag order
maxOrder = 4;
optimList = [];
for r=1:nrep
    ydata = Ysim{r};
    ydata = ydata(:, [1 2]);
    disp(['In ', num2str(r), ' data has the rank of ', num2str(rank(ydata))]);
    try
        oo = orderSelect(ydata, maxOrder);
%         oo
    catch
        oo = -1;
    end
    optimList = [optimList, oo];
end

a=0;
b=0;
c=0;
d=0;
for i=1:size(optimList,2)
    if optimList(i)==1
        a=a+1;
    elseif optimList(i)==2
        b=b+1;
    elseif optimList(i)==3
        c=c+1;
    else
        d=d+1;
    end
end
[~, OO] = max([a,b,c,d]);
OO

% IRF
IRFmat = zeros(nrep,25*(h+1));
IRF_aux = [];
aux_systemNames = {'GDP', 'Inflation','Investment', 'Consumption', 'ShadowRate'};

tic
for r=1:nrep
    ydata = Ysim{r}; 
    Yy = ydata(:, [1 2 3 4 6]);
    Ss =  ydata(:, [5]);
    pr = optimList(r);
    
    % GPFA
    if PFsign
        [~, ~, ~, ~, sigma, ~, mVAR] = estimVARCoeff(Yy, pr, iscon, istr, aux_systemNames);    % aux_systemNames
        [Ar,SIGMAr] = olsvarc(Yy,pr);
        if ~ any(abs(eig(Ar))>=1)
            [Ar] = asybc(Ar,SIGMAr,t, pr);
        end
%         initDelta = 0;    
%         [GAMMA, delta, ~] = solver(mVAR, {{'GDP', 'Inflation'}, {'GDP', 'Inflation'}}, [[1,1]; [1,1]], [[1,1]; [1,-1]], initDelta, [1,2], 1, [], 1);  
        [GAMMA, ~] = identifyPFA(mVAR, {{'GDP', 'Inflation'}, {'GDP', 'Inflation'}}, [[1,1]; [1,1]], [[1,1]; [1,-1]], delta, [1,2], 1);
        IRFr = irfvar_PF(Ar, SIGMAr(1:5, 1:5), pr, h, GAMMA);
%         delta        
    else
        [Ar,SIGMAr] = olsvarc(yr,pr);
        if ~ any(abs(eig(Ar))>=1)
            [Ar] = asybc(Ar,SIGMAr,t,pr);
        end
        IRFr = irfvar(Ar,SIGMAr(1:5,1:5),pr,h);
    end
    
    % aux reg
    par = auxorderSelect(Yy, Ss, maxOrder);
    [Aa, ~] = olsvaraux(Ss, Yy, par);
    IRF_new = [];
    for s=1:size(shocknames,2) 
        if s==1
            sv = [1, 2, 3, 4, 5];  % demand shock
        else
            sv = [6 ,7 ,8 ,9, 10];  % supply shock
        end

        for hh = 1: h+1
            if hh <= par+1
                aaa =  IRFr(sv, 1:hh);
                bbb = [];
                for idx=1:hh
                    bbb = [aaa(:, idx); bbb];
                end
                irff = Aa(:, 1:5*hh) * [bbb];
                IRF_new = [IRF_new, irff];
            else
                aaa =  IRFr(sv, hh-par:hh);
                bbb = [];
                for idx=1:par+1
                    bbb = [aaa(:, idx); bbb];
                end
                irff = Aa(:,1:end) * [bbb];
                IRF_new = [IRF_new, irff];
            end
        end   
    end
    
    IRF_temp = [];
    for i=1:size(IRF_new,1)
        IRF_temp = [IRF_temp, IRF_new(i,:)];
    end

    % each line 
    % {y1[shock1(t=1~41) shock2(t=1~41)], y2[shock1(t=1~41) shock2(t=1~41)], ...}    
    IRF_aux = [IRF_aux; IRF_temp];
    IRFmat(r,:) = vec(IRFr)';
    
    disp(['r= ', num2str(r), ' , pr= ', num2str(pr), ' par= ', num2str(par)])
end
toc

% Calculate 68/95 percent interval endpoints
%CI=prctile(IRFmat,[16 84]);
% CI0   = prctile(IRFmat,[0 100]);
CI1   = prctile(IRFmat,[20 80]);
CI11   = prctile(IRFmat, 50);
% CI2   = prctile(IRFmat,[5.0 95.0]);
CILO1 = reshape(CI1(1,:),25,h+1);
CIUP1 = reshape(CI1(2,:),25,h+1);
CImid = reshape(CI11, 25, h+1);
% CILO2 = reshape(CI2(1,:),q^2,h+1);
% CIUP2 = reshape(CI2(2,:),q^2,h+1);


CI1_aux = prctile(IRF_aux,[20 80]);
CI2_aux = prctile(IRF_aux, 50);
CILO1_aux = [];
CIUP1_aux = [];
CImid_aux = [];
for i = 1:1*size(shocknames,2)
    CILO1_aux = [CILO1_aux; CI1_aux(1, (i-1)*41+1:i*41)];
    CIUP1_aux = [CIUP1_aux; CI1_aux(2, (i-1)*41+1:i*41)];
    CImid_aux = [CImid_aux; CI2_aux(1, (i-1)*41+1:i*41)];
end

CILO1 = [CILO1(1:5, :); CILO1_aux(1,:); CILO1(6:10, :);  CIUP1_aux(2,:)];
CIUP1 = [CIUP1(1:5, :); CIUP1_aux(1,:); CIUP1(6:10, :); CIUP1_aux(2,:)];
CIMID = [CImid(1:5, :); CImid_aux(1,:); CImid(6:10,  :);  CIUP1_aux(2,:)];

% plot
periods = 0:h-1;
color = {'b','r'};
for i_shock = 1:2
    if i_shock == 1
        hfig =  findobj('type','figure');
        nfig = length(hfig);
        figure(nfig + i_shock)
    else
        figure(nfig+ i_shock)
    end
    for i_var=1:6
        subplot(2,3,i_var)
        set(gcf,'color','w'); % sets white background color
        if i_shock == 1
            set(gcf, 'Position', [1 41 1920 963]); % sets the figure fullscreen left
        else
            set(gcf, 'Position', [1 41 1920 963]); % sets the figure fullscreen right
        end
        varname = yfnames{i_var};
        shockname = shocknames{i_shock};
        %name = names{i_shock};
        hold on
        x2 = [periods, fliplr(periods)];
        % The bigger CI
        if i_shock == 1
            i_this = i_var;
        else
            i_this = 6 + i_var;
        end
        
       inBetween = [CILO1(i_this,1:h), fliplr(CIUP1(i_this,1:h))];
       hh1 = fill(x2, inBetween, [0.8 0.8 0.8],'LineStyle','none');
       set(hh1,'facealpha',.5)
        
%         inBetween = [CILO2(i_this,1:h), fliplr(CIUP2(i_this,1:h))];
%         hh2 = fill(x2, inBetween, [0.6 0.6 0.6],'LineStyle','none');
%         set(hh2,'facealpha',.5)        
        
        plot(periods,CIMID(i_this,1:h),'linewidth',4,'Color',color{i_shock})
        plot(periods,IRF(i_this,1:h),'linewidth',4,'Color',color{i_shock},'marker','o', 'linestyle','-')
        plot(periods,zeros(1,h), 'Color','k')
        xt = get(gca,'XTick');
        set(gca, 'FontSize', 14)
        title(varname,'fontsize',20,'interpreter','latex') % 24
        pause(0.1)

        num_ticks = get(gca,'ytick')';
        for i = 1:length(num_ticks)
            if abs(num_ticks(i)) < 10^(-12)
                num_ticks(i) = 0;
            end
        end
        
%         a = [cellstr(num2str(num_ticks*100))];
%         % Create a vector of '%' signs
%         pct = char(ones(size(a,1),1)*'\%');
%         % Append the '%' signs after the percentage values
%         new_yticks = [char(a),pct];
%         % 'Reflect the changes on the plot
%         set(gca,'yticklabel',new_yticks,'TickLabelInterpreter','latex')
%         hold off
        
        %grid on
        clear new_yticks num_ticks a pct
        ylabel('Percent','fontsize',15,'interpreter','latex')
    end

    if i_shock == 1
        lgd = legend('95\% confidence interval $ \ \ \ $', 'AuxReg IRF (D)', 'GPFA IRF $ \ \ \ $');
    else
        lgd = legend('95\% confidence interval $ \ \ \ $', 'AuxReg IRF (S)', 'GPFA IRF $ \ \ \ $');
    end
   
    set(lgd,'Position',...
        [0.126041666666667 -0.00975710738586667 0.777788036318183 0.0564541461102402],...
        'Orientation','horizontal','FontSize',20,'interpreter','latex');
    %[0.742708333333333 0.0933940774487471 0.174143115621851 0.359908883826879]
    %[0.126041666666667 -0.00975710738586667 0.777788036318183 0.0564541461102402]
    %[0.7625 0.0882019694529008 0.174143115621851 0.359908883826879]
    pause(1)
    legend boxoff
%     cd(['/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Empirics/thisOne' '/Export_Fig']) 
%     export_fig(shockname)
%     cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Empirics/thisOne'
end

asd


%% OVD
%%
J  = [eye(q,q) zeros(q,q*(p-1))];
TH1 = J*A^0*J'; 
TH = TH1*chol(SIGMA)' * GAMMA_; 
TH = TH'; 
TH2 = (TH.*TH); 
TH3(1,:,:) = TH2;
for i = 2:h
    TH = J*A^(i-1)*J'*chol(SIGMA)' * GAMMA_; 
    TH = TH'; 
    TH2  = (TH.*TH); 
    TH3(i,:,:) = squeeze(TH3(i-1,:,:)) + TH2;
end
TH4 = sum(TH3,2);

% First dimension refers to horizon h
% Second dimension refers to shocks j=1,...,q that explain any given variable
% Third dimension refers to variables whose variation is to be explained
VC  = zeros(h,q,q);
for j = 1:q
    VC(:,j,:) = TH3(:,j,:)./TH4;
end

fig = figure(2); % name figure
fig.Color = 'w'; % white backgroung
fig.WindowState = 'maximized'; % maximize figure size
for i = 1: q^2
    
    if PFsign
        if i > q
            break
        end
    end
    
    subplot(2,3,i)
    %varname = varnames{i_var};
    %shockname = names{i_shock};
    %name = names{i_shock};

    hold on
    x2 = [0:h-1, fliplr(0:h-1)];

    demand_shock = VC(:, 1, i);
    supply_shock = VC(:, 2, i);
    subtotal = demand_shock + supply_shock;
    plot(0:h-1, demand_shock, 'linewidth', 2.5, 'Color', 'b')
    hold on
    plot(0:h-1, supply_shock, 'linewidth', 2.5, 'Color', 'r')
    hold on
    plot(0:h-1, subtotal, 'k--', 'linewidth', 1)

    xt = get(gca,'XTick');
    set(gca, 'FontSize', 12)
    %title([varname],'fontsize',20,'interpreter','latex') % 24
    pause(0.1)
    
    % Convert y-axis values to percentage values by multiplication
%     num_ticks = get(gca,'ytick')';
%     for is = 1:length(num_ticks)
%         if abs(num_ticks(is)) < 10^(-12)
%             num_ticks(is) = 0;
%         end
%     end
%     a = [cellstr(num2str(num_ticks*100))];
%     
%     % Create a vector of '%' signs
%     pct = char(ones(size(a,1),1)*'\%');
%     % Append the '%' signs after the percentage values
%     new_yticks = [char(a),pct];
%     % 'Reflect the changes on the plot
%     set(gca,'yticklabel',new_yticks,'TickLabelInterpreter','latex')
%     %clear new_yticks
%     %axis tight
%     hold off
%     %grid on
%     clear new_yticks num_ticks a pct
    xlabel('Quarter','fontsize',10,'interpreter','latex')
    if i <= q
        title(yfnames{i},'fontsize',18,'interpreter','latex') 
    end
    
end

lgd = legend(' Demand shocks', 'Supply shocks', 'Together');
set(lgd,'Position',[0.241230597012285 0.008080837103293 0.56686861342278 0.046511627906977],...
    'Orientation','horizontal','FontSize',15,'interpreter','latex');
pause(1)
legend boxoff

cd(['/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Empirics/thisOne' '/Export_Fig']) 
export_fig('Variance Decomposition')
cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Empirics/thisOne'