%**************************************************************
% MAIN_PROG - Solves the NK model with capital
%
% Code by Marco Brianti, University of Alberta
% January, 1 2022
% currently working on
%**************************************************************
cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Empirics/thisOne'
clear
close all
clc
disp('NK model with physical capital')

%% pertubation
%%
% y c n w ii m lam pii a t chi betv k r
% col1: output; col8: inflation; col5: investment; col2: consumption;
% col3: total hours; col14: interest rate
po6 = [1 8 5 2 3 14];
po2 = [1 8];
po = po2;
po_aux = [5 2 3 14];

% Load Parameters
[flexp, setp] = paramset(1);

% Compute state-space representation of model
[fyn, fxn, fypn, fxpn, ss]  = model(flexp,setp);
[gx,hx]                      = gx_hx_alt(fyn,fxn,fypn,fxpn);


%% IRFs
%%

% we have only four shocks for now: tech shock, government spending shock,
% monetary policy shock, demand shock
nss                 = 4;

% Construct eta matrix
nx                  = length(hx);
ny                  = size(gx,1);
eta                 = zeros(nx,nss);
eta(1:nss,1:nss)    = eye(nss);   % same variance?
% eta(1,1) = setp.siga;
% eta(2,2) = setp.sigt;
% eta(3,3) = setp.sigc;
% eta(4,4) = setp.sigb;

% Naming variables and figure(s)
subplottitles = {'Output, $y$', 'Consumption, $c$', 'Hours, $n$', 'Wage, $w$', 'Investment, $i$',...
                        'Stochastic discount factor, $m$', 'Multiplier, $\lambda$','Inflation, $\pi$',...
                        'Technology, $a$', 'Goverment spending, $t$','Monetary policy, $\chi$',...
                        'Discount factor, $\beta$', 'Capital, k', 'Interest rate, $r$'};
shocktitles = {'Technology shock','Goverment spending shock','Monetary policy shock', 'Demand Shock'};

% Impulse Responses to shock i in detrended levels
horz   = 41; % length of the impulse response functions

for is = 1:nss
    shockmat     = zeros(1,nss);
    shockmat(is) = 1;
    irfs         = ir(gx,hx,eta*shockmat',horz);
    if is==3
        demandShock = irfs';
        demandShock = demandShock(po6, :);
    elseif is==1
        supplyShock = irfs';
        supplyShock = supplyShock(po6, :);
    end
    figure
    for j=1:length(subplottitles)
        subplot(4,4,j)
        plot([0:horz-1],irfs(:,j))
        title(subplottitles{j},'interprete','latex','fontsize',12)
    end
    sgtitle(shocktitles{is},'interprete','latex','fontsize',24)
end

IRFs_model=[demandShock; supplyShock];

close all

%% simulation
%%
rng(1);
nsim = 2000;
Ysim = cell(nsim,1);
T = 214;
for isim = 1:nsim
    x0 = eta *  [randn(1,1); 0; randn(1,1); 0];    % [randn(1,1); randn(1,1); randn(1,1); randn(1,1)];  %
    pd = length(x0);
    MX = [gx; eye(pd)];
    x= x0;
    for t=1:T
        Y(t, :) = (MX*x)';
        x = hx*x+eta*[randn(1,1); 0; randn(1,1); 0];    % [randn(1,1); randn(1,1); randn(1,1); randn(1,1)];  %
    end
    Ysim{isim} = Y;
end

% selecting the optim order
maxOrder = 4;
optimList = [];
for isim=1:nsim
    ydata = Ysim{isim};
    ydata = ydata(:, po);
    disp(['In ', num2str(isim), ' data has the rank of ', num2str(rank(ydata))]);
    try
        oo = orderSelect(ydata, maxOrder);
    catch
        oo = -1;
    end
    optimList = [optimList, oo];
end     

% ******************************************************************************
% 1. Only two shocks, no matter 4 or 6 var system in deviation, the simulated data are not full rank in col
%     It warns: X is rank deficient to within machine precision. 
% 2. Only two shocks, no matter 4 or 6 var system in level, the simulated data are full rank in col
%     But It also warns the same thing
% 3. Four shocks, 4 var sysetm, no matter in deviation or in level, no problem (full rank and no warning)
% 4. Four shocks, 6 var system,  in deviation, NOT full rank;
%                                              in level, full rank but it warns.
% ******************************************************************************

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
OO % optim lag order is param0=1, param1=4, param2=1

% B) estimation
q = size(po,2);
t = 214;
h = 40;     % should be horiz-1, not include the t=0 period
p = OO;
iscon = 1; 
istr = 0; 
shocknames = {'demand shock', 'supply shock'}; 
% yfnames = {'GDP', 'Inflation', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate', 'TFP'};
% system_names0 = {'GDP', 'Inflation', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate'};
system_names0 = {'GDP', 'Inflation'};
IRFmat = zeros(nsim,(q^2)*(h+1));
rowAvail = [];

for isim=1:nsim
    isim
    Y0 = Ysim{isim};
    Y0 = Y0(:, po);
    
    % Estimation
    y_data = Y0;
    [A,SIGMA,Uhat,V,X] = olsvarc(y_data, p);	
    if ~ any(abs(eig(A))>=1)
        [A] = asybc(A,SIGMA,t,p);
    end
    SIGMA = SIGMA(1:q,1:q);
    
%     SIGMA = round(SIGMA, 5);
%     % VAR Impulse response analysis (OIR)
%     [IRF]= irfvar(A,SIGMA,p,h);
%     eig(SIGMA)

    % switch for PF indetification
    PFsign = 1;

    % PF
    if PFsign
        try
            [~, ~, ~, ~, sigma, ~, mVAR] = estimVARCoeff(y_data, p, iscon, istr, system_names0);  % lag order selection
            if isim == 1 
                initDelta = 0;    
                [GAMMA_, delta, ~] = solver(mVAR, {{'GDP', 'Inflation'}, {'GDP', 'Inflation'}}, [[1,1]; [1,1]], [[1,1]; [1,-1]], initDelta, [1,2], 1, [], 1);
            else
                [GAMMA, ~] = identifyPFA(mVAR, {{'GDP', 'Inflation'}, {'GDP', 'Inflation'}}, [[1,1]; [1,1]], [[1,1]; [1,-1]], delta, [1,2], 1);
            end
            IRF = irfvar_PF(A, SIGMA, p, h, GAMMA_);
            rowAvail = [rowAvail, isim];
        catch
            disp(['In isim = ', num2str(isim), ' covariance matrix is not PD'])
            continue
        end
    end
    
    IRFmat(isim,:) = vec(IRF)';
    % each line:
    % {t=1[shock1(y1 y2), shock2(y1 y2)], t=2[shock1(y1 y2), shock2(y1 y2)], ...}
    
end

% save simulated0209_param1.mat IRFmat;
% save rows0209_param1.mat rowAvail;
% load simulated0209_param1.mat;
% load rows0209_param1.mat;

% axuilary regression
IRF_aux = [];
cl = size(po_aux,2)*size(shocknames,2)+1;
for isim = 1:nsim
    isim
    Ydata = Ysim{isim};
    Yy = Ydata(:, po);
    Ss = Ydata(:, po_aux);
    Aa = olsvaraux(Ss, Yy, p);
    IRFdata = IRFmat(isim, :);
    IRF_ = reshape(IRFdata, q^2, h+1);
    IRF_new = [];
    
    for s=1:size(shocknames,2)
     
        if s==1
            sv = [1, 2];  % demand shock
        else
            sv = [3, 4];  % supply shock
        end
        
        for hh = 1: h+1
            if hh <= p+1
                aaa =  IRF_(sv, 1:hh);
                bbb = [];
                for idx=1:hh
                    bbb = [aaa(:, idx); bbb];
                end
                irff = Aa(:, 1:1+2*hh) * [1; bbb];
                IRF_new = [IRF_new, irff];
            else
                aaa =  IRF_(sv, hh-p:hh);
                bbb = [];
                for idx=1:p+1
                    bbb = [aaa(:, idx); bbb];
                end
                irff = Aa * [1; bbb];
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
end

CI1_aux = prctile(IRF_aux,[30 70]);
CI2_aux = prctile(IRF_aux, 50);
CILO1_aux = [];
CIUP1_aux = [];
CImid_aux = [];
for i = 1:size(po_aux,2)*size(shocknames,2)
    CILO1_aux = [CILO1_aux; CI1_aux(1, (i-1)*41+1:i*41)];
    CIUP1_aux = [CIUP1_aux; CI1_aux(2, (i-1)*41+1:i*41)];
    CImid_aux = [CImid_aux; CI2_aux(1, (i-1)*41+1:i*41)];
end


%% IRF
IRFmat = IRFmat(rowAvail,:);
CI1   = prctile(IRFmat,[30 70]);
CI2   = prctile(IRFmat, 50);
CILO1 = reshape(CI1(1,:),q^2,h+1);
CIUP1 = reshape(CI1(2,:),q^2,h+1);
CImid = reshape(CI2,q^2,h+1);

CILO1 = [CILO1(1:2, :); CILO1_aux([1 3 5 7],:); CILO1(3:4, :); CILO1_aux([2 4 6 8],:)];
CIUP1 = [CIUP1(1:2, :); CIUP1_aux([1 3 5 7],:); CIUP1(3:4, :); CIUP1_aux([2 4 6 8],:)];
CIMID = [CImid(1:2, :); CImid_aux([1 3 5 7],:); CImid(3:4, :); CImid_aux([2 4 6 8],:)];

periods = 0:h-1;
color = {'b','r'};
numSys = 6;
rows = floor(sqrt(numSys));
cols = ceil(numSys/rows);
for i_shock = 1:2
    if i_shock == 1
        hfig =  findobj('type','figure');
        nfig = length(hfig);
        figure(nfig + i_shock)
    else
        figure(nfig+ i_shock)
    end
    for i_var=1:6
%         subplot(rows, cols, i_var)
        subplot(rows, cols, i_var)
        set(gcf,'color','w'); % sets white background color
        if i_shock == 1
            set(gcf, 'Position', [1 41 1920 963]); % sets the figure fullscreen left
        else
            set(gcf, 'Position', [1 41 1920 963]); % sets the figure fullscreen right
        end
        varname = subplottitles{po6(i_var)};
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
%         
%         inBetween = [CILO2(i_this,1:h), fliplr(CIUP2(i_this,1:h))];
%         hh2 = fill(x2, inBetween, [0.6 0.6 0.6],'LineStyle','none');
%         set(hh2,'facealpha',.5)        
        
        plot(periods,CIMID(i_this,1:h),'linewidth',4,'Color',color{i_shock})
        plot(periods,IRFs_model(i_this,1:h),'linewidth',4,'Color',color{i_shock},'marker','o', 'linestyle','none')
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
        lgd = legend('confidence interval $\ \ \ $', 'median IRF $ \ \ \ $','model Implied IRF $ \ \ \ $','Demand shock $ \ \ \ $');
%         lgd = legend('95\% confidence interval $ \ \ \ $','90\% confidence interval $ \ \ \ $','Demand shock $ \ \ \ $');
    else
        lgd = legend('confidence interval  $\ \ \ $', 'median IRF $ \ \ \ $','model Implied IRF $ \ \ \ $','Supply shock $ \ \ \ $');
%         lgd = legend('95\% confidence interval $ \ \ \ $','90\% confidence interval $ \ \ \ $','Supply shock $ \ \ \ $');
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

%% VD
%%
% J  = [eye(q,q) zeros(q,q*(p-1))];
% TH1 = J*A^0*J'; 
% TH = TH1*chol(SIGMA)' * GAMMA_; 
% TH = TH'; 
% TH2 = (TH.*TH); 
% TH3(1,:,:) = TH2;
% for i = 2:h
%     TH = J*A^(i-1)*J'*chol(SIGMA)' * GAMMA_; 
%     TH = TH'; 
%     TH2  = (TH.*TH); 
%     TH3(i,:,:) = squeeze(TH3(i-1,:,:)) + TH2;
% end
% TH4 = sum(TH3,2);
% 
% % First dimension refers to horizon h
% % Second dimension refers to shocks j=1,...,q that explain any given variable
% % Third dimension refers to variables whose variation is to be explained
% VC  = zeros(h,q,q);
% for j = 1:q
%     VC(:,j,:) = TH3(:,j,:)./TH4;
% end
% 
% fig = figure(2); % name figure
% fig.Color = 'w'; % white backgroung
% fig.WindowState = 'maximized'; % maximize figure size
% for i = 1: q^2
%     
%     if PFsign
%         if i > q
%             break
%         end
%     end
%     
%     subplot(2,3,i)
%     %varname = varnames{i_var};
%     %shockname = names{i_shock};
%     %name = names{i_shock};
% 
%     hold on
%     x2 = [0:h-1, fliplr(0:h-1)];
% 
%     demand_shock = VC(:, 1, i);
%     supply_shock = VC(:, 2, i);
%     subtotal = demand_shock + supply_shock;
%     plot(0:h-1, demand_shock, 'linewidth', 2.5, 'Color', 'b')
%     hold on
%     plot(0:h-1, supply_shock, 'linewidth', 2.5, 'Color', 'r')
%     hold on
%     plot(0:h-1, subtotal, 'k--', 'linewidth', 1)
% 
%     xt = get(gca,'XTick');
%     set(gca, 'FontSize', 12)
%     %title([varname],'fontsize',20,'interpreter','latex') % 24
%     pause(0.1)
%     
%     % Convert y-axis values to percentage values by multiplication
% %     num_ticks = get(gca,'ytick')';
% %     for is = 1:length(num_ticks)
% %         if abs(num_ticks(is)) < 10^(-12)
% %             num_ticks(is) = 0;
% %         end
% %     end
% %     a = [cellstr(num2str(num_ticks*100))];
% %     
% %     % Create a vector of '%' signs
% %     pct = char(ones(size(a,1),1)*'\%');
% %     % Append the '%' signs after the percentage values
% %     new_yticks = [char(a),pct];
% %     % 'Reflect the changes on the plot
% %     set(gca,'yticklabel',new_yticks,'TickLabelInterpreter','latex')
% %     %clear new_yticks
% %     %axis tight
% %     hold off
% %     %grid on
% %     clear new_yticks num_ticks a pct
%     xlabel('Quarter','fontsize',10,'interpreter','latex')
%     if i <= q
%         title(yfnames{i},'fontsize',18,'interpreter','latex') 
%     end   
% end