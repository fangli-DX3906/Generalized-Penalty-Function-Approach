%%
%%
% Sign-restriction identified VAR with Bayesian estimation using natural conjugate Gaussian
% Inverse Wishart Prior
%
% Diffuse prior for redu-form VAR. See Furlanetto et al 2017 for references
%
% The model can be partially idenfified
%
% Code by Brianti, M. (see Furlanetto et al. 2017, Appendix A)
% Today 06/24/2022
%

%% warm ups
%%
clear
close all
clc
cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/code'
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/code/data')
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/code/tools');   

% Data read
data_read;

% Endogenous variables in the system
ynames                       = {'logGDP','logTFP','logI','logSP500','logCASH',...
    'EBP','MU3','logGDPDEF'}; % names must be the same as in "Variables Transformation"
yfnames                      = {'GDP','TFP','Investment','Stock prices','Corporate cash',...
    'EBP','Uncertainty','Price index'}; % names for the figure
% Position of some variables. Used later on for the
locI = find(strcmp(ynames,'logI')==1);
locY = find(strcmp(ynames,'logGDP')==1);
locZ = find(strcmp(ynames,'logTFP')==1);

% Prepare system of endogenous variables
for i = 1:length(ynames)
    y(:,i) = eval(ynames{i});
end
sdate        = 1982; % first observation of the EBP
edate        = 2020; % Lenza & Primiceri, remove covid
sloc         = find(QUARTER == sdate);
eloc         = find(QUARTER == edate);
y            = y(sloc:eloc,:);
time         = QUARTER(sloc:eloc,:);
rnan         = ~any(isnan(y),2); % remove nan, if any
if sum(rnan) ~= length(rnan)
    warning(['There are some NaN between ',num2str(sdate),' and ',num2str(edate)])
    y            = y(rnan,:);    % remove nan
    time         = time(rnan,:); % data range
end

%%
%%
% Reduced-form VAR with frequentist approach
[T,N]                    = size(y);
p                        = 1;                % VAR lag order
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
else
    disp('*************************')
    disp('Everything seems good')
    disp('*************************')
end

%% Bayesian starts
%%
% Posterior distribution parameters
S        = T*sigmautilde;  % Inverse Wishart parameter 1
if p > 1
    nu       = (T-N)*(p-1);    % Inverse Wishart parameter 2
else
    nu       = T-N;    % Inverse Wishart parameter 2
end

% Simulate impulse response functions
nburn          = 2000; % nburn for the reduced-form parameters to achieve convergence
% nsim           = 2500; % number of simulations
nsim           = 1000; % number of simulations
ntot           = nsim + nburn;
M              = 2; % for each simulation, required number of accepted rotations
H              = 21; % horizon IRFs
signs          = zeros(N,N);
signs(end-3,:) = [-1,+1,-1,+1,+0,+0,+0,+0];   % investment --> since they are sorted, it is useless to add other shocks if interested in uncertainty and financial. It saves time
signs(end-2,:) = [-1,+1,-1,-1,+1,+1,+1,+0];   % uncertainty, positive on corporate cash
signs(end-1,:) = [-1,+1,-1,-1,-1,+1,+1,+0];   % financial, negative on corpoarte cash
signs(end,:)   = [-1,-1,+0,+0,+0,+0,+0,+1];   % supply
h1             = 1; % first horizon to implement the restriction
h2             = 1; % last horizon to implement the restriction. Cumulate from h1 to h2
nunrestricted  = sum(sum(signs==0));
sigmau         = sigmautilde; % initial value for sigmau. Then it will use the drawns
counter2       = 0;

% Pre-allocate matrices
IRF            = zeros(N^2,T-p+1,M*nsim);
VD             = zeros(N^2,H+1,M*nsim);
vd             = zeros(N,N,H+1);
yhat           = zeros(N,T-p,N,M*nsim);

for isim = 1:ntot

    disp(['Simulation ',num2str(isim)])

    % This is the covariance for the posterior density of alpha
    COV    = kron(sigmau,(X'*X)^(-1));

    % Posterior of beta|sigmau,Y,X ~ Normal(bethat,kron(sigmau,(X'X)^(1)))
    bet    = bethat + chol(COV)'*randn(N*(1+N*p),1);
    B      = (reshape(bet,1+N*p,N));
    usim   = Y - X*B;

    % Draw from the posterior of sigmau|Y,X ~ IW(S,nu)
    sigmau = inv(wish(inv(S),nu));

    if isim > nburn
        % Calculate IRFs
        id            = eye(N*p);
        BB            = [B(2:end,:)'; id(1:end-N,:)];
        counter       = 0;
        while counter < M 
            [Q,R]               = qr(randn(N,N)); % draw orthogonal matrix Q
            Q                   = sign(diag(R)).*Q; % flip the sign of the columns of Q if elements along the main diagonal of R are negative
            candidate           = irfvar(BB,sigmau,p,H,Q); % candidate impulse response functions
            candsim             = sum(candidate(:,h1:h2),2); % sign of the candidate impulse response functions. Also take care of the horizon
            candsimMAT          = (reshape(candsim,N,N))';
            candsimMAT(:,locI)  = candsimMAT(:,locI) - candsimMAT(:,locY);
            candsimMAT(:,locZ)  = candsimMAT(:,locZ) - candsimMAT(:,locY);
            signsimMAT          = sign(candsimMAT);
            [signsimMAT,id]     = sortrows(signsimMAT,'descend'); % sort signsimMAT as "signs" above ???
            if sum(sum((signsimMAT - signs).^2)) == nunrestricted
                counter           = counter + 1;
                counter2          = counter2 + 1;
                Q                 = Q(:,id); % sort column-wise Q, accordingly

                % Calculate impulse responses
                IRF(:,:,counter2) = irfvar(BB,sigmau,p,T-p,Q); % re-calculare IRFs with Q sorted

                % Calculate variance decomposition
                IRsq              = reshape(IRF(:,1:H+1,counter2).^2,N,N,H+1);
                IRsqsumH          = cumsum(IRsq,3);
                totFEV            = squeeze(sum(IRsqsumH,2));
                for ih = 1:size(totFEV,2)
                    vd(:,:,ih) = IRsqsumH(:,:,ih)./totFEV(:,ih);
                end
                VD(:,:,counter2)  = reshape(vd,N^2,H+1);

                % Calculate historical decomposition
                shocksim          = (chol(sigmau)'*Q) \ (usim'); %compute structural shocks shocksim from reduced form shocks usim
                for it = 1:T-p
                    for is = 1:N
                        for iv = 1:N
                            yhat(iv,it,is,counter2) = dot(IRF(iv+(is-1)*N,1:it,counter2),shocksim(is,it:-1:1));
                        end
                    end
                end


            end
        end
    end
end

% Select the credible intervals and the median
IR   = quantile(100*IRF,[0.10 0.16 0.5 0.84 0.90],3);
FEVD = quantile(100*VD,[0.10 0.16 0.5 0.84 0.90],3);
YHAT = median(yhat,4);

load IR.mat
load FEVD.mat
load YHAT.mat
% Calculate the variance decomposition for the median as Furlanetto et al.
IRsq              = reshape(IR(:,1:H+1,3).^2,N,N,H+1);
IRsqsumH          = cumsum(IRsq,3);
totFEV            = squeeze(sum(IRsqsumH,2));
for ih = 1:size(totFEV,2)
    vd(:,:,ih) = IRsqsumH(:,:,ih)./totFEV(:,ih);
end
FEVDmed           = reshape(100*vd,N*N,H+1);
%FEVD              = FEVD - FEVD(:,:,3) + FEVDmed;

% Figure impulse response functions
ssnames   = {'demand','investment','uncertainty','financial','supply'};
sizelabel = 24;
sizeticks = 18;
for ifig = 3:4
    fig             = figure(ifig); % name figure
    fig.Color       = 'w'; % white backgroung
    %fig.WindowState = 'maximized'; % maximize figure size
    fig.Position    = [-1681 40 1431 954]; % 
    if ifig == 3
        Colrs = [1 0 0];
    else Colrs = [0 0 1];
    end
    for in = 1:N
        subplot(2,ceil(N/2),in)
        plot_xtick = [0 5:5:H-1];   % x-axis ticks
        set(gca,'XTick',plot_xtick)
        grid on
        %varname = varnames{i_var};
        %shockname = names{i_shock};
        %name = names{i_shock};
        ic = in + N*(ifig-1) + (N-5)*N;
        hold on
        x2 = [0:H-1, fliplr(0:H-1)];
        % The bigger CI
        inBetween = [IR(ic,1:H,1), fliplr(IR(ic,1:H,end))];
        hh1 = fill(x2, inBetween,Colrs,'LineStyle','none');
        set(hh1,'facealpha',.05)
        % The smaller CI
        inBetween = [IR(ic,1:H,2), fliplr(IR(ic,1:H,end-1))];
        hh2 = fill(x2, inBetween,Colrs,'LineStyle','none');
        set(hh2,'facealpha',.2)
        plot(0:H-1,IR(ic,1:H,3),'-','linewidth',3,'Color',Colrs)
        plot(0:H-1,zeros(1,H),'-','Color','k')
        xt = get(gca,'XTick');
        set(gca, 'FontSize', sizeticks)
        %title([varname],'fontsize',20,'interpreter','latex') % 24
        set(gca,'TickLabelInterpreter','latex')
        axis tight
        hold off
        %grid on
        clear new_yticks num_ticks a pct
        if ic == 1 + N*(ifig-1) + (N-5)*N
            xlabel('Quarter','fontsize',sizeticks,'interpreter','latex')
            ylabel('Percent','fontsize',sizeticks,'interpreter','latex')
        end
        title(yfnames{in},'fontsize',sizelabel,'interpreter','latex') % 24
    end

    % Print figure
    tt         = now;
    DateString = datestr(tt);
    newStr = strrep(DateString,'-','_');
    newStr = strrep(newStr,' ','_');
    newStr = strrep(newStr,':','_');
    %print([pwd, '/figures/','figure',num2str(ifig),DateString],'-depsc','-r0')
    print([pwd, '\figures\','figure_IRFs_',ssnames{ifig},'_shock_',newStr],'-dpng','-r0')

end

%save workspace_13JUL2022_baseline_IRF_only_F_U
save 'IR.mat' IR
save 'FEVD.mat' FEVD
save 'YHAT.mat' YHAT