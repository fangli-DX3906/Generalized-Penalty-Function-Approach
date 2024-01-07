clear
close all
clc

cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1'
addpath ('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/code');
addpath ('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/code');
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/code/tools');   
addpath('/Users/fangli/Dropbox/Marco_Fang_project/SimpleNKGPFA/NK model without capital/data')
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/data')
addpath('/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/results')

totalVars = {'Output', 'Consumption', 'Hours', 'Wage', 'SDF', 'Multiplier',... 
                   'Inflation', 'TFP','Goverment spending', ...
                   'Monetary policy shock', 'Discount factor', 'Interest rate'};
aux_names = {'Consumption', 'Hours', 'Wage', 'Interest rate', 'SDF', 'Multiplier'};   
var_names  = {'Output', 'Inflation'};
shock_names = {'demand shock', 'supply shock'};
yfnames = [var_names, aux_names];
var_idx = findingVars(var_names, totalVars);
aux_idx = findingVars(aux_names, totalVars);
neconomy =2000;
maxOrder = 8;
maxAuxOrder =  12;
n_shock = length(shock_names);
H = 40;

% Data read
load economies2000.mat
load IRFs_model.mat
IRFs_model = IRFs_model([1 7 2 3 4 12 5 6 13 19 14 15 16 24 17 18], :);
load IRF_simulated_bayesian.mat

% merge data
% load IRF_3906_1_200.mat
% a = IRF_bayesian;
% load IRF_3908_201_850.mat
% b = IRF_bayesian;
% load IRF_3909_851_1800.mat
% c = IRF_bayesian;
% load IRF_3910_1801_1876.mat
% d =IRF_bayesian;
% load IRF_3907_2000_1875.mat
% e = IRF_bayesian;
% 
% IRF_simulated = cell(neconomy, 1);
% for i=1:neconomy
%     if i<=200
%         IRF_simulated{i} = a{i};
%     elseif 201<=i && i<=850
%         IRF_simulated{i} = b{i};
%     elseif 851<=i && i<=1800
%         IRF_simulated{i} = c{i};
%     elseif 1801<=i && i<=1874
%         IRF_simulated{i} = d{i};
%     else 
%         IRF_simulated{i} = e{i};
%     end
% end
% save '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/IRF_simulated_bayesian.mat' IRF_simulated

% plot
x = randperm(neconomy);
disp(['----------', num2str(x(1)), '------------'])
IRFxxx = IRF_simulated{x(1)};

periods = 0:H;
numVars = size(yfnames,2);
numShocks = length(shock_names);
rows = floor(sqrt(numVars));
cols = ceil(numVars/rows);
horz = H + 1;
colorScheme = {'b','r'};
legendString = {'95\% Confidence Interval $ \ \ \ $', 'BSignRestrction IRF $\ \ \ $', 'Model Implied IRF $ \ \ \ $', 'Demand shock $ \ \ \ $';...
                         '95\% Confidence Interval $ \ \ \ $', 'BSignRestrction IRF $\ \ \ $', 'Model Implied IRF $ \ \ \ $', 'Supply shock $ \ \ \ $'};

for i_shock = 1:numShocks
    if i_shock == 1
        hfig =  findobj('type','figure');
        nfig = length(hfig);
        figure(nfig + i_shock)
    else
        figure(nfig+ i_shock)
    end
    
    for i_var=1: numVars
        subplot(rows, cols, i_var)
        set(gcf,'color','w'); % sets white background color
        if i_shock == 1
            set(gcf, 'Position', [1           1        1440         849]); % sets the figure fullscreen left
        else
            set(gcf, 'Position', [1           1        1440         849]); % sets the figure fullscreen right
        end
        varname = yfnames{i_var};
        shockname = shock_names{i_shock};
        hold on
        x2 = [periods, fliplr(periods)];
        % The bigger CI
        if i_shock == 1
            i_this = i_var;
        else
            i_this = numVars + i_var;
        end
       
        inBetween = [IRFxxx(i_this,1:H+1,4), fliplr(IRFxxx(i_this,1:H+1,16))];
        hh1 = fill(x2, inBetween, [0.5 0.5 0.5], 'LineStyle','none');
        set(hh1,'facealpha',.5) 
        plot(periods, IRFxxx(i_this,1:H+1,10),'linewidth',2,'Color', colorScheme{i_shock})
        plot(periods, IRFs_model(i_this,1:H+1),'-o','linewidth',1,'Color','k')
        plot(periods,zeros(1,H+1), 'Color','k')
        xt = get(gca,'XTick');
        set(gca, 'FontSize', 14)
        title(varname,'fontsize',21,'interpreter','latex') % 24
        pause(0.1)

        num_ticks = get(gca,'ytick')';
        for i = 1:length(num_ticks)
            if abs(num_ticks(i)) < 10^(-12)
                num_ticks(i) = 0;
            end
        end
        
        %grid on
        clear new_yticks num_ticks a pct
        ylabel('Percent','fontsize',15,'interpreter','latex')
    end

    if i_shock == 1
        lgd = legend(legendString{i_shock,:});
    else
        lgd = legend(legendString{i_shock,:});
    end
   
    set(lgd,'Position',...
        [0.126041666666667 -0.00975710738586667 0.777788036318183 0.0564541461102402],...
        'Orientation','horizontal','FontSize',20,'interpreter','latex');
    pause(1)
    legend boxoff
end