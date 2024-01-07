function plots(subsequentPeriods, colorScheme, varNames, shockNames, IRFmat, plotCI, CImat, legendString)

periods = 0:subsequentPeriods;
numVars = size(varNames,2);
numShocks = size(shockNames, 2);
rows = floor(sqrt(numVars));
cols = ceil(numVars/rows);
horz = subsequentPeriods + 1;

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
        varname = varNames{i_var};
        shockname = shockNames{i_shock};
        hold on
        x2 = [periods, fliplr(periods)];
        % The bigger CI
        if i_shock == 1
            i_this = i_var;
        else
            i_this = numVars + i_var;
        end
        
        if plotCI
            for n = 1: size(CImat, 1)
               thisCILO = CImat{n,1};
               thisCIUP = CImat{n,2};
               inBetween = [thisCILO(i_this, 1:horz), fliplr(thisCIUP(i_this, 1:horz))];
               hh1 = fill(x2, inBetween, [0.8 0.8 0.8], 'LineStyle','none');
               set(hh1,'facealpha',.5) 
            end
        end
        
        N = size(IRFmat, 2);
        for n = 1:N
            IRF_temp = IRFmat{n};
            if n==1
                plot(periods, IRF_temp(i_this, 1:horz),'linewidth',5,'Color', colorScheme{i_shock})
            else
                plot(periods, IRF_temp(i_this, 1:horz),'linewidth',5,'Color','k', 'linestyle','--')
            end
        end
        
        plot(periods,zeros(1,subsequentPeriods+1), 'Color','k')
        xt = get(gca,'XTick');
        set(gca, 'FontSize', 14)
        title(varname,'fontsize',30,'interpreter','latex') % 24
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
%     cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Empirics/SimpleNKGPFA/NK model without capital/Export_Fig'
%     export_fig(shockname)
%     cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Empirics/SimpleNKGPFA/NK model without capital'
end

end

