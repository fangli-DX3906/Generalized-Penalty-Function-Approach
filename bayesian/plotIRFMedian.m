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
        
        % GPFA 95% CI
        inBetween = [IRFs_GPFA1(ic,1:H+1,1), fliplr(IRFs_GPFA1(ic,1:H+1,3))];
        hh1 = fill(x2, inBetween,Colrs,'LineStyle','none');
        set(hh1,'facealpha', .25)
        
        inBetween = [IRFs_SIGN1(ic,1:H+1,1), fliplr(IRFs_SIGN1(ic,1:H+1,3))];
        hh2 = fill(x2, inBetween,Colrs,'LineStyle','none', 'EdgeColor', 'k', 'Marker', 'o');
        set(hh2,'facealpha', 0)
        
%         inBetween = [IRFs_sign(ic,1:H+1,1), fliplr(IRFs_sign(ic,1:H+1,3))];
%         hh2 = fill(x2, inBetween,Colrs,'LineStyle','none');
%         set(hh2,'facealpha', .5)

%         plot(0:H, IRFs_sign(ic,1:H+1,1),'o','linewidth',3,'Color','k')   
%         plot(0:H, IRFs_sign(ic,1:H+1,3),'o','linewidth',3,'Color','k')   

        plot(0:H, IRFs_GPFA1(ic,1:H+1,2),'-','linewidth',3,'Color',Colrs)   % median GPFA IRF
        plot(0:H, IRFs_model(ic,1:H+1),'-','linewidth',3,'Color', 'k')   % model-implied IRF
        plot(0:H, IRFs_SIGN1(ic,1:H+1,2),'-o','linewidth',3,'Color',Colrs)    % median sign restriction IRF

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
    lgd = legend(' 95\% GPFA confid intervl $ \ \ \ $', ' 95\% SR confid intervl $ \ \ \ $', ' median GPFA IRF $ \ \ \ $', ' model-implied IRF $ \ \ \ $', ' SR median IRF $ \ \ \ $');
    set(lgd,'Position', [0.126041666666667 -0.00975710738586667 0.777788036318183 0.0564541461102402],...
         'Orientation','horizontal','FontSize',20,'interpreter','latex');
    pause(1)
    legend boxoff
%     cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/bayesian GPFA/figures'
%     export_fig(shock_names{ifig})
%     cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/bayesian GPFA'
end
