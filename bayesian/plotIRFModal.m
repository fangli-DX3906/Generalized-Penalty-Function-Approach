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
        
        % GPFA 95% Credible set
        inBetween = [IRFs_GPFA2(ic,1:H+1,2), fliplr(IRFs_GPFA2(ic,1:H+1,3))];
        hh1 = fill(x2, inBetween,Colrs,'LineStyle','none');
        set(hh1,'facealpha', .25)
        
        % SR 95% Credible set
        inBetween = [IRFs_SIGN2(ic,1:H+1,2), fliplr(IRFs_SIGN2(ic,1:H+1,3))];
        hh2 = fill(x2, inBetween,Colrs,'LineStyle','none', 'EdgeColor', 'k', 'Marker', 'o');
        set(hh2,'facealpha', 0)
        
%         inBetween = [IRFs_sign(ic,1:H+1,1), fliplr(IRFs_sign(ic,1:H+1,3))];
%         hh2 = fill(x2, inBetween,Colrs,'LineStyle','none');
%         set(hh2,'facealpha', .5)

%         plot(0:H, IRFs_sign(ic,1:H+1,1),'o','linewidth',3,'Color','k')   
%         plot(0:H, IRFs_sign(ic,1:H+1,3),'o','linewidth',3,'Color','k')   

        plot(0:H, IRFs_GPFA2(ic,1:H+1,1),'-','linewidth',3,'Color',Colrs)   % modal GPFA IRF
        plot(0:H, IRFs_model(ic,1:H+1),'-','linewidth',3,'Color', 'k')   % model-implied IRF
        plot(0:H, IRFs_SIGN2(ic,1:H+1,1),'-o','linewidth',3,'Color',Colrs)    % modal sign restriction IRF

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
    lgd = legend(' 95\% GPFA credible sets $ \ \ \ $', ' 95\% SR credible sets $ \ \ \ $', ' modal GPFA IRF $ \ \ \ $', ' model-implied IRF $ \ \ \ $', ' SR modal IRF $ \ \ \ $');
    set(lgd,'Position', [0.126041666666667 -0.00975710738586667 0.777788036318183 0.0564541461102402],...
         'Orientation','horizontal','FontSize',20,'interpreter','latex');
    pause(1)
    legend boxoff
%     cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/bayesian GPFA/figures'
%     export_fig([shock_names{ifig}, 'IK'])
%     cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/bayesian GPFA'
end
