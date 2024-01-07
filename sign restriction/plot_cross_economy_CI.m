function  plot_cross_economy_CI(IRF, var_name, shock_name, color_scheme, H, figtitle, legend_string, ifsave, save_name)

sizelabel = 20;
sizeticks = 18;
n_shock = length(shock_name);
n_var = length(var_name);
ci_lower_1 = 1;
ci_upper_1 = 2;
ci_lower_2 = 3;
ci_upper_2 = 4;
theoretical = 5;
median = 6;
point_est  = 7;

for ifig = 1:n_shock
    fig             = figure(ifig); % name figure
    fig.Color       = 'w'; % white backgroung
    %fig.WindowState = 'maximized'; % maximize figure size
    fig.Position    = [-1681 40 1431 954]; 
    
    for in = 1: n_var
        subplot(1, n_var, in)
        if H>=20
            xx = 5;
        elseif H>=10 && H<20
            xx = 2;
        else
            xx = 1;
        end
        plot_xtick = [0 xx:xx:H];   % x-axis ticks
        set(gca,'XTick',plot_xtick)
        grid on
        ic = in + n_var*(ifig-1);
        hold on
        x2 = [0:H, fliplr(0:H)];
        
        % The bigger CI
        inBetween = [IRF(ic,1:H+1, ci_lower_1), fliplr(IRF(ic,1:H+1, ci_upper_1))];
        hh1 = fill(x2, inBetween,'k','LineStyle','none');
        set(hh1,'facealpha',.35)
%         plot(0:H, IRF(ic,:, min_pos), '-','linewidth',1,'Color',color_scheme{ifig})
%         plot(0:H, IRF(ic,:, max_pos), '-','linewidth',1,'Color',color_scheme{ifig})

        % The smaller CI
        inBetween = [IRF(ic,1:H+1, ci_lower_2), fliplr(IRF(ic,1:H+1, ci_upper_2))];
        hh1 = fill(x2, inBetween,color_scheme{ifig},'LineStyle','none');
        set(hh1,'facealpha',.25)
%         plot(0:H, IRF(ic,:, ci95_max_pos), '-','linewidth',1,'Color',color_scheme{ifig})
%         plot(0:H, IRF(ic,:, ci95_min_pos), '-','linewidth',1,'Color',color_scheme{ifig})

%         plot(0:H, IRF(ic,:, median_pos),'-','linewidth',2,'Color',color_scheme{ifig})
%         plot(0:H, IRF(ic,:, mean_pos), '-o','linewidth',2,'Color',color_scheme{ifig})
        plot(0:H, IRF(ic,1:H+1, theoretical), '-d','linewidth',2,'Color',color_scheme{ifig})        
%         plot(0:H, IRF(ic,:, 4), '-o','linewidth',2,'Color','k')
        plot(0:H, IRF(ic,1:H+1, median), '--','linewidth',2,'Color','k')        
        plot(0:H, IRF(ic,1:H+1, point_est), '-','linewidth',2,'Color','k')        
        plot(0:H, zeros(1,H+1),'-','Color','k')
        
        xt = get(gca,'XTick');
        set(gca, 'FontSize', sizeticks)
        %title([varname],'fontsize',20,'interpreter','latex') % 24
        set(gca,'TickLabelInterpreter','latex')
        axis tight
        hold off
        %grid on
        clear new_yticks num_ticks a pct
        
        if ic == 1 + n_var*(ifig-1)
            xlabel('Quarter','fontsize',sizeticks,'interpreter','latex')
            ylabel('Percent','fontsize',sizeticks,'interpreter','latex')
        end
        title(var_name{in},'fontsize',sizelabel,'interpreter','latex')
    end
    sgtitle(figtitle{ifig}, 'fontsize', sizelabel+7,'interpreter','latex')
    
    lgd = legend(legend_string{ifig,:});
    set(lgd,'Position',...
        [0.126041666666667 -0.00975710738586667 0.777788036318183 0.0564541461102402],...
        'Orientation','horizontal','FontSize',18,'interpreter','latex');
    pause(1)
    legend boxoff
    
    if ifsave==true
        cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1'
        export_fig([shock_name{ifig}, '_', save_name])
        cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Dissertation/CHAPTER_1/code'
    end
end

end