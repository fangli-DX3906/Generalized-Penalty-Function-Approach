% Figure impulse response functions
ssnames = {'monetary','demand','investment','uncertainty','financial','supply'};
for ifig = 1:6
    fig             = figure(ifig); % name figure
    fig.Color       = 'w'; % white backgroung
    fig.WindowState = 'maximized'; % maximize figure size
    for in = 1:N
        subplot(2,ceil(N/2),in)
        %varname = varnames{i_var};
        %shockname = names{i_shock};
        %name = names{i_shock};
        ic = in + N*(ifig-1) + (N-6)*N;
        hold on
        x2 = [0:H-1, fliplr(0:H-1)];
        % The bigger CI
        inBetween = [IR(ic,1:H,1), fliplr(IR(ic,1:H,end))];
        hh1 = fill(x2, inBetween, [0 0 0],'LineStyle','none');
        set(hh1,'facealpha',.05)
        % The smaller CI
        inBetween = [IR(ic,1:H,2), fliplr(IR(ic,1:H,end-1))];
        hh2 = fill(x2, inBetween, [0 0 0],'LineStyle','none');
        set(hh2,'facealpha',.2)
        plot(0:H-1,IR(ic,1:H,3),'-','linewidth',3,'Color','k')
        plot(0:H-1,zeros(1,H),'-','Color','b')
        xt = get(gca,'XTick');
        set(gca, 'FontSize', 12)
        %title([varname],'fontsize',20,'interpreter','latex') % 24
        set(gca,'TickLabelInterpreter','latex')
        axis tight
        hold off
        %grid on
        clear new_yticks num_ticks a pct
        if ic == 1 + N*(N-2+ifig-1)
            xlabel('Quarter','fontsize',12,'interpreter','latex')
            ylabel('Percent','fontsize',12,'interpreter','latex')
        end
        title(yfnames{in},'fontsize',16,'interpreter','latex') % 24
    end
    lgd = legend('80\% confidence interval $ \ \ \ $','68\% confidence interval $ \ \ \ $','Point estimate $ \ \ \ $');
    set(lgd,'Position',[0.126041666666667 -0.00975710738586667 0.777788036318183 0.0564541461102402],...
        'Orientation','horizontal','FontSize',20,'interpreter','latex');
    pause(1)
    legend boxoff

    tt         = now;
    DateString = datestr(tt);
    newStr = strrep(DateString,'-','_');
    newStr = strrep(newStr,' ','_');
    newStr = strrep(newStr,':','_');
    %print([pwd, '/figures/','figure',num2str(ifig),DateString],'-depsc','-r0')
    print([pwd, '\figures\','figure_IRFs_',ssnames{ifig},'_shock_',newStr],'-dpng','-r0')

end