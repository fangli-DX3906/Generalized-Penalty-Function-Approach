function plotOneIR(IR, CI, irHoriz, varName, shockName, isGrid)

xaxis = 1:irHoriz;

plot(xaxis, IR, 'linewidth', 2, 'Color', 'r');
hold on 
plot(xaxis, zeros(size(1:irHoriz)), 'linewidth', 2, 'Color', 'b');

if ~isnan(CI)
    hold on
    plot(xaxis, CI,  'r--');
    hold on 
    plot(xaxis, zeros(size(xaxis)), 'k--')
end

xlabel('Quarter','fontsize',12,'interpreter','latex');
if ~isempty(varName)
    str = strcat(varName);
    title(str, 'interpreter', 'latex');
end

if isGrid == 1
    grid on
end

end

