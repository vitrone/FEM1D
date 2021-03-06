function logerr_plot(t, err)
PX  =   16;
PY  =   16;
lw  =   2;
ms  =   8;
fs  =   14;

figure1 =   figure('PaperPosition',[0 0 PX PY],...
    'PaperUnits','centimeters',...
    'PaperPositionMode','manual',...
    'PaperType','<custom>',...
    'PaperSize',[PX PY]);
axes1   =   axes('Parent',figure1,...
    'FontName','Times',...
    'FontSize',fs,...
    'YScale','log');
box(axes1,'on');
hold(axes1,'all');
plot(t, err, 'r', 'MarkerSize', ms, 'LineWidth', lw)

ylabel('error','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
axis([0, max(t), min(err), max(err)])