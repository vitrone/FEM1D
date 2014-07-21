clc
clear
Nt = 2^14;
c = 10;
fileStr = 'CQ2_TBC_c10';
[err, t] = LSE.test_quadW_CQ2(c, Nt);
%%
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
plot(t, err(1,:), 'r', 'LineWidth', lw)
plot(t, err(2,:), 'g', 'LineWidth', lw)
plot(t, err(3,:), 'b', 'LineWidth', lw)

ylabel('error','Interpreter','latex')
xlabel('$\tau$','Interpreter','latex')
legend('Trapezoidal',...
       'BDF-1',...
       'BDF-2')
   
axis([0,max(t),1e-8,1e-2])
legend1 =   legend(axes1,'show');
set(legend1, 'Interpreter','latex','Location','NorthEast');

myfig = gcf;
tmpStr = [fileStr, '.eps'];
print(myfig, '-r1000', '-tiff', '-depsc', tmpStr);
tmpStr = [fileStr, '.pdf'];
print(myfig, '-r1000', '-dpdf', tmpStr);
%%
K = (1:100:Nt);
M = [(t(K))', (err(:, K))'];
file_name = sprintf('%s.dat',fileStr);
dlmwrite(file_name, M, 'delimiter', '\t', 'precision', '%.9g')


