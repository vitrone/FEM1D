function curve_3d(x,t,PSI,K)
% 3d curve display
PX  =   9;
PY  =   8;
fs  =   7;

% specs for slides
% PX  =   20;
% PY  =   20;
% lw  =   2;
% ms  =   10;
% fs  =   20;
figure1 =   figure('PaperPosition',[0 0 PX PY],...
    'PaperUnits','centimeters',...
    'PaperPositionMode','manual',...
    'PaperType','<custom>',...
    'PaperSize',[PX PY]);
axes1   =   axes('Parent',figure1,...
    'FontName','Times',...
    'FontSize',fs);
box(axes1,'off');
hold(axes1,'all');
set(figure1,'CurrentAxes',axes1)

ylabel('t','Interpreter','latex')
xlabel('x','Interpreter','latex')
zlabel('Intensity','Interpreter','latex');

view(40,80);
Z   =   ones(size(x));
for k=1:max(size(K)),
    plot3(x,t(K(k))*Z,abs(PSI(:,K(k))).^2,'b');
end
hold off