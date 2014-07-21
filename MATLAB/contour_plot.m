function contour_plot(x,t,PSI,options)
% Filled contour plot template.
% options =   struct('levels',(-5:5),...
%                 'magnification_factor',1e3,...
%                 'tolerance',1e-9);
% data to be plotted --> x,t,PSI

[X,T]   =   ndgrid(x,t);
PX  =   9;
PY  =   8;
fs  =   7;

% specs for slides
% PX  =   20;
% PY  =   20;
% lw  =   2;
% ms  =   10;
% fs  =   20;
figure1=figure('PaperPosition',[0 0 PX PY],...
    'PaperUnits','centimeters',...
    'PaperPositionMode','manual',...
    'PaperType','<custom>',...
    'PaperSize',[PX PY],...
    'Visible', 'off');
axes1 = axes('Parent',figure1,...
    'FontName','Times',...
    'FontSize',fs);
box(axes1,'on');
hold(axes1,'all');
set(figure1,'CurrentAxes',axes1)
view(2),
v       =   options.levels;
mag_f   =   options.magnification_factor;
tol     =   options.tolerance;
% v   =   (-5:5);
map =   jetm(length(v)-1);
set(axes1,'CLim',[v(1) v(end)]);
% mag_f   =   1e3;
contour(X,T,mag_f*log10(abs(PSI)+tol),'LineColor',[0 0 0],...
    'LevelList',v,...
    'Fill','on',...
    'Parent',axes1);
colormap(map);
colorbar('peer',axes1,'FontSize',7,...
    'FontName','Times',...
    'YTick',v);
% title(['log-contour plot, Magnication factor' num2str(mag_f)],...
%       'Interpreter','latex')
ylabel('t','Interpreter','latex')
xlabel('x','Interpreter','latex')
axis([min(x),max(x),0,max(t)])



% print -f<figure number> -r<resolution in dpi> -tiff -depsc <file name>
% -tiff attaches 72 dpi tiff preview
% -depsc : colored eps vector format
% -deps  : black and white eps vector format
