function fig=draw_burden_rho90_contours( ...
        catalog,results,cfg,ax)

if nargin<4 || ~isgraphics(ax,'axes')
    error('An axes handle is required.');
end
fig=ancestor(ax,'figure');
hold(ax,'on');

gamma = logspace(log10(0.2),log10(cfg.coupled.gammaPlotMax),241);
gB = linspace(0,0.8,161);
[GAMMA,GB] = meshgrid(gamma,gB);
levels = cfg.coupled.rho90ContourLevels;

theta = cfg.analyticalTheta(:);
t90E = [analytical_t90(theta,cfg.coupled.piE,'scarce'), ...
    analytical_t90(theta,cfg.coupled.piE,'parallel')];
t90B = [analytical_t90(theta,cfg.coupled.piB,'scarce'), ...
    analytical_t90(theta,cfg.coupled.piB,'parallel')];
rhoFrontier = GAMMA*(t90B(1)/t90E(1));
rhoParallel = GAMMA*(t90B(2)/t90E(2));

[~,h1] = contour(ax,GAMMA,GB,rhoFrontier,levels,'LineWidth',1.8, ...
    'LineStyle','-','LineColor',[0.12 0.32 0.55]);
[~,h2] = contour(ax,GAMMA,GB,rhoParallel,levels,'LineWidth',1.1, ...
    'LineStyle','--','LineColor',[0.85 0.33 0.10]);
for j = 1:numel(levels)
    text(ax,levels(j),0.615,sprintf('%g',levels(j)), ...
        'FontName','Arial','FontSize',7,'Color',[0.08 0.08 0.08], ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'Rotation',90,'BackgroundColor','w','Margin',0.5, ...
        'Clipping','on');
end

idx = find(strcmp({catalog.group},'core grid'));
x = nan(numel(idx),1);
y = x;
z = x;
for k = 1:numel(idx)
    met = results{idx(k)}.metrics;
    x(k) = median(met.gamma,'omitnan');
    y(k) = median(met.gB,'omitnan');
    z(k) = median(met.rho90_0,'omitnan');
end
sc = scatter(ax,x,y,42,z,'filled','MarkerEdgeColor','k','LineWidth',0.6);
colormap(ax,parula(256));
caxis(ax,[min(gamma) max(gamma)]);
set(ax,'ColorScale','log');
cb = colorbar(ax);
cb.Ticks = [0.25 0.5 1 2 4 8];
cb.TickLabels = {'0.25','0.5','1','2','4','8'};
cb.FontName = 'Arial';
cb.FontSize = 9;
cb.Label.Interpreter = 'latex';
cb.Label.String = 'Simulation median, $\rho_{90}^{0}$';
cb.Label.FontSize = 9;

set(ax,'XScale','log','XTick',[0.25 0.5 1 2 4 8], ...
    'XTickLabel',{'0.25','0.5','1','2','4','8'}, ...
    'FontName','Arial','FontSize',9,'LineWidth',0.8, ...
    'Box','off','TickDir','out','Layer','top');
xlim(ax,[0.2 cfg.coupled.gammaPlotMax]);
ylim(ax,[0 0.9]);
xlabel(ax,'Relative intrinsic repair burden, $\gamma=\Gamma_B/\Gamma_E$', ...
    'Interpreter','latex','FontSize',9);
ylabel(ax,'Effective electricity-dependent service share, $g_B$', ...
    'Interpreter','latex','FontSize',9);
grid(ax,'off');
ax.TitleFontSizeMultiplier = 1;
ax.LabelFontSizeMultiplier = 1;
lgd = legend(ax,[h1 h2 sc],{'Frontier-limited','Fully parallel', ...
    'Simulation medians'},'Location','northwest','NumColumns',1, ...
    'FontName','Arial','FontSize',8,'Box','off');
lgd.NumColumns = 3;
lgd.ItemTokenSize = [10 8];

end

function t90 = analytical_t90(theta,piVec,regime)
model=struct('meanDuration',1,'releaseSpec',struct('type','identity'), ...
    'pi',piVec(:)');
if strcmpi(regime,'parallel')
    limitName='parallel'; model.distName='exponential'; model.durationCV=1;
else
    limitName='frontier';
end
out=analytical_response(limitName,theta,0:numel(piVec)-1,model);
t90=milestone_time(theta,out.systemRestoredFraction,0.9,'linear');
end
