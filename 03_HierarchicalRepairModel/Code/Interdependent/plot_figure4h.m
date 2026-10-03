function plot_figure4h(catalog, results, coupling, cfg)

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 8.1 8.0]);
ax = axes(fig, 'Position', [.14 .15 .59 .66]);
hold(ax, 'on');

[~, h2] = contour(ax, coupling.GAMMA, coupling.GB, ...
    coupling.limit(2).epsilonUB, coupling.levels, ...
    'LineWidth', 1.1, 'LineStyle', '--', 'ShowText', 'off');
[~, h1] = contour(ax, coupling.GAMMA, coupling.GB, ...
    coupling.limit(1).epsilonUB, coupling.levels, ...
    'LineWidth', 1.1, 'LineStyle', '-', 'ShowText', 'off');
h1.LineColor = [.15 .15 .15];
h2.LineColor = [.15 .15 .15];

idx = find(strcmp({catalog.group}, 'core grid'));
x = nan(numel(idx), 1);
y = x;
z = x;
for k = 1:numel(idx)
    output = results{idx(k)};
    x(k) = median(output.metrics.gamma, 'omitnan');
    y(k) = median(output.metrics.gB, 'omitnan');
    z(k) = median(output.metrics.epsilonDec, 'omitnan');
end

sc = scatter(ax, x, y, 42, z, 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', .6);
colormap(ax, parula(256));
caxis(ax, [0 max([cfg.coupled.mainContourLevels(:); z(:); .01])]);

set(ax, 'XScale', 'log', 'XTick', [.25 .5 1 2 4 8], ...
    'XTickLabel', {'0.25', '0.5', '1', '2', '4', '8'}, ...
    'FontName', 'Arial', 'FontSize', 9, ...
    'TickLabelInterpreter', 'none', 'LineWidth', .8, 'Layer', 'top');
xlim(ax, [.2 cfg.coupled.gammaPlotMax]);
ylim(ax, [0 .9]);
xlabel(ax, '$\gamma$', ...
    'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 9);
ylabel(ax, 'g_B', ...
    'Interpreter', 'tex', 'FontName', 'Arial', 'FontSize', 9);
grid(ax, 'on');
box(ax, 'on');

cb = colorbar(ax);
cb.FontName = 'Arial';
cb.FontSize = 9;
cb.Label.String = 'Relative error, ε_{dec}^{sim}';
cb.Label.Interpreter = 'tex';
cb.Label.FontName = 'Arial';
cb.Label.FontSize = 9;

lgd = legend(ax, [h1 h2 sc], ...
    {'Frontier-limited', 'Fully parallel', 'Simulation medians'}, ...
    'Orientation', 'horizontal', 'NumColumns', 3, ...
    'FontName', 'Arial', 'FontSize', 6.5, ...
    'Interpreter', 'none', 'Box', 'off');
lgd.ItemTokenSize = [12 8];
lgd.AutoUpdate = 'off';

drawnow;
ax.Position = [.14 .15 .59 .66];
cb.Position = [.77 .15 .035 .66];
lgd.Units = 'normalized';
lgd.Position = [.11 .86 .86 .06];
export_publication_figure(fig, ...
    fullfile(cfg.outputDir, 'Figure4h_coupled_network_contours'), ...
    cfg.figureResolution);
close(fig);
end
