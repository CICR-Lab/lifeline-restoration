function file = figure4g_capacity_priority_matlab(rootDir, mode, dataRoot, outputRoot)

if nargin < 1
    rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
if nargin < 2
    mode = 'full';
end

if nargin < 3 || isempty(dataRoot), dataRoot = resolve_figure_data_root(rootDir,mode); end
if nargin < 4 || isempty(outputRoot), outputRoot = resolve_figure_output_root(rootDir,mode); end
sourceFile = fullfile(dataRoot, 'Simulation', 'source_data', ...
    'Fig4g_capacity_priority.csv');
outDir = fullfile(outputRoot, 'main');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

t = readtable(sourceFile);
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 8.0 7.6]);
ax = axes(fig, 'Position', [.24 .22 .71 .66]);

x = [.01 .03 .1 .3 1];
q = [0 .5 1];
colors = [
    .906 .435 .318
    .424 .353 .682
    .173 .498 .722];
labels = {'$q_F=0$, unstructured', '$q_F=0.5$, mixed', ...
    '$q_F=1$, frontier'};
lines = gobjects(3, 1);

hold(ax, 'on');
for i = 1:3
    meanValue = nan(1, 5);
    lower = meanValue;
    upper = meanValue;
    for j = 1:5
        row = abs(t.responseAxis2Value - q(i)) < 1e-10 & ...
            abs(t.responseAxis1Value - x(j)) < 1e-10;
        meanValue(j) = 100 * t.mean(row);
        lower(j) = 100 * t.p05(row);
        upper(j) = 100 * t.p95(row);
    end

    fill(ax, [x fliplr(x)], [lower fliplr(upper)], colors(i, :), ...
        'FaceAlpha', .12, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
    lines(i) = plot(ax, x, meanValue, '-o', ...
        'Color', colors(i, :), 'LineWidth', 1.6, ...
        'MarkerSize', 4.2, 'MarkerFaceColor', colors(i, :), ...
        'MarkerEdgeColor', 'w');
end
hold(ax, 'off');

set(ax, 'XScale', 'log', 'XLim', [.0085 1.18], 'YLim', [0 100], ...
    'XTick', x, 'XTickLabel', {'0.01', '0.03', '0.10', '0.30', '1'}, ...
    'FontName', 'Arial', 'FontSize', 9, ...
    'TickLabelInterpreter', 'none', 'LineWidth', .8, ...
    'Box', 'off', 'TickDir', 'out');
grid(ax, 'on');
ax.GridColor = [.85 .87 .88];
ax.GridAlpha = .8;

mixed_font_xlabel(ax, 'Relative repair capacity,', '$C/N$', 'Arial', 9);
ylabel(ax, 'Mean marginal consistency (%)', ...
    'Interpreter', 'none', 'FontName', 'Arial', 'FontSize', 9);

lgd = legend(ax, lines, labels, 'Interpreter', 'latex', ...
    'Location', 'southwest', 'Box', 'off');
lgd.FontSize = 7.0;
stem=fullfile(outDir,'Figure4g_capacity_priority_matlab');
export_publication_figure(fig,stem,450);
file=[stem '.fig'];
close(fig);
end

