function file = figure4f_task_service_matlab(rootDir, mode, dataRoot, outputRoot)

if nargin < 1
    rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
if nargin < 2
    mode = 'full';
end

if nargin < 3 || isempty(dataRoot), dataRoot = resolve_figure_data_root(rootDir,mode); end
if nargin < 4 || isempty(outputRoot), outputRoot = resolve_figure_output_root(rootDir,mode); end
sourceFile = fullfile(dataRoot, 'Simulation', 'source_data', ...
    'Fig4f_task_service.csv');
outDir = fullfile(outputRoot, 'main');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

t = readtable(sourceFile);
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 8.0 7.6]);
ax = axes(fig, 'Position', [.27 .50 .68 .38]);

tasks = {'shallow', 'uniform', 'deep'};
services = {'shallow', 'uniform', 'intermediate', 'deep', 'deep_tail'};
taskNames = cellstr(t.responseAxis1Label);
serviceNames = cellstr(t.responseAxis2Label);
M = nan(3, 5);

for i = 1:3
    for j = 1:5
        row = strcmp(taskNames, tasks{i}) & strcmp(serviceNames, services{j});
        M(i, j) = 100 * t.mean(row);
    end
end

imagesc(ax, M, [0 100]);
colormap(ax, viridis(256));

for i = 1:3
    for j = 1:5
        textColor = 'w';
        if M(i, j) >= 58
            textColor = [.06 .09 .12];
        end
        text(ax, j, i, sprintf('%.1f', M(i, j)), ...
            'HorizontalAlignment', 'center', ...
            'FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold', ...
            'Color', textColor);
    end
end

set(ax, 'XTick', 1:5, ...
    'XTickLabel', {'Shallow', 'Uniform', 'Intermediate', 'Deep', 'Deep-tail'}, ...
    'YTick', 1:3, 'YTickLabel', {'Shallow', 'Uniform', 'Deep'}, ...
    'TickLength', [0 0], 'FontName', 'Arial', 'FontSize', 9, ...
    'TickLabelInterpreter', 'none', 'Box', 'off');
xtickangle(ax, 35);

xlabel(ax, 'Layer service-share profile', ...
    'Interpreter', 'none', 'FontName', 'Arial', 'FontSize', 9);
ylabel(ax, 'Task-layer profile', ...
    'Interpreter', 'none', 'FontName', 'Arial', 'FontSize', 9);

cb = colorbar(ax, 'southoutside');
cb.Label.String = 'Mean marginal consistency (%)';
cb.Label.Interpreter = 'none';
cb.Label.FontName = 'Arial';
cb.Label.FontSize = 9;
cb.FontName = 'Arial';
cb.FontSize = 9;
cb.Ticks = 0:20:100;

ax.Position = [.27 .50 .68 .38];
cb.Position = [.27 .18 .68 .035];
stem=fullfile(outDir,'Figure4f_task_service_matlab');
export_publication_figure(fig,stem,450);
file=[stem '.fig'];
close(fig);
end

function map = viridis(n)
anchors = [
    .267 .005 .329
    .283 .141 .458
    .254 .265 .530
    .207 .372 .553
    .164 .471 .558
    .128 .567 .551
    .135 .659 .518
    .267 .749 .441
    .478 .821 .318
    .741 .873 .150
    .993 .906 .144];
map = interp1(linspace(0, 1, size(anchors, 1)), anchors, ...
    linspace(0, 1, n), 'pchip');
map = max(0, min(1, map));
end

