function file=plot_supplementary_figure_s15(root,mode,dataRoot,outputRoot)
if nargin<3 || isempty(dataRoot), dataRoot=resolve_figure_data_root(root,mode); end
if nargin<4 || isempty(outputRoot), outputRoot=resolve_figure_output_root(root,mode); end
sourceFile=fullfile(dataRoot,'Analytical','Analytical_results.mat');
if ~isfile(sourceFile), error('Missing analytical results. Run run_hierarchical_repair_analysis first.'); end
loaded=load(sourceFile,'s15','cfg'); grid=loaded.s15; cfg=loaded.cfg;
cfg.output.supp=fullfile(outputRoot,'supplementary');
if ~isfolder(cfg.output.supp), mkdir(cfg.output.supp); end

fig = figure('Color','w','Units','centimeters', ...
    'Position',[1 1 24 20.5]);
tl = tiledlayout(fig,5,17,'TileSpacing','tight','Padding','tight');
letters = 'abcdefghij';

for m = 1:numel(cfg.sensitivity.indicatorNames)
    valuePair = squeeze(grid.valuesByLimit(:,:,m,:));
    finiteValues = valuePair(isfinite(valuePair));
    sharedLimits = [min(finiteValues) max(finiteValues)];
    if sharedLimits(1) == sharedLimits(2)
        sharedLimits = sharedLimits + [-0.5 0.5];
    end

    rowStart = 17*(m-1);
    spacer = nexttile(tl,rowStart+9);
    axis(spacer,'off');

    for r = 1:numel(grid.limitNames)
        if r == 1
            valueTile = rowStart+1;
            classTile = rowStart+5;
        else
            valueTile = rowStart+10;
            classTile = rowStart+14;
        end
        valueMap = grid.valuesByLimit(:,:,m,r);
        classMap = grid.classesByLimit(:,:,m,r);

        axValue = nexttile(tl,valueTile,[1 4]);
        plot_value_map(axValue,valueMap,sharedLimits,grid,cfg,m,r, ...
            letters(2*(m-1)+r));

        axClass = nexttile(tl,classTile,[1 4]);
        plot_class_map(axClass,classMap,grid,cfg,m,r);
    end
end

stem=fullfile(cfg.output.supp,'Supplementary_Figure_S15');
export_publication_figure(fig,stem,cfg.figure.dpi); file=[stem '.fig'];
close(fig);
end

function plot_value_map(ax,values,colorLimits,grid,cfg,metricIndex, ...
        limitIndex,panelLetter)
imagesc(ax,grid.depthNorm,grid.L,values);
axis(ax,'xy');
colormap(ax,parula);
caxis(ax,colorLimits);
colorbar(ax);
format_map_axes(ax,grid,cfg,metricIndex,true);
title(ax,panel_title(grid,cfg,metricIndex,limitIndex,'Model value', ...
    panelLetter), ...
    'FontWeight','normal','Interpreter','latex');
end

function plot_class_map(ax,classes,grid,cfg,metricIndex,limitIndex)
displayClasses=classes;
displayClasses(classes==2)=1;
displayClasses(classes==3)=2;
imagesc(ax,grid.depthNorm,grid.L,displayClasses);
axis(ax,'xy');
colormap(ax,cfg.colors.classMap([1 3],:));
caxis(ax,[0.5 2.5]);
cb = colorbar(ax);
cb.Ticks = 1:2;
cb.TickLabels = {'In 5-95','Out'};
cb.FontSize = cfg.figure.fontSize - 2;
format_map_axes(ax,grid,cfg,metricIndex,false);
title(ax,panel_title(grid,cfg,metricIndex,limitIndex,'Empirical class',''), ...
    'FontWeight','normal','Interpreter','latex');
if all(~isfinite(classes(:)))
    text(ax,0.5,0.5,'Empirical ranges not supplied','Units','normalized', ...
        'HorizontalAlignment','center','FontSize',7,'Color',cfg.colors.gray);
end
end

function format_map_axes(ax,grid,cfg,metricIndex,showYLabels)
if showYLabels
    ylabel(ax,'Maximum depth, $L$','Interpreter','latex');
else
    ax.YTickLabel = [];
end
if metricIndex == numel(cfg.sensitivity.indicatorNames)
    xlabel(ax,'$\bar{l}_{\pi}/L$','Interpreter','latex');
else
    ax.XTickLabel = [];
end
xlim(ax,[grid.depthNorm(1) grid.depthNorm(end)]);
ylim(ax,[grid.L(1)-0.5 grid.L(end)+0.5]);
apply_axes_style(ax,cfg);
set(ax,'TickLabelInterpreter','latex');
end

function titleText = panel_title(grid,cfg,metricIndex,limitIndex,mapType, ...
        panelLetter)
metric = cfg.sensitivity.indicatorLabels{metricIndex};
if isempty(panelLetter)
    metricLine = sprintf('$%s$: %s',metric,mapType);
else
    metricLine = sprintf('$\\mathbf{%s}\\quad %s$: %s', ...
        panelLetter,metric,mapType);
end
if metricIndex == 1
    titleText = {grid.limitLabels{limitIndex},metricLine};
else
    titleText = metricLine;
end
end
