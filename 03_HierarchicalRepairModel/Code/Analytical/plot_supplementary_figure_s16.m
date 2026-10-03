function files=plot_supplementary_figure_s16(root,mode,dataRoot,outputRoot)
if nargin<3 || isempty(dataRoot), dataRoot=resolve_figure_data_root(root,mode); end
if nargin<4 || isempty(outputRoot), outputRoot=resolve_figure_output_root(root,mode); end
sourceFile=fullfile(dataRoot,'Analytical','Analytical_results.mat');
if ~isfile(sourceFile), error('Missing analytical results. Run run_hierarchical_repair_analysis first.'); end
loaded=load(sourceFile,'s16','cfg'); results=loaded.s16; cfg=loaded.cfg;
cfg.output.supp=fullfile(outputRoot,'supplementary');
cfg.figure.fontSize=10;
cfg.figure.analytical03HeightCm=16.5;
if ~isfolder(cfg.output.supp), mkdir(cfg.output.supp); end

families = {'dependency','rates','release_duration'};
familyTitles = {'dependency structures','layer-specific rates', ...
    'service release/duration'};
letters = 'abcdefghijkl';
outputNames = {'Supplementary_Figure_S16a_dependency_structures', ...
    'Supplementary_Figure_S16b_layer_specific_rates', ...
    'Supplementary_Figure_S16c_service_release_duration'};
files=cell(3,1);

for r = 1:3
    fig = figure('Color','w','Units','centimeters', ...
        'Position',[1 1 cfg.figure.widthCm cfg.figure.analytical03HeightCm]);
    tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
    frontier = results{r,1};
    parallel = results{r,2};
    base = 4*(r-1);

    axFrontier = nexttile(tl,1);
    plot_curve_family(axFrontier,frontier,cfg, ...
        letters(base+1),['Frontier repair: ' familyTitles{r}]);

    axParallel = nexttile(tl,2);
    plot_curve_family(axParallel,parallel,cfg, ...
        letters(base+2),['Parallel repair: ' familyTitles{r}]);

    axValues = nexttile(tl,3);
    axis(axValues,'off');

    axClasses = nexttile(tl,4);
    axis(axClasses,'off');
    drawnow;
    plot_split_values(fig,axValues,frontier,parallel,cfg,letters(base+3));
    plot_split_classes(fig,axClasses,frontier,parallel,cfg,letters(base+4));

    stem=fullfile(cfg.output.supp,outputNames{r});
    export_publication_figure(fig,stem,cfg.figure.dpi); files{r}=[stem '.fig'];
    close(fig);
end

end

function plot_curve_family(ax,res,cfg,panelLetter,titleText)
hold(ax,'on');
cols = lines(numel(res.names));
for k = 1:numel(res.names)
    curve = res.curves{k};
    fit = res.fits{k};
    keep = curve.u <= cfg.fit.uMax;
    label = sprintf('%s (MAE = %.3f)',res.names{k},fit.mae);
    plot(ax,curve.u(keep),curve.R(keep),'-','LineWidth',1.35, ...
        'Color',cols(k,:),'DisplayName',label);
end
xlim(ax,[0 cfg.fit.uMax]);
ylim(ax,[0 1]);
xlabel(ax,'Normalized recovery time, $\tau/\tau_{50}$', ...
    'Interpreter','latex', ...
    'FontSize',cfg.figure.fontSize);
ylabel(ax,'Restored fraction','FontName',cfg.figure.fontName, ...
    'FontSize',cfg.figure.fontSize);
legend(ax,'Location','southeast','Box','off','FontSize',8);
apply_axes_style(ax,cfg);
titleSize=cfg.figure.fontSize;
ax.TitleFontSizeMultiplier=titleSize/ax.FontSize;
title(ax,titleText,'FontWeight','normal','Interpreter','tex');
text(ax,-0.12,1.08,panelLetter,'Units','normalized', ...
    'FontName',cfg.figure.fontName,'FontSize',11,'FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'Clipping','off');
end

function plot_split_values(fig,container,frontier,parallel,cfg,panelLetter)
[axFrontier,axParallel] = split_container_axes(fig,container);
scaledFrontier = zeros(size(frontier.values));
scaledParallel = zeros(size(parallel.values));
for m = 1:5
    block = [frontier.values(m,:) parallel.values(m,:)];
    lo = min(block(:));
    hi = max(block(:));
    if hi > lo
        scaledFrontier(m,:) = (frontier.values(m,:)-lo)/(hi-lo);
        scaledParallel(m,:) = (parallel.values(m,:)-lo)/(hi-lo);
    else
        scaledFrontier(m,:) = 0.5;
        scaledParallel(m,:) = 0.5;
    end
end
plot_value_heatmap(axFrontier,frontier.values,scaledFrontier, ...
    frontier.names,cfg,false,'Indicator values: frontier');
plot_value_heatmap(axParallel,parallel.values,scaledParallel, ...
    parallel.names,cfg,true,'Parallel');
text(axFrontier,-0.15,1.30,panelLetter,'Units','normalized', ...
    'FontName',cfg.figure.fontName,'FontSize',11,'FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'Clipping','off');
end

function plot_value_heatmap(ax,values,scaled,names,cfg,showX,titleText)
imagesc(ax,scaled,[0 1]);
colormap(ax,parula);
format_heatmap_axes(ax,names,cfg,showX);
for row = 1:5
    for col = 1:numel(names)
        textColor = contrast_color(scaled(row,col));
        text(ax,col,row,sprintf('%.2f',values(row,col)), ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontSize',5.6,'Color',textColor);
    end
end
set_heatmap_title(ax,titleText,cfg,showX);
end

function plot_split_classes(fig,container,frontier,parallel,cfg,panelLetter)
[axFrontier,axParallel] = split_container_axes(fig,container);
plot_class_heatmap(axFrontier,frontier.classes,frontier.names,cfg,false, ...
    'Empirical classification: frontier');
plot_class_heatmap(axParallel,parallel.classes,parallel.names,cfg,true, ...
    'Parallel');
text(axFrontier,-0.15,1.30,panelLetter,'Units','normalized', ...
    'FontName',cfg.figure.fontName,'FontSize',11,'FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'Clipping','off');
add_shared_class_colorbar(fig,axFrontier,axParallel,cfg);
end

function plot_class_heatmap(ax,classes,names,cfg,showX,titleText)
displayClasses=classes;
displayClasses(classes==2)=1;
displayClasses(classes==3)=2;
imagesc(ax,displayClasses,[0.5 2.5]);
colormap(ax,cfg.colors.classMap([1 3],:));
format_heatmap_axes(ax,names,cfg,showX);
set_heatmap_title(ax,titleText,cfg,showX);
end

function add_shared_class_colorbar(fig,axTop,axBottom,cfg)
drawnow;
topPosition = axTop.Position;
bottomPosition = axBottom.Position;
colorbarGap = 0.010;
colorbarWidth = 0.016;
rightMargin = 0.060;

topPosition(3) = topPosition(3)-rightMargin;
bottomPosition(3) = bottomPosition(3)-rightMargin;
axTop.Position = topPosition;
axBottom.Position = bottomPosition;

legendAxes = axes(fig,'Position',[topPosition(1) bottomPosition(2) ...
    0.001 0.001],'Visible','off','HandleVisibility','off');
colormap(legendAxes,cfg.colors.classMap([1 3],:));
caxis(legendAxes,[0.5 2.5]);
cb = colorbar(legendAxes);
cb.Position = [topPosition(1)+topPosition(3)+colorbarGap ...
    bottomPosition(2) colorbarWidth ...
    topPosition(2)+topPosition(4)-bottomPosition(2)];
cb.Ticks = 1:2;
cb.TickLabels = {'In 5-95','Out'};
cb.Ruler.TickLabelRotation = 0;
cb.FontSize = 10;
end

function set_heatmap_title(ax,titleText,cfg,~)
titleSize = cfg.figure.fontSize;
ax.TitleFontSizeMultiplier=titleSize/ax.FontSize;
title(ax,titleText,'FontWeight','normal','Interpreter','tex');
end

function [axTop,axBottom] = split_container_axes(fig,container)
position = container.Position;
delete(container);
interRowClearance = 0.035;
bottomLabelClearance = 0.100;
position(2) = position(2)+bottomLabelClearance;
position(4) = max(position(4)-interRowClearance-bottomLabelClearance,0.1);
gap = 0.20*position(4);
subHeight = 0.5*(position(4)-gap);
axBottom = axes(fig,'Position',[position(1) position(2) ...
    position(3) subHeight]);
axTop = axes(fig,'Position',[position(1) position(2)+subHeight+gap ...
    position(3) subHeight]);
end

function format_heatmap_axes(ax,names,cfg,showX)
metricLabels = cfg.sensitivity.indicatorLabels;
rowLabels = cell(5,1);
for m = 1:5
    rowLabels{m} = ['$' metricLabels{m} '$'];
end
xticks(ax,1:numel(names));
if showX
    xticklabels(ax,names);
    xtickangle(ax,35);
else
    xticklabels(ax,[]);
end
yticks(ax,1:5);
yticklabels(ax,rowLabels);
set(ax,'TickLabelInterpreter','latex');
axis(ax,'tight');
apply_axes_style(ax,cfg);
set(ax,'FontSize',10);
end

function color = contrast_color(value)
map = parula(256);
index = 1+round(min(max(value,0),1)*(size(map,1)-1));
background = map(index,:);
luminance = 0.2126*background(1)+0.7152*background(2)+ ...
    0.0722*background(3);
if luminance >= 0.52
    color = [0 0 0];
else
    color = [1 1 1];
end
end
