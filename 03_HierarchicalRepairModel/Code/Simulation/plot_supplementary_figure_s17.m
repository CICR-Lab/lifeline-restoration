function file=plot_supplementary_figure_s17(rootDir,mode,dataRoot,outputRoot)

if nargin<3 || isempty(dataRoot), dataRoot=resolve_figure_data_root(rootDir,mode); end
if nargin<4 || isempty(outputRoot), outputRoot=resolve_figure_output_root(rootDir,mode); end
sourceDir=fullfile(dataRoot,'Simulation','source_data');
outputDir=fullfile(outputRoot,'supplementary');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end
fontSize = 9;

consistency=readtable(fullfile(sourceDir,'FigS17a_consistency.csv'));
allFive=readtable(fullfile(sourceDir,'FigS17a_allfive.csv'));
audit=readtable(fullfile(sourceDir,'FigS17b_model_data_audit.csv'));
variance=readtable(fullfile(sourceDir,'FigS17c_variance_structure.csv'));
sobol=readtable(fullfile(sourceDir,'FigS17d_group_sobol.csv'));

fig = figure('Color','w','Units','centimeters', ...
    'Position',[2 2 23.9 17.4]);
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

draw_a(nexttile(layout,1),consistency,allFive,fontSize);
draw_b(nexttile(layout,2),audit,fontSize);
draw_c(nexttile(layout,3),variance,fontSize);
draw_d(nexttile(layout,4),sobol,fontSize);

stem=fullfile(outputDir,'Supplementary_Figure_S17');
export_publication_figure(fig,stem,300);
file=[stem '.fig'];
close(fig);
end

function draw_a(ax,t,allFive,fontSize)
metricNames = metrics();
metricLabels = labels();
names = cellstr(t.metricName);
values = nan(1,5);
for i = 1:5
    values(i) = 100*t.estimate(strcmp(names,metricNames{i}));
end
meanMarginal = 100*t.estimate(strcmp(names,'mean_five'));
jointValue = 100*allFive.estimate(1);

hold(ax,'on');
for i = 1:5
    plot(ax,[0 values(i)],[i i],'-','Color',[.86 .92 .96], ...
        'LineWidth',5);
    plot(ax,values(i),i,'o','MarkerSize',6, ...
        'MarkerFaceColor',[.20 .50 .73],'MarkerEdgeColor','w');
    text(ax,values(i)+1.4,i,sprintf('%.1f',values(i)), ...
        'FontName','Arial','FontSize',fontSize,'VerticalAlignment','middle');
end
plot(ax,meanMarginal,7,'d','MarkerSize',7, ...
    'MarkerFaceColor',[.12 .13 .14]);
plot(ax,jointValue,8,'s','MarkerSize',7,'MarkerFaceColor','w', ...
    'LineWidth',1.2);
text(ax,meanMarginal+1.4,7,sprintf('%.1f',meanMarginal), ...
    'FontName','Arial','FontSize',fontSize,'VerticalAlignment','middle');
text(ax,jointValue+1.4,8,sprintf('%.1f',jointValue), ...
    'FontName','Arial','FontSize',fontSize,'VerticalAlignment','middle');
plot(ax,[0 100],[6.25 6.25],'-','Color',[.75 .77 .78]);
hold(ax,'off');

set(ax,'XLim',[0 100],'YLim',[.6 8.4],'YDir','reverse', ...
    'YTick',[1:5 7 8], ...
    'YTickLabel',[metricLabels {'Mean marginal','All five simultaneously'}], ...
    'TickLabelInterpreter','latex');
style_axes(ax,fontSize);
xlabel(ax,'Fraction of base-ensemble trajectories (%)');
text(ax,.98,.03,'n = 409,600','Units','normalized', ...
    'HorizontalAlignment','right', ...
    'FontName','Arial','FontSize',fontSize);
panel_label_s17(ax,'a',fontSize+1);
end

function draw_b(ax,t,fontSize)
metricNames = metrics();
metricLabels = labels();
names = cellstr(t.metricName);
readouts = cellstr(t.readout);

hold(ax,'on');
for i = 1:5
    row = strcmp(names,metricNames{i}) & strcmp(readouts,'primary');
    current = t(row,:);
    scale = current.empiricalMedian;
    empirical = [current.empiricalP05 current.empiricalMedian ...
        current.empiricalP95]/scale;
    simulation = [current.simulationP05 current.simulationMedian ...
        current.simulationP95]/scale;
    plot(ax,empirical([1 3]),[i-.12 i-.12],'-', ...
        'Color',[.82 .82 .83],'LineWidth',9);
    plot(ax,empirical(2),i-.12,'o','MarkerSize',5, ...
        'MarkerFaceColor',[.16 .16 .16]);
    plot(ax,simulation([1 3]),[i+.12 i+.12],'-', ...
        'Color',[.20 .50 .73],'LineWidth',1.8);
    plot(ax,simulation(2),i+.12,'o','MarkerSize',5, ...
        'MarkerFaceColor',[.20 .50 .73],'MarkerEdgeColor','w');
end
empiricalLegend = plot(ax,nan,nan,'-','Color',[.82 .82 .83], ...
    'LineWidth',9);
simulationLegend = plot(ax,nan,nan,'-','Color',[.20 .50 .73], ...
    'LineWidth',1.8);
plot(ax,[1 1],[.45 5.55],'--','Color',[.7 .72 .73]);
hold(ax,'off');

set(ax,'XScale','log','XLim',[.005 11], ...
    'XTick',[.005 .01 .03 .1 .3 1 3 10], ...
    'YLim',[-.1 5.55],'YDir','reverse','YTick',1:5, ...
    'YTickLabel',metricLabels,'TickLabelInterpreter','latex');
style_axes(ax,fontSize);
xlabel(ax,'Indicator value / empirical median (log scale)');
lgd = legend(ax,[empiricalLegend simulationLegend], ...
    {'Empirical 5th--95th percentile', ...
    'Simulation 5th--95th percentile'}, ...
    'Location','northwest','Box','off');
lgd.FontName = 'Arial';
lgd.FontSize = fontSize;
lgd.ItemTokenSize = [22 8];
panel_label_s17(ax,'b',fontSize+1);
end

function draw_c(ax,t,fontSize)
outcomeNames = outcomes();
metricLabels = labels();
names = cellstr(t.outcomeName);
designShare = nan(1,5);
stochasticShare = nan(1,5);
lower = nan(1,5);
upper = nan(1,5);
hasCI = all(ismember({'designShareP025','designShareP975'}, ...
    t.Properties.VariableNames));

for i = 1:5
    row = strcmp(names,outcomeNames{i});
    designShare(i) = 100*t.designShare(row);
    stochasticShare(i) = 100*t.stochasticShare(row);
    if hasCI
        lower(i) = 100*t.designShareP025(row);
        upper(i) = 100*t.designShareP975(row);
    end
end

bars = bar(ax,1:5,[designShare' stochasticShare'],'stacked', ...
    'BarWidth',.68);
bars(1).FaceColor = [.20 .50 .73];
bars(2).FaceColor = [.79 .79 .81];
hold(ax,'on');
if hasCI
    errorbar(ax,1:5,designShare,designShare-lower,upper-designShare, ...
        'k.','LineWidth',.9,'CapSize',5);
end
for i = 1:5
    text(ax,i,designShare(i)/2,sprintf('%.1f',designShare(i)), ...
        'Color','w','FontWeight','bold','HorizontalAlignment','center', ...
        'FontName','Arial','FontSize',fontSize);
end
hold(ax,'off');

set(ax,'XTick',1:5,'XTickLabel',metricLabels,'YLim',[0 115], ...
    'YTick',0:20:100, ...
    'TickLabelInterpreter','latex');
style_axes(ax,fontSize);
ylabel(ax,'Share of transformed run-level variance (%)');
lgd = legend(ax,bars,{'Input-driven','Stochastic realization'}, ...
    'Location','northoutside','Orientation','horizontal','Box','off');
lgd.FontName = 'Arial';
lgd.FontSize = fontSize;
panel_label_s17(ax,'c',fontSize+1);
end

function draw_d(ax,t,fontSize)
outcomeNames = outcomes();
metricLabels = labels();
names = cellstr(t.outcomeName);
groupShare = nan(5,3);
for i = 1:5
    for group = 1:3
        row = t.groupID == group & strcmp(names,outcomeNames{i});
        groupShare(i,group) = 100*t.firstOrderRaw(row);
    end
end
values = [groupShare 100-sum(groupShare,2)];
colors = [.31 .50 .68; .35 .63 .31; .90 .56 .21; .45 .46 .48];

bars = bar(ax,1:5,values,'stacked','BarWidth',.68);
for group = 1:4
    bars(group).FaceColor = colors(group,:);
end
hold(ax,'on');
base = zeros(5,1);
for group = 1:4
    for i = 1:5
        if values(i,group) >= 8
            text(ax,i,base(i)+values(i,group)/2, ...
                sprintf('%.1f',values(i,group)), ...
                'HorizontalAlignment','center','Color','w', ...
                'FontWeight','bold','FontName','Arial','FontSize',fontSize);
        end
    end
    base = base+values(:,group);
end
hold(ax,'off');

upperLimit=max(105,5*ceil(max(sum(max(values,0),2))/5));
set(ax,'XTick',1:5,'XTickLabel',metricLabels,'YLim',[0 upperLimit], ...
    'YTick',0:20:100, ...
    'TickLabelInterpreter','latex');
style_axes(ax,fontSize);
ylabel(ax,'Share of input-driven variance (%)');
lgd = legend(ax,bars,{'Task-network architecture', ...
    'Task and service attributes','Repair organization', ...
    'Cross-group interactions'}, ...
    'Location','northoutside','NumColumns',2,'Box','off');
lgd.FontName = 'Arial';
lgd.FontSize = fontSize;
panel_label_s17(ax,'d',fontSize+1);
end

function style_axes(ax,fontSize)
set(ax,'FontName','Arial','FontSize',fontSize,'LineWidth',.8, ...
    'Box','off','TickDir','out','Layer','top');
grid(ax,'off');
ax.TitleFontSizeMultiplier = 1;
end

function panel_label_s17(ax,labelText,fontSize)
text(ax,-0.12,1.11,labelText,'Units','normalized','FontWeight','bold', ...
    'FontName','Arial','FontSize',fontSize,'HorizontalAlignment','left', ...
    'VerticalAlignment','top','Clipping','off');
end

function names = metrics()
names = {'tau80_tau50','tau90_tau80','tau95_tau90','kappa','eta90'};
end

function names = outcomes()
names = {'log_tau80_tau50','log_tau90_tau80','log_tau95_tau90', ...
    'log_kappa','logit_eta90'};
end

function textLabels = labels()
textLabels = {'$\tau_{80}/\tau_{50}$','$\tau_{90}/\tau_{80}$', ...
    '$\tau_{95}/\tau_{90}$','$\kappa$','$\eta_{90}$'};
end

