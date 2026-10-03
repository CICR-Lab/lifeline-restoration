function files=plot_analytical_figure4(root,mode,dataRoot,outputRoot)
if nargin<3 || isempty(dataRoot), dataRoot=resolve_figure_data_root(root,mode); end
if nargin<4 || isempty(outputRoot), outputRoot=resolve_figure_output_root(root,mode); end
sourceFile=fullfile(dataRoot,'Analytical','Analytical_results.mat');
if ~isfile(sourceFile)
    error('Missing analytical source data. Run run_hierarchical_repair_analysis first.');
end
s=load(sourceFile,'b','c','d','cfg');
cfg=s.cfg; cfg.output.main=fullfile(outputRoot,'main');
if ~isfolder(cfg.output.main), mkdir(cfg.output.main); end
files.b=plot_b(s.b,cfg);
files.c=plot_c(s.c,cfg);
files.d=plot_d(s.d,cfg);
end

function file=plot_b(data,cfg)
[fig,ax]=panel_axes(cfg); hold(ax,'on');
cols=[cfg.colors.blue;cfg.colors.orange;cfg.colors.green;cfg.colors.red];
for k=1:numel(data.depths)
    plot(ax,data.uFrontier(:,k),data.frontier(:,k),'-', ...
        'LineWidth',cfg.figure.lineWidth,'Color',cols(k,:), ...
        'DisplayName',sprintf('$l=%d$ (F)',data.depths(k)));
end
for k=1:numel(data.depths)
    plot(ax,data.uParallel(:,k),data.parallel(:,k),'--', ...
        'LineWidth',1.4,'Color',cols(k,:), ...
        'DisplayName',sprintf('$l=%d$ (P)',data.depths(k)));
end
xlim(ax,[0 7]); ylim(ax,[0 1]);
mixed_font_xlabel(ax,'Normalized recovery time,','$\tau/\tau_{50}$', ...
    cfg.figure.fontName,cfg.figure.fontSize);
ylabel(ax,'Layer-associated restored fraction','Interpreter','tex');
lgd=legend(ax,'Location','southeast','Box','off','NumColumns',1, ...
    'Interpreter','latex','FontSize',7.5);
lgd.Title.String='';
drawnow;
lgd.Units='normalized';
legendPosition=lgd.Position;
annotation(fig,'textbox',[legendPosition(1), ...
    legendPosition(2)+legendPosition(4)+0.008,legendPosition(3),0.050], ...
    'String',{'F: frontier-limited';'P: parallel repair'}, ...
    'EdgeColor','none','Margin',0,'FontName',cfg.figure.fontName, ...
    'FontSize',7.5,'Interpreter','latex','HorizontalAlignment','left', ...
    'VerticalAlignment','bottom');
apply_axes_style(ax,cfg);
file=finish(fig,cfg,'Figure4b_depth_curves');
end

function file=plot_c(data,cfg)
[fig,ax]=panel_axes(cfg); hold(ax,'on');
cols=[cfg.colors.blue;cfg.colors.orange;cfg.colors.green];
markers={'o','s','d'};
metricNames={'tau80_tau50','tau90_tau80','tau95_tau90'};
metricLabels={'\tau_{80}/\tau_{50}','\tau_{90}/\tau_{80}', ...
    '\tau_{95}/\tau_{90}'};
xMin=-.3; pooledX=8.5; xMax=9;
bracketX=[7.78 8.18 8.58];
bracketNotch=[0.50 0.78 0.50];
for m=1:3
    range=get_empirical_range(data.empirical,metricNames{m});
    patch(ax,[xMin xMax xMax xMin],[range.centralLow range.centralLow ...
        range.centralHigh range.centralHigh],cols(m,:), ...
        'FaceAlpha',.10,'EdgeColor','none','HandleVisibility','off');
    plot(ax,data.depths,data.frontier(:,m),['-' markers{m}], ...
        'Color',cols(m,:),'LineWidth',1.4,'MarkerSize',4, ...
        'DisplayName',sprintf('$%s$ (F)',metricLabels{m}));
    plot(ax,data.depths,data.parallel(:,m),['--' markers{m}], ...
        'Color',cols(m,:),'LineWidth',1.2,'MarkerSize',4, ...
        'DisplayName',sprintf('$%s$ (P)',metricLabels{m}));
    draw_ratio_bracket(ax,bracketX(m),range.centralLow,range.centralHigh, ...
        cols(m,:),metricLabels{m},cfg.figure.fontSize-2,bracketNotch(m));
end
xlim(ax,[xMin xMax]); ylim(ax,[1 2.45]); xticks(ax,data.depths);
apply_axes_style(ax,cfg);
mixed_font_xlabel(ax,'Initial dependency depth,','$l$', ...
    cfg.figure.fontName,cfg.figure.fontSize);
ylabel(ax,'Milestone ratio','Interpreter','tex');
lgd=legend(ax,'Location','northeast','Box','off','NumColumns',2, ...
    'Interpreter','latex','FontSize',7.5);
drawnow; lgd.Units='normalized'; posL=lgd.Position;
lgd.Position=[posL(1),posL(2)+.025,posL(3),posL(4)];
file=finish(fig,cfg,'Figure4c_phase_ratios');
end

function draw_ratio_bracket(ax,x,y0,y1,colorValue,labelText,fontSize,notchFraction)
yNotch=y0+notchFraction*(y1-y0);
tick=.11;
waist=.065;
yNotch=min(max(yNotch,y0+1.3*waist),y1-1.3*waist);
plot(ax,[x-tick x],[y0 y0],'-','Color',colorValue, ...
    'LineWidth',1.2,'HandleVisibility','off','Clipping','off');
plot(ax,[x-tick x],[y1 y1],'-','Color',colorValue, ...
    'LineWidth',1.2,'HandleVisibility','off','Clipping','off');
plot(ax,[x x],[y0 yNotch-waist],'-','Color',colorValue, ...
    'LineWidth',1.2,'HandleVisibility','off','Clipping','off');
plot(ax,[x x],[yNotch+waist y1],'-','Color',colorValue, ...
    'LineWidth',1.2,'HandleVisibility','off','Clipping','off');
plot(ax,[x x+waist x],[yNotch-waist yNotch yNotch+waist],'-', ...
    'Color',colorValue,'LineWidth',1.2,'HandleVisibility','off', ...
    'Clipping','off');
end

function file=plot_d(data,cfg)
[fig,ax]=panel_axes(cfg); hold(ax,'on');
cols=[cfg.colors.blue;cfg.colors.orange;cfg.colors.green;cfg.colors.purple];
for k=1:numel(data.profiles)
    lo=min(data.frontier(:,k),data.parallel(:,k));
    hi=max(data.frontier(:,k),data.parallel(:,k));
    fill(ax,[data.uGrid;flipud(data.uGrid)],[lo;flipud(hi)],cols(k,:), ...
        'FaceAlpha',.10,'EdgeColor','none','HandleVisibility','off');
    plot(ax,data.uGrid,data.frontier(:,k),'-','LineWidth',1.6, ...
        'Color',cols(k,:),'DisplayName',sprintf('%s (F), $\\eta_{90}=%.3f$', ...
        data.profiles(k).displayName,data.metrics(k).frontier.eta90));
    plot(ax,data.uGrid,data.parallel(:,k),'--','LineWidth',1.6, ...
        'Color',cols(k,:),'DisplayName',sprintf('%s (P), $\\eta_{90}=%.3f$', ...
        data.profiles(k).displayName,data.metrics(k).parallel.eta90));
end
xlim(ax,[0 4.2]); ylim(ax,[0 1]);
mixed_font_xlabel(ax,'Normalized recovery time,','$\tau/\tau_{50}$', ...
    cfg.figure.fontName,cfg.figure.fontSize);
ylabel(ax,'System restored fraction','Interpreter','tex');
lgd=legend(ax,'Location','southeast','Box','off','Interpreter','latex', ...
    'NumColumns',1,'FontSize',5.8);
drawnow; lgd.Units='normalized'; posL=lgd.Position; lgd.Position=[posL(1)+.015,posL(2)-.020,posL(3),posL(4)];
apply_axes_style(ax,cfg); drawnow;
pos=ax.Position;
inset=axes(fig,'Position',[pos(1)+.57*pos(3),pos(2)+.52*pos(4), ...
    .30*pos(3),.24*pos(4)]);
B=vertcat(data.profiles.pi)'; bh=bar(inset,0:3,B,'grouped','BarWidth',.82);
for k=1:numel(bh), bh(k).FaceColor=cols(k,:); bh(k).EdgeColor='none'; end
xlim(inset,[-.5 3.5]); ylim(inset,[0 .7]); xticks(inset,0:3);
set(inset,'FontName',cfg.figure.fontName,'FontSize',7,'Box','off', ...
    'TickDir','out');
insetLabelSize=cfg.figure.fontSize-2;
mixed_font_xlabel(inset,'Initial dependency depth,','$l$', ...
    cfg.figure.fontName,insetLabelSize,0.020);
ylabel(inset,'$\pi_l$','Interpreter','latex','FontSize',insetLabelSize);
file=finish(fig,cfg,'Figure4d_service_distribution');
end

function [fig,ax]=panel_axes(cfg)
fig=figure('Color','w','Units','centimeters', ...
    'Position',[2 2 cfg.figure.mainPanelSizeCm]);
ax=axes(fig);
end

function file=finish(fig,cfg,name)
s=cfg.figure.mainPanelSizeCm;
set(fig,'PaperUnits','centimeters','PaperSize',s, ...
    'PaperPosition',[0 0 s],'PaperPositionMode','manual');
file=fullfile(cfg.output.main,[name '.fig']);
export_publication_figure(fig,fullfile(cfg.output.main,name),cfg.figure.dpi);
close(fig);
end
