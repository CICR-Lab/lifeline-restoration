function file=plot_supplementary_figure_s18(root,mode,dataRoot,outputRoot)
if nargin<3 || isempty(dataRoot), dataRoot=resolve_figure_data_root(root,mode); end
if nargin<4 || isempty(outputRoot), outputRoot=resolve_figure_output_root(root,mode); end
sourceFile=fullfile(dataRoot,'Interdependent', ...
    'Interdependent_results.mat');
if ~isfile(sourceFile), error('Missing interdependent results. Run run_hierarchical_repair_analysis first.'); end
loaded=load(sourceFile,'catalog','results','cfg');
catalog=loaded.catalog; results=loaded.results; cfg=loaded.cfg;
cfg.figureDir=fullfile(outputRoot,'supplementary');
if ~isfolder(cfg.figureDir), mkdir(cfg.figureDir); end
fontSize = 9;
legendFontSize = 8;

fig = figure('Color','w','Units','centimeters', ...
    'Position',[2 2 23.9 17.4]);
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(layout,1);
draw_burden_rho90_contours( ...
    catalog,results,cfg,ax);
panel_label(ax,'a',fontSize);

ax = nexttile(layout,2);
hold(ax,'on');
[~,ref] = get_coupled_result_by_id(catalog,results,'core_gam1_gb0p5');
ex = ref.example;
[tFill,y0,yc] = common_curve_grid(ex.tB0,ex.RB0,ex.tBc,ex.RBc);
hFill = fill(ax,[tFill;flipud(tFill)],[y0;flipud(yc)], ...
    [0.85 0.85 0.85],'EdgeColor','none','FaceAlpha',0.45, ...
    'DisplayName','$\Delta L_{E\rightarrow B}^{\mathrm{sim}}$');
hE = stairs(ax,ex.tE,ex.RE,'-','LineWidth',1.7, ...
    'DisplayName','Electricity, $R_E^0$');
hB0 = stairs(ax,ex.tB0,ex.RB0,'--','LineWidth',1.7, ...
    'DisplayName','Dependent ungated, $R_B^0$');
hBc = stairs(ax,ex.tBc,ex.RBc,'-','LineWidth',1.7, ...
    'DisplayName','Dependent coupled, $R_B^c$');
xline(ax,ex.TEstar,':','$T_E^*$','LineWidth',1.2, ...
    'LabelVerticalAlignment','middle','Interpreter','latex', ...
    'HandleVisibility','off');
xlabel(ax,'Recovery-phase time, $\tau$', ...
    'Interpreter','latex','FontSize',fontSize);
ylabel(ax,'Restored service fraction', ...
    'FontName','Arial','FontSize',fontSize);
ylim(ax,[0 1]);
style_axes(ax);
legend(ax,[hE hB0 hBc hFill], ...
    {'Electricity, $R_E^0$','Ungated B, $R_B^0$', ...
    'Coupled B, $R_B^c$','$\Delta L_{E\rightarrow B}^{\mathrm{sim}}$'}, ...
    'Location','northwest','Interpreter','latex', ...
    'FontSize',legendFontSize,'Box','off');
panel_label(ax,'b',fontSize);

ax = nexttile(layout,3);
hold(ax,'on');
allUB = [];
allE = [];
medUB = [];
medE = [];
for i = 1:numel(catalog)
    output = results{i};
    if strcmp(catalog(i).id,'sens_alternative_paths')
        continue;
    end
    keep = output.metrics.boundComparable > 0.5;
    if ~any(keep)
        continue;
    end
    ub = output.metrics.epsilonUB(keep);
    realized = output.metrics.epsilonDec(keep);
    allUB = [allUB;ub];
    allE = [allE;realized];
    medUB(end+1,1) = median(ub,'omitnan');
    medE(end+1,1) = median(realized,'omitnan');
end
scatter(ax,allUB,allE,10,'filled','MarkerFaceAlpha',0.16, ...
    'MarkerEdgeAlpha',0.16);
scatter(ax,medUB,medE,45,'filled','MarkerEdgeColor','k');
maximum = max([allUB;allE;0.01]);
plot(ax,[0 maximum],[0 maximum],'k--','LineWidth',1.1);
xlim(ax,[0 1.03*maximum]);
ylim(ax,[0 1.03*maximum]);
xlabel(ax,'Analytical upper bound, $\varepsilon_{\mathrm{dec}}^{\mathrm{ub}}$', ...
    'Interpreter','latex','FontSize',fontSize);
ylabel(ax,'Realized error, $\varepsilon_{\mathrm{dec}}^{\mathrm{sim}}$', ...
    'Interpreter','latex','FontSize',fontSize);
style_axes(ax);
fractionBelow = mean(allE <= allUB+1e-10,'omitnan');
medianRatio = median(allE./max(allUB,eps),'omitnan');
text(ax,0.04,0.94, ...
    sprintf('Below bound: %.1f%%\nMedian realized/bound: %.2f', ...
    100*fractionBelow,medianRatio), ...
    'Units','normalized','VerticalAlignment','top', ...
    'FontName','Arial','FontSize',fontSize);
panel_label(ax,'c',fontSize);

ax = nexttile(layout,4);
hold(ax,'on');
ids = {'core_gam1_gb0p5','sens_qE_0p5','sens_qE_0p8','sens_qE_0p95', ...
    'sens_gate_shallow','sens_gate_intermediate','sens_gate_deep', ...
    'sens_capacity_0p01','sens_capacity_1','sens_task_specific', ...
    'sens_high_CVd','sens_high_CVw','sens_qF_0', ...
    'sens_qF_0p5','sens_small_network','sens_alternative_paths'};
scenarioLabels = {'Reference','$q_E=0.5$','$q_E=0.8$','$q_E=0.95$', ...
    'Shallow gates','Intermediate gates','Deep gates','$C/N=0.01$', ...
    '$C/N=1$','Task-specific links','$CV_{d,B}=2$','$CV_{w,B}=2$', ...
    '$q_{F,B}=0$','$q_{F,B}=0.5$','$N=250$', ...
    'Alternative paths'};
[~,reference] = get_coupled_result_by_id(catalog,results,ids{1});
referenceRho = median(reference.metrics.rho90_0,'omitnan');
referenceError = median(reference.metrics.epsilonDec,'omitnan');
yPosition = (numel(ids):-1:1)';
hRho = [];
hError = [];
xMaximum = 1;
for i = 1:numel(ids)
    [~,output] = get_coupled_result_by_id(catalog,results,ids{i});
    rho = quantile_local(output.metrics.rho90_0,[0.05 0.5 0.95])/referenceRho;
    errorValue = quantile_local(output.metrics.epsilonDec,[0.05 0.5 0.95]) ...
        /max(referenceError,eps);
    rhoHandle = errorbar(ax,rho(2),yPosition(i)+0.13, ...
        rho(2)-rho(1),rho(3)-rho(2),'horizontal');
    set(rhoHandle,'LineStyle','none','Marker','o','LineWidth',1.0, ...
        'Color',[0.20 0.45 0.75],'MarkerFaceColor',[0.20 0.45 0.75]);
    errorHandle = errorbar(ax,errorValue(2),yPosition(i)-0.13, ...
        errorValue(2)-errorValue(1),errorValue(3)-errorValue(2),'horizontal');
    set(errorHandle,'LineStyle','none','Marker','s','LineWidth',1.0, ...
        'Color',[0.85 0.40 0.20],'MarkerFaceColor',[0.85 0.40 0.20]);
    if isempty(hRho)
        hRho = rhoHandle;
        hError = errorHandle;
    end
    finiteUpper = [rho(3),errorValue(3)];
    finiteUpper = finiteUpper(isfinite(finiteUpper));
    if ~isempty(finiteUpper)
        xMaximum = max([xMaximum,finiteUpper]);
    end
end
xline(ax,1,'k--','HandleVisibility','off');
set(ax,'YTick',sort(yPosition), ...
    'YTickLabel',scenarioLabels(end:-1:1), ...
    'TickLabelInterpreter','latex');
ylim(ax,[-0.8 numel(ids)+0.6]);
xlabel(ax,'Outcome / reference median', ...
    'FontName','Arial','FontSize',fontSize);
if ~isfinite(xMaximum) || xMaximum <= 0
    xMaximum = 2.2;
end
xlim(ax,[0 1.08*max(2.2,xMaximum)]);
style_axes(ax);
legend(ax,[hRho hError], ...
    {'$\rho_{90}^{0}$','$\varepsilon_{\mathrm{dec}}^{\mathrm{sim}}$'}, ...
    'Location','southeast','NumColumns',1, ...
    'FontSize',legendFontSize,'Box','off','Interpreter','latex');
panel_label(ax,'d',fontSize);

sizeCm = fig.Position(3:4);
set(fig,'PaperUnits','centimeters','PaperSize',sizeCm, ...
    'PaperPosition',[0 0 sizeCm],'PaperPositionMode','manual', ...
    'InvertHardcopy','off');
stem=fullfile(cfg.figureDir,'Supplementary_Figure_S18');
export_publication_figure(fig,stem, ...
    cfg.figureResolution);
file=[stem '.fig']; close(fig);
end

function style_axes(ax)
set(ax,'FontName','Arial','FontSize',9,'LineWidth',0.8, ...
    'Box','off','TickDir','out','Layer','top');
grid(ax,'off');
ax.TitleFontSizeMultiplier = 1;
ax.LabelFontSizeMultiplier = 1;
end

function [t,y0,yc] = common_curve_grid(t0,R0,tc,Rc)
t = unique([t0(:);tc(:)]);
y0 = interp1(t0,R0,t,'previous','extrap');
yc = interp1(tc,Rc,t,'previous','extrap');
y0 = max(min(y0,1),0);
yc = max(min(yc,1),0);
end
