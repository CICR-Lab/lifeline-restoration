function file=plot_supplementary_figure_s14(root,mode,dataRoot,outputRoot)
if nargin<3 || isempty(dataRoot), dataRoot=resolve_figure_data_root(root,mode); end
if nargin<4 || isempty(outputRoot), outputRoot=resolve_figure_output_root(root,mode); end
sourceFile=fullfile(dataRoot,'Analytical','Analytical_results.mat');
if ~isfile(sourceFile), error('Missing analytical results. Run run_hierarchical_repair_analysis first.'); end
loaded=load(sourceFile,'s14','cfg'); data=loaded.s14; cfg=loaded.cfg;
cfg.output.supp=fullfile(outputRoot,'supplementary');
if ~isfolder(cfg.output.supp), mkdir(cfg.output.supp); end
t=data.caseTable; summary=data.fitSummary;
fig=figure('Color','w','Units','centimeters','Position',[2 2 12.5 9]);
ax=axes(fig); hold(ax,'on'); colors=[cfg.colors.blue;cfg.colors.orange];
markers={'o','s'};
for r=1:numel(data.repairLimits)
    keep=t.RepairLimit==string(data.repairLimits{r});
    scatter(ax,t.DmaxTau50(keep),t.LossRest(keep),22,colors(r,:), ...
        markers{r},'filled','MarkerFaceAlpha',.48, ...
        'DisplayName',data.limitLabels{r});
end
xMax=max(t.DmaxTau50); yMax=max(t.LossRest);
xLine=linspace(0,1.03*xMax,200)';
plot(ax,xLine,summary.K(1)*xLine,'k-','LineWidth',1.5, ...
    'HandleVisibility','off');
xlabel(ax,'$D_{\max}\tau_{50}$','Interpreter','latex');
ylabel(ax,'$L_{\mathrm{rest}}$','Interpreter','latex');
xlim(ax,[0 1.04*xMax]); ylim(ax,[0 1.04*yMax]);
legend(ax,'Location','northwest','Orientation','vertical','Box','off', ...
    'FontSize',cfg.figure.fontSize);
apply_axes_style(ax,cfg);
stem=fullfile(cfg.output.supp,'Supplementary_Figure_S14');
export_publication_figure(fig,stem,cfg.figure.dpi); file=[stem '.fig']; close(fig);
end
