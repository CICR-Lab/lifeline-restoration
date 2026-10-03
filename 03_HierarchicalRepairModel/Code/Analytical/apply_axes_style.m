function apply_axes_style(ax,cfg)
set(ax,'FontName',cfg.figure.fontName,'FontSize',cfg.figure.fontSize, ...
    'LineWidth',0.8,'Box','off','TickDir','out','Layer','top');
grid(ax,'off');
end
