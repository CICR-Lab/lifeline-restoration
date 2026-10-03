function labels=mixed_font_xlabel(ax,plainText,latexText,fontName,fontSize,yOffset)
if nargin < 6, yOffset = 0; end
placeholder=xlabel(ax,' ','Interpreter','none', ...
    'FontName',fontName,'FontSize',fontSize);
placeholder.Units='normalized'; y=placeholder.Position(2)+yOffset;
plain=text(ax,0,y,plainText,'Units','normalized','Interpreter','none', ...
    'FontName',fontName,'FontSize',fontSize,'HorizontalAlignment','left', ...
    'VerticalAlignment','top','Clipping','off');
symbol=text(ax,0,y,latexText,'Units','normalized','Interpreter','latex', ...
    'FontSize',fontSize,'HorizontalAlignment','left', ...
    'VerticalAlignment','top','Clipping','off');
drawnow;
gap=0.018; totalWidth=plain.Extent(3)+gap+symbol.Extent(3);
plain.Position(1)=0.5-totalWidth/2;
symbol.Position(1)=plain.Position(1)+plain.Extent(3)+gap;
labels=[plain symbol];
end
