function export_publication_figure(fig,basePath,resolution)
[folder,~,~]=fileparts(basePath);
if ~isfolder(folder), mkdir(folder); end
set(fig,'Color','w'); drawnow;
if exist('exportgraphics','file')==2
    exportgraphics(fig,[basePath '.png'],'Resolution',resolution);
else
    oldMode=get(fig,'PaperPositionMode');
    cleanup=onCleanup(@() set(fig,'PaperPositionMode',oldMode)); %#ok<NASGU>
    set(fig,'PaperPositionMode','auto');
    print(fig,[basePath '.png'],'-dpng',sprintf('-r%d',resolution));
end
save_fig_copy(fig,[basePath '.fig']);
end

function save_fig_copy(fig,figPath)
figCopy=copyobj(fig,groot);
cleanup=onCleanup(@() close(figCopy));
prepare_fig_file_window(figCopy,fig);
savefig(figCopy,figPath,'-v7');
end

function prepare_fig_file_window(fig,sourceFig)
oldUnits=sourceFig.Units;
sourceFig.Units='pixels';
pos=sourceFig.Position;
sourceFig.Units=oldUnits;
screen = get(groot,'ScreenSize');
maxW = 0.65 * screen(3);
maxH = 0.72 * screen(4);
scale = min([1, maxW / pos(3), maxH / pos(4)]);
pos(3:4) = max([260 210], pos(3:4) * scale);
pos(1) = max(40, screen(1) + 0.04 * screen(3));
pos(2) = max(40, screen(2) + 0.08 * screen(4));
fig.Units = 'pixels';
fig.Position = pos;
fig.OuterPosition = pos;
fig.Visible = 'on';
fig.Resize = 'on';
fig.DockControls = 'off';
fig.WindowStyle = 'normal';
fig.WindowState = 'normal';
fig.CreateFcn = make_resize_callback(pos);
drawnow;
end

function callback = make_resize_callback(pos)
callback = sprintf([ ...
    'set(gcbf,''WindowStyle'',''normal'',''WindowState'',''normal'',' ...
    '''Units'',''pixels'',''Resize'',''on'',''DockControls'',''off'',' ...
    '''OuterPosition'',[%g %g %g %g],''Position'',[%g %g %g %g]);'], ...
    pos, pos);
end
