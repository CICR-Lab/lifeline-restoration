function figures=run_figure_4(varargin)
% Render Figure 4b-h from packaged or generated figure-data inputs.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'Code','Analytical'));
addpath(fullfile(root,'Code','Simulation'));
addpath(fullfile(root,'Code','Interdependent'));
addpath(fullfile(root,'Code','Common'));
[mode,dataSource,dataRoot,outputRoot] = parseFigureOptions(root,varargin{:});

dataRoot=resolve_figure_data_root(root,mode,dataSource,dataRoot);
figures.analytical=plot_analytical_figure4(root,mode,dataRoot,outputRoot);
figures.figure4e=figure4e_global_sensitivity_matlab(root,mode,dataRoot,outputRoot);
figures.figure4f=figure4f_task_service_matlab(root,mode,dataRoot,outputRoot);
figures.figure4g=figure4g_capacity_priority_matlab(root,mode,dataRoot,outputRoot);
figures.figure4h=plot_interdependent_figure4h(root,mode,dataRoot,outputRoot);
fprintf('\nFigure 4 panels completed in:\n  %s\n', ...
    fullfile(outputRoot,'main'));
end

function [mode,dataSource,dataRoot,outputRoot] = parseFigureOptions(root,varargin)
mode='full';
if ~isempty(varargin) && isTextScalar(varargin{1}) && ...
        ismember(lower(char(varargin{1})),{'quick','full'})
    mode=lower(char(varargin{1}));
    varargin=varargin(2:end);
end
p=inputParser;
p.addParameter('Mode',mode,@(x)isTextScalar(x));
p.addParameter('DataSource','packaged',@(x)isTextScalar(x));
p.addParameter('DataRoot','',@(x)isTextScalar(x));
p.addParameter('OutputDir','',@(x)isTextScalar(x));
p.parse(varargin{:});
mode=lower(char(p.Results.Mode));
dataSource=lower(char(p.Results.DataSource));
if ~ismember(mode,{'quick','full'}), error('Mode must be ''quick'' or ''full''.'); end
if ~ismember(dataSource,{'packaged','generated'})
    error('DataSource must be ''packaged'' or ''generated''.');
end
dataRoot=char(p.Results.DataRoot);
outputRoot=char(p.Results.OutputDir);
if isempty(outputRoot), outputRoot=resolve_figure_output_root(root,mode); end
end

function tf=isTextScalar(value)
tf=(ischar(value) && (isrow(value) || isempty(value))) || ...
    (isstring(value) && isscalar(value));
end

function outPng=assembleFigure4Png(mainDir,~)
outPng=fullfile(mainDir,'Figure4.png');
fig=figure('Color','w','Units','pixels','Position',[80 80 1500 1120], ...
    'Visible','off','InvertHardcopy','off');
panelFiles=struct( ...
    'b',"Figure4b_depth_curves.png", ...
    'c',"Figure4c_phase_ratios.png", ...
    'd',"Figure4d_service_distribution.png", ...
    'e',"Figure4e_global_sensitivity_matlab.png", ...
    'f',"Figure4f_task_service_matlab.png", ...
    'g',"Figure4g_capacity_priority_matlab.png", ...
    'h',"Figure4h_coupled_network_contours.png");
rect=struct( ...
    'a',[0.035 0.675 0.285 0.285], ...
    'b',[0.335 0.680 0.275 0.280], ...
    'c',[0.615 0.680 0.315 0.280], ...
    'd',[0.035 0.365 0.285 0.280], ...
    'e',[0.300 0.365 0.675 0.280], ...
    'f',[0.035 0.060 0.285 0.255], ...
    'g',[0.335 0.060 0.275 0.255], ...
    'h',[0.615 0.060 0.315 0.255]);
names=fieldnames(panelFiles);
for i=1:numel(names)
    drawTrimmedPngPanel(fig,fullfile(mainDir,panelFiles.(names{i})),rect.(names{i}));
end
addPanelLabels(fig,rect);
set(fig,'PaperPositionMode','auto');
print(fig,outPng,'-dpng','-r300');
close(fig);
end

function drawTrimmedPngPanel(fig,pngPath,targetRect)
img=imread(pngPath);
img=trimWhiteBorder(img,12);
ax=axes(fig,'Position',targetRect);
image(ax,img);
axis(ax,'image');
axis(ax,'off');
end

function out=trimWhiteBorder(img,pad)
mask=any(img<248,3);
[r,c]=find(mask);
if isempty(r)
    out=img;
    return
end
r1=max(1,min(r)-pad);
r2=min(size(img,1),max(r)+pad);
c1=max(1,min(c)-pad);
c2=min(size(img,2),max(c)+pad);
out=img(r1:r2,c1:c2,:);
end

function drawFigure4ePanel(fig,sourceFile,targetRect)
t=readtable(sourceFile);
fontName='Arial';
cellFontSize=9.0;
rowHeaderFontSize=9.0;
columnHeaderFontSize=10;
titleFontSize=9.0;
ax=axes(fig,'Position',[ ...
    targetRect(1)+0.24*targetRect(3), ...
    targetRect(2)+0.04*targetRect(4), ...
    0.74*targetRect(3), ...
    0.88*targetRect(4)]);

outcomeNames={ ...
    'log_tau80_tau50', ...
    'log_tau90_tau80', ...
    'log_tau95_tau90', ...
    'log_kappa', ...
    'logit_eta90'};
columnLabels={ ...
    '$\tau_{80}/\tau_{50}$', ...
    '$\tau_{90}/\tau_{80}$', ...
    '$\tau_{95}/\tau_{90}$', ...
    '$\kappa$', ...
    '$\eta_{90}$'};
rowText={ ...
    'Repair-task count, '; ...
    'Maximum dependency depth, '; ...
    'Task-layer profile'; ...
    'Task-count concentration'; ...
    'Reconnection-prerequisite configuration'; ...
    'Parent-choice concentration, '; ...
    'Layer service-share profile, '; ...
    'Deep-to-shallow duration ratio'; ...
    'Within-layer repair-duration, '; ...
    'Within-layer service-weight, '; ...
    'Duration-service association, '; ...
    'Relative repair capacity, '; ...
    'Frontier-priority probability, '};
rowSymbols={ ...
    '$N$'; ...
    '$L$'; ...
    ''; ...
    ''; ...
    ''; ...
    '$\alpha_p$'; ...
    '$\pi_l$'; ...
    ''; ...
    '$CV_d$'; ...
    '$CV_w$'; ...
    '$\rho_{dw}$'; ...
    '$C/N$'; ...
    '$q_F$'};

blue=[0.173 0.498 0.722];
paleBlue=0.82+0.18*blue;
outcomeColumn=cellstr(t.outcomeName);
totalOrder=nan(13,5);
firstOrder=nan(13,5);
for row=1:13
    for column=1:5
        match=t.factorID==row & strcmp(outcomeColumn,outcomeNames{column});
        totalOrder(row,column)=100*t.totalOrderRaw(match);
        firstOrder(row,column)=100*t.firstOrderRaw(match);
    end
end

hold(ax,'on');
for row=1:13
    faceColor=[1 1 1];
    if ismember(row,[3 7 12])
        faceColor=paleBlue;
    end
    for column=1:5
        rectangle(ax,'Position',[column-0.5 row-0.5 1 1], ...
            'FaceColor',faceColor,'EdgeColor',[0.84 0.86 0.88], ...
            'LineWidth',0.55);
        totalValue=totalOrder(row,column);
        firstValue=firstOrder(row,column);
        if round(totalValue,1)==0, totalValue=0; end
        if round(firstValue,1)==0, firstValue=0; end
        text(ax,column,row,sprintf('{\\bf %.1f} | %.1f',totalValue,firstValue), ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'Interpreter','tex','FontName',fontName,'FontSize',cellFontSize);
    end
end
plot(ax,[0.5 5.5],[6.5 6.5],'-','Color',[0.35 0.38 0.40],'LineWidth',1);
plot(ax,[0.5 5.5],[11.5 11.5],'-','Color',[0.35 0.38 0.40],'LineWidth',1);
hold(ax,'off');
set(ax,'XLim',[0.5 5.5],'YLim',[0.5 13.5],'YDir','reverse', ...
    'XTick',1:5,'XTickLabel',columnLabels,'YTick',1:13, ...
    'YTickLabel',[],'TickLength',[0 0],'FontName',fontName, ...
    'FontSize',rowHeaderFontSize,'TickLabelInterpreter','latex','Box','off');
ax.XAxis.FontSize=columnHeaderFontSize;
for row=1:13
    drawEPanelRowLabel(ax,row,rowText{row},rowSymbols{row},fontName,rowHeaderFontSize);
end
title(ax,'Global sensitivity: total-order | first-order (% of input-driven variance)', ...
    'FontName',fontName,'FontSize',titleFontSize,'FontWeight','normal', ...
    'Interpreter','none');
end

function drawEPanelRowLabel(ax,row,plainText,latexText,fontName,fontSize)
y=1-(row-0.5)/13;
rightEdge=-0.018;
if isempty(latexText)
    text(ax,rightEdge,y,plainText,'Units','normalized', ...
        'HorizontalAlignment','right','VerticalAlignment','middle', ...
        'Interpreter','none','FontName',fontName,'FontSize',fontSize, ...
        'Clipping','off');
    return
end
symbol=text(ax,rightEdge,y,latexText,'Units','normalized', ...
    'HorizontalAlignment','right','VerticalAlignment','middle', ...
    'Interpreter','latex','FontSize',fontSize,'Clipping','off');
drawnow;
text(ax,rightEdge-symbol.Extent(3)-0.005,y,plainText,'Units','normalized', ...
    'HorizontalAlignment','right','VerticalAlignment','middle', ...
    'Interpreter','none','FontName',fontName,'FontSize',fontSize, ...
    'Clipping','off');
end

function copyFigPanel(figPath,targetFig,targetRect)
src=openfig(figPath,'invisible');
cleanup=onCleanup(@()close(src));
objs=panelObjects(src);
box=objectUnion(objs);
newObjs=copyobj(objs,targetFig);
for i=1:numel(objs)
    placeCopiedObject(newObjs(i),objs(i),box,targetRect);
end
removeCopiedPanelLabels(newObjs);
scalePanelGraphics(newObjs,1.00);
end

function objs=panelObjects(src)
ax=findall(src,'Type','axes');
leg=findall(src,'Type','legend');
cb=findall(src,'Type','colorbar');
ax=sortAxesByArea(ax);
objs=[ax(:); cb(:); leg(:)];
end

function ax=sortAxesByArea(ax)
areas=zeros(numel(ax),1);
for i=1:numel(ax)
    set(ax(i),'Units','normalized');
    p=get(ax(i),'Position');
    areas(i)=p(3)*p(4);
end
[~,idx]=sort(areas,'descend');
ax=ax(idx);
end

function box=objectUnion(objs)
x1=inf;
y1=inf;
x2=-inf;
y2=-inf;
for i=1:numel(objs)
    set(objs(i),'Units','normalized');
    p=get(objs(i),'Position');
    x1=min(x1,p(1));
    y1=min(y1,p(2));
    x2=max(x2,p(1)+p(3));
    y2=max(y2,p(2)+p(4));
end
box=[x1 y1 x2-x1 y2-y1];
end

function placeCopiedObject(obj,srcObj,box,targetRect)
set(srcObj,'Units','normalized');
p=get(srcObj,'Position');
newPos=[ ...
    targetRect(1)+(p(1)-box(1))/box(3)*targetRect(3), ...
    targetRect(2)+(p(2)-box(2))/box(4)*targetRect(4), ...
    p(3)/box(3)*targetRect(3), ...
    p(4)/box(4)*targetRect(4)];
set(obj,'Units','normalized','Position',newPos);
end

function removeCopiedPanelLabels(objs)
letters=["a","b","c","d","e","f","g","h"];
for i=1:numel(objs)
    txt=findall(objs(i),'Type','text');
    for j=1:numel(txt)
        s=string(get(txt(j),'String'));
        if isscalar(s) && any(s==letters)
            delete(txt(j));
        end
    end
end
end

function scalePanelGraphics(objs,scale)
items=findall(objs);
for i=1:numel(items)
    if isprop(items(i),'FontSize')
        items(i).FontSize=items(i).FontSize*scale;
    end
    if isprop(items(i),'LineWidth'), items(i).LineWidth=items(i).LineWidth*1.05; end
    if isprop(items(i),'MarkerSize'), items(i).MarkerSize=items(i).MarkerSize*1.05; end
end
end

function addPanelLabels(fig,rect)
labels={'a','b','c','d','e','f','g','h'};
for i=1:numel(labels)
    r=rect.(labels{i});
    annotation(fig,'textbox',[r(1)-0.015 r(2)+r(4)-0.012 0.03 0.03], ...
        'String',labels{i},'EdgeColor','none','FitBoxToText','off', ...
        'Margin',0,'FontName','Arial','FontWeight','bold', ...
        'FontSize',20,'HorizontalAlignment','left', ...
        'VerticalAlignment','top');
end
end
