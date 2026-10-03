function figs = render_contextual_supplement(projectRoot,runDir,dataDir,figDir)
% Render editable Supplementary Figures 2-5 for the contextual analysis.
if nargin<1||isempty(projectRoot)||strlength(string(projectRoot))==0
    moduleRoot=ctx_module_root();
    projectRoot=string(fileparts(fileparts(fileparts(moduleRoot))));
end
if nargin<2||isempty(runDir)||strlength(string(runDir))==0
    if nargin<3||isempty(dataDir)||strlength(string(dataDir))==0
        error("ContextualSupplementaryInference:RunDirRequired","runDir or dataDir is required.");
    end
    runDir = "";
end
projectRoot=string(projectRoot);runDir=string(runDir);
if nargin<3||isempty(dataDir)||strlength(string(dataDir))==0,dataDir=fullfile(runDir,"figure_data");end
if nargin<4||isempty(figDir)||strlength(string(figDir))==0,figDir=fullfile(runDir,"figures");end
dataDir=string(dataDir);figDir=string(figDir);
if ~isfolder(figDir),mkdir(figDir);end
figs=gobjects(4,1);
figs(1)=renderS1(dataDir);saveOutputs(figs(1),figDir,"Supplementary_Figure_2_fitted_R2_alternative_definitions");
figs(2)=renderS3(dataDir);saveOutputs(figs(2),figDir,"Supplementary_Figure_3_adjusted_T90_factors_alternative_definitions");
figs(3)=renderS2(dataDir);saveOutputs(figs(3),figDir,"Supplementary_Figure_4_LOEO_CV_R2");
figs(4)=renderS4(dataDir);saveOutputs(figs(4),figDir,"Supplementary_Figure_5_D0_T90_model_diagnostics");
end

function fig=renderS1(dataDir)
A=readtable(fullfile(dataDir,"s1a_D0_fitted_R2_sensitivity.csv"),"TextType","string");
B=readtable(fullfile(dataDir,"s1b_T90_fitted_R2_sensitivity.csv"),"TextType","string");
C=readtable(fullfile(dataDir,"s1c_paired_fitted_delta_R2_sensitivity.csv"),"TextType","string");
A.label=cleanLabels(A.label);B.label=cleanLabels(B.label);C.label=cleanLabels(C.label);
fig=baseFigure("Supplementary Figure 2",[1 1 22.17 15.80]);
axA=axes(fig,"Position",[.27 .57 .20 .35]);axB=axes(fig,"Position",[.27 .10 .20 .35]);axC=axes(fig,"Position",[.67 .10 .29 .82]);
drawForest(axA,A,"Fitted \itR\rm^{2} for \itD\rm_{0}",false,[0.18 0.50]);
drawForest(axB,B,"Fitted \itR\rm^{2} for ln(\itT\rm_{90})",false,[0.28 0.63]);
drawForest(axC,C,"Increase in fitted \itR\rm^{2} after adding \itD\rm_{0}",true,[-.005 .14]);
panelLetter(fig,axA,"a");panelLetter(fig,axB,"b");panelLetter(fig,axC,"c");
end

function fig=renderS2(dataDir)
P=readtable(fullfile(dataDir,"s2a_primary_LOEO_performance.csv"),"TextType","string");
R=readtable(fullfile(dataDir,"s2b_paired_LOEO_delta_sensitivity.csv"),"TextType","string");
fig=baseFigure("Supplementary Figure 4",[1 1 20.90 14.24]);
axA=axes(fig,"Position",[.08 .18 .41 .67]);axB=axes(fig,"Position",[.69 .10 .27 .82]);
hold(axA,"on");axisColor=[.18 .18 .18];context=[.70 .75 .76];plus=[.29 .43 .46];x=[1 2.25 3.55 4.20];cols=[context;context;context;plus];
for i=1:height(P)
    bar(axA,x(i),P.r2_loeo_cv(i),.43,"FaceColor",cols(i,:),"EdgeColor",axisColor,"LineWidth",.55);
    errorbar(axA,x(i),P.r2_loeo_cv(i),P.r2_loeo_cv(i)-P.r2_ci95_low(i),P.r2_ci95_high(i)-P.r2_loeo_cv(i),"o","Color",axisColor,"MarkerFaceColor",cols(i,:),"MarkerSize",5,"LineWidth",.85,"CapSize",6);
end
inc=P.r2_loeo_cv(4)-P.r2_loeo_cv(3);y=.59;plot(axA,[x(3) x(3) x(4) x(4)],y+[-.015 0 0 -.015],"Color",axisColor,"LineWidth",.8);
text(axA,mean(x(3:4)),y+.025,sprintf("\\Delta\\itR\\rm^{2} = %.3f",inc),"Interpreter","tex","FontName","Arial","FontSize",7,"HorizontalAlignment","center");
axA.XLim=[.55 4.6];axA.YLim=[0 .68];axA.XTick=[x(1) x(2) mean(x(3:4))];axA.XTickLabel=["Full \itD\rm_{0} record set" "Full \itT\rm_{90} record set" "\itD\rm_{0}-aligned \itT\rm_{90} record set"];axA.XTickLabelRotation=20;
ylabel(axA,"LOEO-CV \itR\rm^{2}","Interpreter","tex");styleAxes(axA);
h1=patch(axA,nan,nan,context,"EdgeColor",axisColor);h2=patch(axA,nan,nan,plus,"EdgeColor",axisColor);
lgd=legend(axA,[h1 h2],["System-and-context model" "System-and-context-plus-\itD\rm_{0} model"],"Interpreter","tex","Location","northoutside","Orientation","horizontal","NumColumns",2,"Box","off","FontName","Arial","FontSize",7);
lgd.Units="normalized";lgd.Position=[.0650447 .856392 .462614 .0312036];
R.label=cleanLabels(arrayfun(@loeoLabel,R.specification));R.estimate=R.delta_r2_loeo_cv;R.ci95_low=R.delta_ci95_low;R.ci95_high=R.delta_ci95_high;drawForest(axB,R,"Increase in LOEO-CV \itR\rm^{2} after adding \itD\rm_{0}",true,[-.005 .14]);
panelLetter(fig,axA,"a");panelLetter(fig,axB,"b");
axA.Position=[.08 .18 .41 .67];
end

function fig=renderS3(dataDir)
S=readtable(fullfile(dataDir,"s3a_system_factor_sensitivity.csv"),"TextType","string");
D=readtable(fullfile(dataDir,"s3b_D0_factor_sensitivity.csv"),"TextType","string");
S.label=cleanLabels(S.label);D.label=cleanLabels(D.label);
fig=baseFigure("Supplementary Figure 3",[1 1 19.92 15.82]);
axA=axes(fig,"Position",[.23 .10 .25 .82]);axB=axes(fig,"Position",[.78 .10 .20 .82]);
hold(axA,"on");n=height(S);y=(n:-1:1)';off=.13;
for i=1:n
    plot(axA,[S.water_ci95_low(i) S.water_ci95_high(i)],[y(i)+off y(i)+off],"k-","LineWidth",.9);plot(axA,S.water_factor(i),y(i)+off,"ko","MarkerFaceColor","k","MarkerSize",4.5);
    plot(axA,[S.gas_ci95_low(i) S.gas_ci95_high(i)],[y(i)-off y(i)-off],"k-","LineWidth",.9);plot(axA,S.gas_factor(i),y(i)-off,"ks","MarkerFaceColor","w","MarkerSize",4.8,"LineWidth",.8);
end
xline(axA,1,"--","Color",[.6 .6 .6],"LineWidth",.7);axA.XScale="log";axA.XLim=[.8 11];axA.XTick=[1 10];axA.YTick=flipud(y);axA.YTickLabel=formatMathLabels(flipud(S.label));axA.YLim=[.4 n+.6];xlabel(axA,"Adjusted factor","Interpreter","tex");styleAxes(axA);
h1=plot(axA,nan,nan,"ko","MarkerFaceColor","k","MarkerSize",4.5);h2=plot(axA,nan,nan,"ks","MarkerFaceColor","w","MarkerSize",4.8);
lgd=legend(axA,[h1 h2],["Water supply vs electric power" "Natural gas vs electric power"],"Location","northoutside","Orientation","vertical","Box","off","FontName","Arial","FontSize",7);
lgd.Units="normalized";lgd.Position=[.257696 .877451 .194798 .0422794];
hold(axB,"on");n=height(D);y=(n:-1:1)';for i=1:n,plot(axB,[D.ci95_low(i) D.ci95_high(i)],[y(i) y(i)],"k-","LineWidth",.9);plot(axB,D.estimate(i),y(i),"ko","MarkerFaceColor","w","MarkerSize",4.7,"LineWidth",.8);end
xline(axB,1,"--","Color",[.6 .6 .6],"LineWidth",.7);axB.XScale="log";axB.XLim=[.98 1.30];axB.YTick=flipud(y);axB.YTickLabel=formatMathLabels(flipud(D.label));axB.YLim=[.4 n+.6];xlabel(axB,"Adjusted factor for \itD\rm_{0}: 0.10 to 0.20","Interpreter","tex");styleAxes(axB);
panelLetter(fig,axA,"a");panelLetter(fig,axB,"b");
axA.Position=[.23 .10 .25 .726472816611781];
end

function fig=renderS4(dataDir)
C=readtable(fullfile(dataDir,"s4b_D0_calibration.csv"));
F=readtable(fullfile(dataDir,"s4_full_T90_residuals.csv"));FQ=readtable(fullfile(dataDir,"s4_full_T90_qq.csv"));
A=readtable(fullfile(dataDir,"s4_aligned_plus_D0_residuals.csv"));AQ=readtable(fullfile(dataDir,"s4_aligned_plus_D0_qq.csv"));
fig=baseFigure("Supplementary Figure 5",[1 1 17.09 11.91]);
pos=[.07 .33 .25 .34;.39 .57 .25 .34;.72 .57 .25 .34;.39 .10 .25 .34;.72 .10 .25 .34];ax=gobjects(5,1);for i=1:5,ax(i)=axes(fig,"Position",pos(i,:));hold(ax(i),"on");end
plot(ax(1),[0 1],[0 1],"--","Color",[.6 .6 .6]);scatter(ax(1),C.mean_fitted_D0,C.mean_observed_D0,25,"k","filled");xlabel(ax(1),"Mean fitted \itD\rm_{0}","Interpreter","tex");ylabel(ax(1),"Mean observed \itD\rm_{0}","Interpreter","tex");axis(ax(1),"square");
scatter(ax(2),F.fitted_log_time,F.residual,8,[.45 .45 .45],"filled","MarkerFaceAlpha",.45);yline(ax(2),0,"--","Color",[.6 .6 .6]);xlabel(ax(2),"Fitted ln(\itT\rm_{90})","Interpreter","tex");ylabel(ax(2),"Residual");
plot(ax(3),FQ.normal_quantile,FQ.residual_quantile,".","Color",[.35 .35 .35],"MarkerSize",6);qqReference(ax(3),FQ.normal_quantile,FQ.residual_quantile);xlabel(ax(3),"Normal quantile");ylabel(ax(3),"Residual quantile");
scatter(ax(4),A.fitted_log_time,A.residual,8,[.45 .45 .45],"filled","MarkerFaceAlpha",.45);yline(ax(4),0,"--","Color",[.6 .6 .6]);xlabel(ax(4),"Fitted ln(\itT\rm_{90})","Interpreter","tex");ylabel(ax(4),"Residual");
plot(ax(5),AQ.normal_quantile,AQ.residual_quantile,".","Color",[.35 .35 .35],"MarkerSize",6);qqReference(ax(5),AQ.normal_quantile,AQ.residual_quantile);xlabel(ax(5),"Normal quantile");ylabel(ax(5),"Residual quantile");
for i=1:5,styleAxes(ax(i));panelLetter(fig,ax(i),char('a'+i-1));end
ax(1).XTick=[0 .5 1];ax(1).YTick=0:.2:1;
ax(2).XLim=[-2 4.25660500177061];ax(2).XTick=-2:2:4;
ax(4).XLim=[-2.29540063019031 4];ax(4).YLim=[-6 4];ax(4).XTick=-2:2:4;ax(4).YTick=-6:2:4;
ax(5).YLim=[-6 4];ax(5).YTick=-6:2:4;
end

function drawForest(ax,T,xLabel,zeroLine,xLimits)
hold(ax,"on");n=height(T);y=(n:-1:1)';markerColor=[.40 .40 .40];
for i=1:n,plot(ax,[T.ci95_low(i) T.ci95_high(i)],[y(i) y(i)],"Color",markerColor,"LineWidth",1);plot(ax,T.estimate(i),y(i),"o","Color",markerColor,"MarkerFaceColor","w","MarkerSize",4.8,"LineWidth",.8);end
if zeroLine,xline(ax,0,"--","Color",[.6 .6 .6],"LineWidth",.7);end
ax.YTick=flipud(y);ax.YTickLabel=formatMathLabels(flipud(T.label));ax.YLim=[.4 n+.6];ax.XLim=xLimits;xlabel(ax,xLabel,"Interpreter","tex");styleAxes(ax);
end
function label=loeoLabel(id)
switch string(id)
    case "log2_d0",label="log2(D0)";case "raw_d0",label="Raw D0";case "log2_d0_complete",label="log2(D0) + complete-disruption indicator";case "spline_log2_d0",label="Spline in log2(D0)";case "log2p1_d0_zero_inclusive",label="log2(1 + D0)";case "T80",label="T80";case "T95",label="T95";case "popWeightedPGA_g",label="Population-weighted PGA";case "popMeanPGA_g",label="Mean PGA across populated cells";case "pgaMedian_g",label="Median PGA across populated cells";case "pgaMax_g",label="Maximum PGA across populated cells";case "ln1p_pga_over_0p1",label="ln(1 + PGA / 0.1)";case "population_density",label="Population density";case "exclude_bottom_1pct",label="Exclude bottom 1% of positive D0";case "exclude_bottom_5pct",label="Exclude bottom 5% of positive D0";case "lognormal",label="Lognormal";case "weibull",label="Weibull";case "loglogistic",label="Log-logistic";otherwise,label=replace(string(id),"_"," ");end
end
function labels=cleanLabels(labels)
labels=string(labels);
labels=replace(labels,["ln(1 + PGA / 0.1)","ln(1 + PGA / 0.1 g)","ln(1 + PGA / (0.1 g))"],"ln[1 + PGA / (0.1 g)]");
labels(labels=="Primary fractional logit")="Fractional logit";
labels(labels=="Primary model")="System-and-context model";
labels(labels=="log2(D0) (primary)")="log2(D0)";
labels(labels=="Population-weighted mean PGA")="Mean PGA across populated cells";
labels(labels=="Median PGA")="Median PGA across populated cells";
labels(labels=="Maximum PGA")="Maximum PGA across populated cells";
end
function labels=formatMathLabels(labels)
labels=string(labels);
labels=replace(labels,"log2(D0)","log_{2}(\itD\rm_{0})");
labels=replace(labels,"log2(1 + D0)","log_{2}(1 + \itD\rm_{0})");
labels=replace(labels,"Raw D0","Raw \itD\rm_{0}");
labels=replace(labels,"positive D0","positive \itD\rm_{0}");
labels=replace(labels,"T80","\itT\rm_{80}");
labels=replace(labels,"T90","\itT\rm_{90}");
labels=replace(labels,"T95","\itT\rm_{95}");
end
function fig=baseFigure(name,position),fig=figure("Name",name,"NumberTitle","off","Color","w","Units","centimeters","Position",position,"Visible","on");end
function styleAxes(ax),ax.FontName="Arial";ax.TickLabelInterpreter="tex";ax.FontSize=7.2;ax.LineWidth=.75;ax.XColor=[.18 .18 .18];ax.YColor=[.18 .18 .18];ax.TickDir="out";ax.Box="off";ax.XGrid="off";ax.YGrid="off";ax.Layer="top";end
function panelLetter(fig,ax,label),p=ax.Position;annotation(fig,"textbox",[p(1)-.035 p(2)+p(4)+.012 .025 .025],"String",label,"LineStyle","none","FontName","Arial","FontSize",10,"FontWeight","bold","Margin",0,"Color",[.18 .18 .18]);end
function qqReference(ax,x,y),q=abs(x)<.7;b=polyfit(x(q),y(q),1);xx=[min(x) max(x)];plot(ax,xx,polyval(b,xx),"--","Color",[.6 .6 .6],"LineWidth",.7);end
function saveOutputs(fig,figDir,name)
drawnow;figPath=fullfile(figDir,name+".fig");pngPath=fullfile(figDir,name+".png");savefig(fig,figPath);exportgraphics(fig,pngPath,"Resolution",300);padPng(pngPath,30,300);check=openfig(figPath,"visible");drawnow;reopenedPath=fullfile(figDir,name+"_reopened.png");exportgraphics(check,reopenedPath,"Resolution",180);padPng(reopenedPath,18,180);close(check);
end
function padPng(path,padding,dpi)
imageData=imread(path);canvas=255*ones(size(imageData,1)+2*padding,size(imageData,2)+2*padding,size(imageData,3),"like",imageData);canvas(padding+(1:size(imageData,1)),padding+(1:size(imageData,2)),:)=imageData;imwrite(canvas,path,"png","ResolutionUnit","meter","XResolution",round(dpi/0.0254),"YResolution",round(dpi/0.0254));
end
