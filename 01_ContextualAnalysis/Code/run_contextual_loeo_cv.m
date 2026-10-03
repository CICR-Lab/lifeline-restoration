function [primary, baseline] = run_contextual_loeo_cv(T, cfg)
d0 = T(T.d0_model_eligible,:); fullT90 = T(T.t90_context_eligible,:); aligned = T(T.t90_primary_eligible,:);
[d0Prediction,~] = d0Oof(d0,cfg);
[fullPrediction,~] = t90Oof(fullT90,cfg.primary_pga,"none","lognormal",nan,cfg);
[alignedContext,alignedPlus] = t90Oof(aligned,cfg.primary_pga,"log2_d0","lognormal",nan,cfg);
primary = [predictionRow("D0","system_and_context",d0.outage0,d0Prediction,d0.eid); ...
    predictionRow("T90","system_and_context_full_record_set",log(fullT90.T90),fullPrediction,fullT90.eid); ...
    predictionRow("T90","system_and_context_d0_aligned_record_set",log(aligned.T90),alignedContext,aligned.eid); ...
    predictionRow("T90","system_and_context_plus_d0",log(aligned.T90),alignedPlus,aligned.eid)];
boot = primaryBootstrap(d0,d0Prediction,fullT90,fullPrediction,aligned,alignedContext,alignedPlus,cfg);
primary = addIntervals(primary,boot);
primary.delta_r2_loeo_cv = [nan; nan; nan; primary.r2_loeo_cv(4)-primary.r2_loeo_cv(3)];
deltaCI = prctile(boot.delta_r2,[2.5 97.5]); primary.delta_ci95_low = [nan; nan; nan; deltaCI(1)];
primary.delta_ci95_high = [nan; nan; nan; deltaCI(2)];
baseline = struct("record_id", aligned.record_id, "eid", aligned.eid, ...
    "observed", log(aligned.T90), "context_prediction", alignedContext, ...
    "plus_prediction", alignedPlus, "summary", primary(4,:), "config", cfg);
baseline.d0 = struct("record_id", d0.record_id, "eid", d0.eid, ...
    "observed", d0.outage0, "prediction", d0Prediction, "summary", primary(1,:));
end
function directPredictionPrimary(cfg,runDir)
if isfile(fullfile(runDir,"models","primary_prediction_performance.csv")), fprintf("Primary prediction stage already complete; reusing saved outputs.\n"); return; end
S = load(fullfile(runDir,"input","canonical_analysis_table.mat"),"T"); T = S.T;
rng(cfg.random_seed,"twister");
d0 = T(T.d0_model_eligible,:); fullT90 = T(T.t90_context_eligible,:); aligned = T(T.t90_primary_eligible,:);
[d0Prediction,~] = d0Oof(d0,cfg); [fullPrediction,~] = t90Oof(fullT90,cfg.primary_pga,"none","lognormal",nan,cfg);
[alignedContext,alignedPlus] = t90Oof(aligned,cfg.primary_pga,"log2_d0","lognormal",nan,cfg);
metrics = [predictionRow("D0","system_and_context",d0.outage0,d0Prediction,d0.eid); ...
    predictionRow("T90","system_and_context_full_record_set",log(fullT90.T90),fullPrediction,fullT90.eid); ...
    predictionRow("T90","system_and_context_d0_aligned_record_set",log(aligned.T90),alignedContext,aligned.eid); ...
    predictionRow("T90","system_and_context_plus_d0",log(aligned.T90),alignedPlus,aligned.eid)];
boot = primaryBootstrap(d0,d0Prediction,fullT90,fullPrediction,aligned,alignedContext,alignedPlus,cfg);
metrics = addIntervals(metrics,boot);
writetable(metrics,fullfile(runDir,"models","primary_prediction_performance.csv"));
paired = pairedSummary(log(aligned.T90),alignedContext,alignedPlus,aligned.eid,boot);
writetable(paired,fullfile(runDir,"models","primary_prediction_increment.csv"));
writePredictionOutputs(d0,d0Prediction,fullT90,fullPrediction,aligned,alignedContext,alignedPlus,runDir);
writeStageStatus(runDir,"prediction_primary","complete",struct("d0_record_count",height(d0),"full_t90_record_count",height(fullT90),"aligned_t90_record_count",height(aligned)));
fprintf("Primary prediction complete: D0 n=%d; full T90 n=%d; aligned T90 n=%d.\n",height(d0),height(fullT90),height(aligned));
end


function directPredictionRobustness(cfg,runDir)
if isfile(fullfile(runDir,"validation","prediction_robustness_unified.csv")), fprintf("Unified prediction robustness stage already complete; reusing saved outputs.\n"); return; end
S=load(fullfile(runDir,"input","canonical_analysis_table.mat"),"T"); robust=ctx_unified_prediction_robustness(S.T,cfg);
writetable(robust,fullfile(runDir,"validation","prediction_robustness_unified.csv"));
writetable(robust,fullfile(runDir,"figure_data","supplementary_figure_sx_data.csv"));
writeStageStatus(runDir,"prediction_robustness","complete",struct("scenario_count",height(robust),"failed_scenarios",sum(robust.status~="complete")));
end

function [prediction,fold] = d0Oof(D,cfg)
events = unique(D.eid,"stable"); prediction = nan(height(D),1); fold = table();
levels = incomeLevels(D,cfg);
for k=1:numel(events)
    test = D.eid==events(k); train = ~test;
    [Xtr,Xte] = fixedDesign(D(train,:),D(test,:),levels,cfg.primary_pga,cfg.primary_pga_transform,"none",false,[]);
    m = fitglm(Xtr,D.outage0(train),"linear","Distribution","binomial","Link","logit");
    prediction(test) = predict(m,Xte);
end
end

function [contextPrediction,plusPrediction] = t90Oof(D,pgaField,d0Kind,distribution,smallQuantile,cfg)
events = unique(D.eid,"stable"); contextPrediction = nan(height(D),1); plusPrediction = nan(height(D),1); levels = incomeLevels(D,cfg); y = log(D.T90);
for k=1:numel(events)
    test = D.eid==events(k); train = ~test; cutoff = -inf;
    if isfinite(smallQuantile), cutoff = quantile(D.outage0(train),smallQuantile); test = test & D.outage0>cutoff; train = train & D.outage0>cutoff; end
    if ~any(test), continue; end
    [XcTr,XcTe] = fixedDesign(D(train,:),D(test,:),levels,pgaField,cfg.primary_pga_transform,"none",false,[]);
    knots = [];
    if d0Kind == "spline_log2_d0", knots = quantile(log2(D.outage0(train)),[.05 .35 .65 .95]); end
    if d0Kind == "none"
        XpTr = XcTr; XpTe = XcTe;
    else
        [XpTr,XpTe] = fixedDesign(D(train,:),D(test,:),levels,pgaField,cfg.primary_pga_transform,d0Kind,true,knots);
    end
    if distribution == "lognormal"
        bc = linearFit(XcTr,y(train)); bp = linearFit(XpTr,y(train));
        contextPrediction(test) = [ones(sum(test),1),XcTe]*bc; plusPrediction(test) = [ones(sum(test),1),XpTe]*bp;
    else
        cMean=mean(XcTr,1); cSd=std(XcTr,0,1); cSd(cSd==0 | ~isfinite(cSd))=1; XcTrScaled=(XcTr-cMean)./cSd; XcTeScaled=(XcTe-cMean)./cSd;
        pMean=mean(XpTr,1); pSd=std(XpTr,0,1); pSd(pSd==0 | ~isfinite(pSd))=1; XpTrScaled=(XpTr-pMean)./pSd; XpTeScaled=(XpTe-pMean)./pSd;
        tc = fitAft(XcTrScaled,y(train),distribution); tp = fitAft(XpTrScaled,y(train),distribution);
        contextPrediction(test) = aftMedian(XcTeScaled,tc,distribution); plusPrediction(test) = aftMedian(XpTeScaled,tp,distribution);
    end
end
end

function [Xtr,Xte] = fixedDesign(train,test,levels,pgaField,pgaTransform,d0Kind,includeD0,knots)
C = [ctx_pga_transform(train.(pgaField),pgaTransform),(train.year-2000)/10,log2(train.population)];
Ct = [ctx_pga_transform(test.(pgaField),pgaTransform),(test.year-2000)/10,log2(test.population)];
if includeD0
    switch d0Kind
        case "log2_d0", C(:,end+1)=log2(train.outage0); Ct(:,end+1)=log2(test.outage0);
        case "raw_d0", C(:,end+1)=train.outage0; Ct(:,end+1)=test.outage0;
        case "log2_d0_complete", C(:,end+1)=log2(train.outage0); C(:,end+1)=double(train.outage0==1); Ct(:,end+1)=log2(test.outage0); Ct(:,end+1)=double(test.outage0==1);
        case "log2p1_d0", C(:,end+1)=log2(1+train.outage0); Ct(:,end+1)=log2(1+test.outage0);
        case "spline_log2_d0", C=[C rcsBasis(log2(train.outage0),knots)]; Ct=[Ct rcsBasis(log2(test.outage0),knots)];
        otherwise, error("ContextualDirect:UnknownD0","Unknown D0 encoding: %s",d0Kind);
    end
end
sys = ["E";"W";"G"]; S = zeros(height(train),2); St=zeros(height(test),2);
for j=2:3, S(:,j-1)=train.sys==sys(j); St(:,j-1)=test.sys==sys(j); end
I=zeros(height(train),numel(levels)-1); It=zeros(height(test),numel(levels)-1);
for j=2:numel(levels), I(:,j-1)=train.development_level==levels(j); It(:,j-1)=test.development_level==levels(j); end
Xtr=[C S I]; Xte=[Ct St It];
end

function T = predictionRow(outcome,model,y,p,eid)
ok=isfinite(p); m=metric(y(ok),p(ok)); T=table(string(outcome),string(model),sum(ok),numel(unique(eid(ok))),m.r2,m.rmse,m.mae, ...
    'VariableNames',{'outcome','result','record_count','earthquake_count','r2_loeo_cv','rmse','mae'});
end

function boot = primaryBootstrap(d0,d0p,full,fullp,aligned,ctx,plus,cfg)
rng(cfg.random_seed+100,"twister"); n=cfg.bootstrap_replicates; values=nan(n,12);
for b=1:n
    a=sampleEventRows(d0.eid); f=sampleEventRows(full.eid); p=sampleEventRows(aligned.eid);
    dm=metric(d0.outage0(a),d0p(a)); fm=metric(log(full.T90(f)),fullp(f)); cm=metric(log(aligned.T90(p)),ctx(p)); pm=metric(log(aligned.T90(p)),plus(p));
    values(b,:)=[dm.r2 dm.rmse dm.mae fm.r2 fm.rmse fm.mae cm.r2 cm.rmse cm.mae pm.r2 pm.rmse pm.mae];
end
boot=array2table(values,'VariableNames',{'d0_r2','d0_rmse','d0_mae','full_r2','full_rmse','full_mae','aligned_context_r2','aligned_context_rmse','aligned_context_mae','aligned_plus_r2','aligned_plus_rmse','aligned_plus_mae'});
boot.delta_r2=boot.aligned_plus_r2-boot.aligned_context_r2; boot.delta_rmse=boot.aligned_plus_rmse-boot.aligned_context_rmse; boot.delta_mae=boot.aligned_plus_mae-boot.aligned_context_mae; boot.replicate=(1:n)';
end

function T = addIntervals(T,boot)
prefix=["d0";"full";"aligned_context";"aligned_plus"];
for j=1:height(T)
    r=prctile(boot.(prefix(j)+"_r2"),[2.5 97.5]); q=prctile(boot.(prefix(j)+"_rmse"),[2.5 97.5]); a=prctile(boot.(prefix(j)+"_mae"),[2.5 97.5]);
    T.r2_ci95_low(j)=r(1); T.r2_ci95_high(j)=r(2); T.rmse_ci95_low(j)=q(1); T.rmse_ci95_high(j)=q(2); T.mae_ci95_low(j)=a(1); T.mae_ci95_high(j)=a(2);
end
end

function T = pairedSummary(y,ctx,plus,eid,boot)
m1=metric(y,ctx); m2=metric(y,plus); q=prctile(boot.delta_r2,[2.5 97.5]); T=table(m2.r2-m1.r2,q(1),q(2),m2.rmse-m1.rmse,m2.mae-m1.mae, ...
    'VariableNames',{'delta_r2_loeo_cv','ci95_low','ci95_high','delta_rmse','delta_mae'});
end

function writePredictionOutputs(d0,d0p,full,fullp,aligned,ctx,plus,runDir)
A=d0(:,["record_id","eid","outage0"]); A.predicted_d0=d0p; writetable(A,fullfile(runDir,"models","d0_oof_predictions.csv"));
B=full(:,["record_id","eid","T90"]); B.predicted_context_log_t90=fullp; B.predicted_context_t90=exp(fullp); writetable(B,fullfile(runDir,"models","t90_full_context_oof_predictions.csv"));
C=aligned(:,["record_id","eid","T90","outage0"]); C.predicted_context_log_t90=ctx; C.predicted_plus_d0_log_t90=plus; C.residual_plus_d0=log(C.T90)-plus; writetable(C,fullfile(runDir,"models","t90_aligned_oof_predictions.csv"));
end

function R = predictionRobustness(T,cfg)
rows=cell(0,11); primary=T(T.t90_primary_eligible,:); zero=T(T.t90_zero_inclusive_eligible,:); D0=T(T.d0_model_eligible,:);
try, [frac,~]=d0Oof(D0,cfg); rows(end+1,:)=d0RobustnessRow("fractional_logit",D0,frac,"complete",""); catch ME, rows(end+1,:)=failedPredictionRow("d0_model","fractional_logit",height(D0),numel(unique(D0.eid)),"D0 original scale",ME.message); end
try, [zoi,~]=zoiOof(D0,cfg); rows(end+1,:)=d0RobustnessRow("zero_one_inflated_beta",D0,zoi,"complete",""); catch ME, rows(end+1,:)=failedPredictionRow("d0_model","zero_one_inflated_beta",height(D0),numel(unique(D0.eid)),"D0 original scale",ME.message); end
for f=["log2_d0" "raw_d0" "log2_d0_complete" "spline_log2_d0"]
 try, [c,p]=t90Oof(primary,cfg.primary_pga,f,"lognormal",nan,cfg); rows(end+1,:)=robRow("d0_representation",f,primary,c,p,cfg); catch ME, rows(end+1,:)=failedPredictionRow("d0_representation",f,height(primary),numel(unique(primary.eid)),"ln(T90)",ME.message); end
end
try, [c,p]=t90Oof(zero,cfg.primary_pga,"log2p1_d0","lognormal",nan,cfg); rows(end+1,:)=robRow("d0_representation","log2p1_d0_zero_inclusive",zero,c,p,cfg); catch ME, rows(end+1,:)=failedPredictionRow("d0_representation","log2p1_d0_zero_inclusive",height(zero),numel(unique(zero.eid)),"ln(T90)",ME.message); end
for outcome=["T80" "T95"]
 D=T(exactOutcomeDirect(T,outcome)&T.valid_context&T.valid_d0&T.outage0>0,:);
 try, [c,p]=t90OofOutcome(D,outcome,cfg.primary_pga,"log2_d0",cfg); rows(end+1,:)=robRowOutcome("restoration_milestone",outcome,D,c,p,outcome); catch ME, rows(end+1,:)=failedPredictionRow("restoration_milestone",outcome,height(D),numel(unique(D.eid)),"ln("+outcome+")",ME.message); end
end
for f=string(cfg.alternative_pga(:))'
 D=primary(isfinite(primary.(f))&primary.(f)>=0,:);
 try, [c,p]=t90Oof(D,f,"log2_d0","lognormal",nan,cfg); rows(end+1,:)=robRow("pga_summary",f,D,c,p,cfg); catch ME, rows(end+1,:)=failedPredictionRow("pga_summary",f,height(D),numel(unique(D.eid)),"ln(T90)",ME.message); end
end
for q=double(cfg.small_d0_exclusion_quantiles(:))'
 spec="exclude_bottom_"+string(round(100*q))+"pct";
 try, [c,p]=t90Oof(primary,cfg.primary_pga,"log2_d0","lognormal",q,cfg); rows(end+1,:)=robRow("small_positive_d0",spec,primary,c,p,cfg); catch ME, rows(end+1,:)=failedPredictionRow("small_positive_d0",spec,height(primary),numel(unique(primary.eid)),"ln(T90)",ME.message); end
end
for dist=["lognormal" "weibull" "loglogistic"]
 try, [c,p]=t90Oof(primary,cfg.primary_pga,"log2_d0",dist,nan,cfg); rows(end+1,:)=robRow("time_distribution",dist,primary,c,p,cfg); catch ME, rows(end+1,:)=failedPredictionRow("time_distribution",dist,height(primary),numel(unique(primary.eid)),"ln(T90)",ME.message); end
end
R=cell2table(rows,'VariableNames',{'family','specification','record_count','earthquake_count','primary_r2','context_r2','plus_d0_r2','delta_r2','outcome_scale','status','error_message'});
end
function row=d0RobustnessRow(spec,D,p,status,message), m=metric(D.outage0,p); row={"d0_model",string(spec),height(D),numel(unique(D.eid)),m.r2,nan,nan,nan,"D0 original scale",string(status),string(message)}; end
function row=failedPredictionRow(family,spec,n,e,scale,message), row={string(family),string(spec),n,e,nan,nan,nan,nan,string(scale),"failed",string(message)}; end

function row=robRow(family,spec,D,c,p,cfg)
row=robRowOutcome(family,spec,D,c,p,"T90");
end
function row=robRowOutcome(family,spec,D,c,p,outcome)
ok=isfinite(c)&isfinite(p); y=log(D.(outcome)(ok)); cm=metric(y,c(ok)); pm=metric(y,p(ok));
row={string(family),string(spec),sum(ok),numel(unique(D.eid(ok))),nan,cm.r2,pm.r2,pm.r2-cm.r2,"ln("+string(outcome)+")","complete",""};
end

function [c,p]=t90OofOutcome(D,outcome,pga,d0,cfg)
events=unique(D.eid,"stable"); c=nan(height(D),1);p=c;levels=incomeLevels(D,cfg);y=log(D.(outcome));
for k=1:numel(events), test=D.eid==events(k);train=~test;[xc,xct]=fixedDesign(D(train,:),D(test,:),levels,pga,cfg.primary_pga_transform,"none",false,[]);[xp,xpt]=fixedDesign(D(train,:),D(test,:),levels,pga,cfg.primary_pga_transform,d0,true,[]);bc=linearFit(xc,y(train));bp=linearFit(xp,y(train));c(test)=[ones(sum(test),1),xct]*bc;p(test)=[ones(sum(test),1),xpt]*bp;end
end

function [pred,fold]=zoiOof(D,cfg)
events=unique(D.eid,"stable");pred=nan(height(D),1);levels=incomeLevels(D,cfg);fold=table();
for k=1:numel(events)
    test=D.eid==events(k);train=~test;[Xtr,Xte]=fixedDesign(D(train,:),D(test,:),levels,cfg.primary_pga,cfg.primary_pga_transform,"none",false,[]);
    scaleMean=mean(Xtr,1); scaleSd=std(Xtr,0,1); scaleSd(scaleSd==0 | ~isfinite(scaleSd))=1; Ztr=(Xtr-scaleMean)./scaleSd; Zte=(Xte-scaleMean)./scaleSd;
    state=2*ones(sum(train),1);state(D.outage0(train)==0)=1;state(D.outage0(train)==1)=3; B=mnrfit(Ztr,state,'model','nominal'); P=mnrval(B,Zte,'model','nominal');
    partial=train&D.outage0>0&D.outage0<1; theta=fitBeta(Ztr(D.outage0(train)>0&D.outage0(train)<1,:),D.outage0(partial)); mu=logistic([ones(sum(test),1),Zte]*theta(1:end-1)); pred(test)=P(:,3)+P(:,2).*mu;
end
end

function directP3Primary(cfg,runDir)
if isfile(fullfile(runDir,"models","p3_fixed_effect_primary.csv")), fprintf("Primary association stage already complete; reusing saved outputs.\n"); return; end
S=load(fullfile(runDir,"input","canonical_analysis_table.mat"),"T"); primary=S.T(S.T.t90_primary_eligible,:);
[model,D,stats,coef,blocks]=fitP3(primary,"T90",cfg.primary_pga,"population","log2_d0",nan,cfg);
writetable(coef,fullfile(runDir,"models","p3_fixed_effect_primary.csv")); writetable(blocks,fullfile(runDir,"models","p3_joint_tests.csv"));
save(fullfile(runDir,"models","p3_primary_model.mat"),"model","D","stats","-v7.3");
writeStageStatus(runDir,"p3_primary","complete",struct("record_count",height(D),"earthquake_count",numel(unique(D.eid)),"marginal_r2",stats.marginal_r2,"conditional_r2",stats.conditional_r2));
end
function directP3R2(cfg,runDir)
r2Path=fullfile(runDir,"models","p3_r2.csv");
if isfile(r2Path), R=readtable(r2Path); if height(R)==2&&all(R.successful_replicates>=cfg.bootstrap_replicates), fprintf("P3 R2 bootstrap stage already complete; reusing saved outputs.\n"); return; end, end
S=load(fullfile(runDir,"input","canonical_analysis_table.mat"),"T"); primary=S.T(S.T.t90_primary_eligible,:); [model,D]=fitP3(primary,"T90",cfg.primary_pga,"population","log2_d0",nan,cfg);
R=parametricR2(model,D,cfg,fullfile(runDir,"validation","p3_r2_parametric")); writetable(R,r2Path);
writeStageStatus(runDir,"p3_r2","complete",struct("requested_replicates",cfg.bootstrap_replicates,"successful_replicates",R.successful_replicates(1)));
end
function directP3Robustness(cfg,runDir)
if isfile(fullfile(runDir,"validation","p3_association_robustness.csv")), fprintf("P3 robustness stage already complete; reusing saved outputs.\n"); return; end
S=load(fullfile(runDir,"input","canonical_analysis_table.mat"),"T"); rob=p3Robustness(S.T,cfg); writetable(rob,fullfile(runDir,"validation","p3_association_robustness.csv"));
principal=["sys_cat_W";"sys_cat_G";"pga_x";"year_x";"context_x";"d0_x"];
figureSy=rob(rob.status=="complete"&ismember(rob.term,principal),:);
writetable(figureSy,fullfile(runDir,"figure_data","supplementary_figure_sy_data.csv"));
writeStageStatus(runDir,"p3_robustness","complete",struct("row_count",height(rob),"failed_scenarios",numel(unique(rob.specification(rob.status~="complete")))));
end
function directP3Diagnostics(cfg,runDir)
if isfile(fullfile(runDir,"validation","p3_diagnostic_summary.json")), fprintf("P3 diagnostics stage already complete; reusing saved outputs.\n"); return; end
S=load(fullfile(runDir,"input","canonical_analysis_table.mat"),"T"); primary=S.T(S.T.t90_primary_eligible,:); [model,D]=fitP3(primary,"T90",cfg.primary_pga,"population","log2_d0",nan,cfg); diag=diagnostics(model,D,cfg);
writetable(diag.residual,fullfile(runDir,"figure_data","p3_residual_diagnostics.csv")); writetable(diag.random,fullfile(runDir,"figure_data","p3_random_intercepts.csv")); writetable(diag.influence,fullfile(runDir,"figure_data","p3_earthquake_influence.csv"));
writeJson(fullfile(runDir,"validation","p3_diagnostic_summary.json"),diag.summary);
writetable(struct2table(diag.summary),fullfile(runDir,"figure_data","supplementary_figure_sz_summary.csv"));
writeStageStatus(runDir,"p3_diagnostics","complete",diag.summary);
end

function [model,D,stats,coef,block]=fitP3(D,outcome,pgaField,contextKind,d0Kind,smallQuantile,cfg)
if isfinite(smallQuantile), cutoff=quantile(D.outage0,smallQuantile); D=D(D.outage0>cutoff,:); end
D=D(isfinite(D.(outcome))&D.(outcome)>0&isfinite(D.(pgaField))&D.(pgaField)>=0&isfinite(D.population)&D.population>0,:);
if contextKind=="density", D=D(isfinite(D.popDensity_km2)&D.popDensity_km2>0,:); end
D.log_outcome=log(D.(outcome)); D.pga_x=ctx_pga_transform(D.(pgaField),cfg.primary_pga_transform); D.year_x=(D.year-2000)/10;
if contextKind=="population",D.context_x=log2(D.population);else,D.context_x=log2(D.popDensity_km2);end
switch d0Kind
    case "log2_d0",D.d0_x=log2(D.outage0); addComplete=false;
    case "raw_d0",D.d0_x=D.outage0; addComplete=false;
    case "log2p1_d0",D.d0_x=log2(1+D.outage0); addComplete=false;
    case "log2_d0_complete",D.d0_x=log2(D.outage0);D.complete_x=double(D.outage0==1);addComplete=true;
    otherwise,error("ContextualDirect:UnknownP3D0","Unknown D0 encoding.");
end
D.sys_cat=categorical(D.sys,["E";"W";"G"]);D.sys_cat=reordercats(D.sys_cat,cellstr(["E";"W";"G"])); levels=incomeLevels(D,cfg);D.income_cat=categorical(D.development_level,levels);D.income_cat=reordercats(D.income_cat,cellstr(levels));D.event_group=categorical(D.eid);
formula="log_outcome ~ 1 + sys_cat + income_cat + pga_x + year_x + context_x + d0_x";if addComplete,formula=formula+" + complete_x";end;formula=formula+" + (1|event_group)";model=fitlme(D,char(formula),"FitMethod","ML");
fixed=predict(model,D,"Conditional",false);[psi,mse]=covarianceParameters(model);stats=struct('fixed_variance',var(fixed,1),'event_variance',psi{1}(1,1),'residual_variance',mse);total=stats.fixed_variance+stats.event_variance+stats.residual_variance;stats.marginal_r2=stats.fixed_variance/total;stats.conditional_r2=(stats.fixed_variance+stats.event_variance)/total;
coef=coefficientTable(model,D,contextKind,d0Kind,cfg.primary_pga_transform,pgaField); block=jointTests(model);
end

function C=coefficientTable(model,D,contextKind,d0Kind,pgaTransform,pgaField)
q=model.Coefficients;n=height(q);ci=coefCI(model,"Alpha",.05,"DFMethod","residual");p=nan(n,1);df=nan(n,1);for j=1:n,H=zeros(1,n);H(j)=1;[p(j),~,~,df(j)]=coefTest(model,H,0,"DFMethod","residual");end
factor=exp(q.Estimate);lo=exp(ci(:,1));hi=exp(ci(:,2));term=string(q.Name);label=term;label(term=="sys_cat_W")="Water versus electric power";label(term=="sys_cat_G")="Natural gas versus electric power";[pgaDelta,pgaLabel]=ctx_pga_contrast(pgaTransform); if string(pgaTransform)=="raw", pgaLabel=pgaSummaryContrastLabel(pgaField); end; pgaTerm=term=="pga_x";factor(pgaTerm)=exp(q.Estimate(pgaTerm)*pgaDelta);lo(pgaTerm)=exp(ci(pgaTerm,1)*pgaDelta);hi(pgaTerm)=exp(ci(pgaTerm,2)*pgaDelta);label(pgaTerm)=pgaLabel;label(term=="year_x")="Ten-year earthquake year";if contextKind=="population",label(term=="context_x")="Twofold matched service-region population";else,label(term=="context_x")="Twofold population density";end
if d0Kind=="log2_d0",label(term=="d0_x")="Twofold positive D0";elseif d0Kind=="log2p1_d0",label(term=="d0_x")="log2(1 + D0)";else,label(term=="d0_x")=d0Kind;end
label(startsWith(term,"income_cat_"))=replace(erase(term(startsWith(term,"income_cat_")),"income_cat_"),"_"," ")+" versus high income";
C=table(term,label,q.Estimate,q.SE,df,ci(:,1),ci(:,2),p,factor,lo,hi,'VariableNames',{'term','reported_contrast','estimate_log_t90','se','residual_df','ci95_low_log','ci95_high_log','p_value','adjusted_t90_factor','factor_ci95_low','factor_ci95_high'});
end

function B=jointTests(model)
names=string(model.Coefficients.Name);sets={find(startsWith(names,"sys_cat_")),find(startsWith(names,"income_cat_"))};blocks=["lifeline_system";"service_region_income_group"];F=nan(2,1);d1=F;d2=F;p=F;
for j=1:2,H=zeros(numel(sets{j}),numel(names));for k=1:numel(sets{j}),H(k,sets{j}(k))=1;end;[p(j),F(j),d1(j),d2(j)]=coefTest(model,H,zeros(numel(sets{j}),1),"DFMethod","residual");end
B=table(blocks,F,d1,d2,p,repmat("residual-degrees-of-freedom Wald F test",2,1),'VariableNames',{'block','f_statistic','numerator_df','denominator_df','p_value','inference_method'});
end

function R=parametricR2(model,D,cfg,outDir)
n=cfg.bootstrap_replicates; checkpoint=fullfile(outDir,"checkpoint.mat"); failurePath=fullfile(outDir,"bootstrap_failures.csv");
if isfolder(outDir)
    if ~isfile(checkpoint), error("ContextualDirect:R2Resume","Existing bootstrap directory has no checkpoint: %s",outDir); end
    S=load(checkpoint,"values","success","error_message");
    values=S.values; success=S.success; error_message=string(S.error_message);
    if size(values,1)~=n || size(values,2)~=2, error("ContextualDirect:R2Checkpoint","Checkpoint dimensions do not match requested replicates."); end
else
    mkdir(outDir); values=nan(n,2); success=false(n,1); error_message=strings(n,1);
end
[psi,mse]=covarianceParameters(model); ev=psi{1}(1,1); fixed=predict(model,D,"Conditional",false); G=findgroups(D.eid);
for b=1:n
    if success(b), continue; end
    % A per-replicate deterministic seed makes checkpoint resume exact.
    rng(cfg.random_seed+200+b,"twister");
    try
        B=D; eventNoise=sqrt(ev)*randn(max(G),1); B.log_outcome=fixed+eventNoise(G)+sqrt(mse)*randn(height(D),1);
        m=fitlme(B,char(model.Formula),"FitMethod","ML"); fp=predict(m,B,"Conditional",false); [p2,mse2]=covarianceParameters(m);
        vf=var(fp,1); total=vf+p2{1}(1,1)+mse2; values(b,:)=[vf/total,(vf+p2{1}(1,1))/total]; success(b)=true; error_message(b)="";
    catch ME
        error_message(b)=string(ME.message);
    end
    if mod(b,cfg.bootstrap_checkpoint_interval)==0||b==n
        save(checkpoint,"values","success","error_message","-v7.3");
        writetable(table((1:n)',success,error_message,'VariableNames',{'replicate','success','error_message'}),failurePath);
    end
end
fp=predict(model,D,"Conditional",false); vf=var(fp,1); total=vf+ev+mse; point=[vf/total,(vf+ev)/total];
if sum(success)<ceil(cfg.bootstrap_min_success_fraction*n), error("ContextualDirect:R2Failure","Only %d of %d bootstrap refits succeeded.",sum(success),n); end
lo=prctile(values(success,:),2.5); hi=prctile(values(success,:),97.5);
R=table(["marginal_r2";"conditional_r2"],point',lo',hi',repmat(sum(success),2,1),repmat(n,2,1),'VariableNames',{'quantity','estimate','ci95_low','ci95_high','successful_replicates','requested_replicates'});
writetable(array2table(values,'VariableNames',{'marginal_r2','conditional_r2'}),fullfile(outDir,"bootstrap_replicates.csv"));
end
function R=p3Robustness(T,cfg)
rows=cell(0,16); primary=T(T.t90_primary_eligible,:); scenarios=cell(0,7);
scenarios(end+1,:)={"pga_summary","primary_service_region_mean","T90",cfg.primary_pga,"population","log2_d0",nan};
for pga=string(cfg.alternative_pga(:))', scenarios(end+1,:)={"pga_summary",pga,"T90",pga,"population","log2_d0",nan}; end
scenarios(end+1,:)={"pga_function",cfg.pga_transform_sensitivity,"T90",cfg.primary_pga,"population","log2_d0",nan};
for outcome=["T80" "T95"], scenarios(end+1,:)={"restoration_milestone",outcome,outcome,cfg.primary_pga,"population","log2_d0",nan}; end
scenarios=[scenarios;{"d0_representation","raw_d0","T90",cfg.primary_pga,"population","raw_d0",nan};{"d0_representation","log2_d0_complete","T90",cfg.primary_pga,"population","log2_d0_complete",nan};{"d0_representation","log2p1_d0_zero_inclusive","T90",cfg.primary_pga,"population","log2p1_d0",-1};{"small_positive_d0","exclude_bottom_1pct","T90",cfg.primary_pga,"population","log2_d0",.01};{"small_positive_d0","exclude_bottom_5pct","T90",cfg.primary_pga,"population","log2_d0",.05};{"population_measure","population_density","T90",cfg.primary_pga,"density","log2_d0",nan}];
for k=1:size(scenarios,1)
 fam=string(scenarios{k,1}); id=string(scenarios{k,2}); outcome=string(scenarios{k,3}); pga=string(scenarios{k,4}); context=string(scenarios{k,5}); d0=string(scenarios{k,6}); q=scenarios{k,7}; D=primary; cfgScenario=cfg; if fam=="pga_function", cfgScenario.primary_pga_transform=cfg.pga_transform_sensitivity; end
 if outcome~="T90", D=T(exactOutcomeDirect(T,outcome)&T.valid_context&T.valid_d0&T.outage0>0,:); end
 if d0=="log2p1_d0", D=T(T.t90_zero_inclusive_eligible,:); end
 if context=="density", D=D(isfinite(D.popDensity_km2)&D.popDensity_km2>0,:); end
 if isfinite(q)&&q<0, q=nan; end
 try
  [~,Dfit,~,C,~]=fitP3(D,outcome,pga,context,d0,q,cfgScenario);
  for j=1:height(C), rows(end+1,:)={fam,id,outcome,pga,context,d0,height(Dfit),numel(unique(Dfit.eid)),C.term(j),C.reported_contrast(j),C.adjusted_t90_factor(j),C.factor_ci95_low(j),C.factor_ci95_high(j),C.p_value(j),"complete",""}; end
 catch ME
  rows(end+1,:)={fam,id,outcome,pga,context,d0,height(D),numel(unique(D.eid)),"","",nan,nan,nan,nan,"failed",string(ME.message)};
 end
end
R=cell2table(rows,'VariableNames',{'family','specification','outcome','pga_field','context_measure','d0_encoding','record_count','earthquake_count','term','reported_contrast','adjusted_t90_factor','ci95_low','ci95_high','p_value','status','error_message'});
end

function out=diagnostics(model,D,cfg)
res=double(residuals(model,"Conditional",true)); fit=double(fitted(model,"Conditional",true)); n=height(D); q=((1:n)'-.5)/n;
R=table(fit,res,sort(res),norminv(q),'VariableNames',{'fitted_log_t90','conditional_residual','residual_qq_observed','residual_qq_normal'});
rv=double(randomEffects(model)); rv=rv(:); m=numel(rv); qr=((1:m)'-.5)/m; random=table(sort(rv),norminv(qr),'VariableNames',{'random_intercept_sorted','normal_quantile'});
terms=["pga_x";"year_x";"context_x";"d0_x";"sys_cat_W";"sys_cat_G"]; rows=cell(0,7); events=unique(D.eid,"stable"); baseNames=string(model.Coefficients.Name);
for e=events'
 try
  m2=fitlme(D(D.eid~=e,:),char(model.Formula),"FitMethod","ML"); names=string(m2.Coefficients.Name);
  for j=1:numel(terms), base=model.Coefficients.Estimate(baseNames==terms(j)); value=m2.Coefficients.Estimate(names==terms(j)); rows(end+1,:)={e,terms(j),base,value,value-base,"complete",""}; end
 catch ME
  rows(end+1,:)={e,"",nan,nan,nan,"failed",string(ME.message)};
 end
end
influence=cell2table(rows,'VariableNames',{'omitted_earthquake','term','full_estimate_log_t90','omit_one_estimate_log_t90','difference_log_t90','status','error_message'});
X=[D.pga_x,D.year_x,D.context_x,D.d0_x]; vif=nan(4,1);
for j=1:4, A=[ones(height(D),1),X(:,setdiff(1:4,j))]; b=A\X(:,j); r2=1-sum((X(:,j)-A*b).^2)/sum((X(:,j)-mean(X(:,j))).^2); vif(j)=1/(1-r2); end
[psi,mse]=covarianceParameters(model); eventVariance=psi{1}(1,1); coeff=double(model.Coefficients.Estimate); logLikelihood=model.LogLikelihood;
converged=isfinite(logLikelihood)&&isfinite(mse)&&mse>0&&all(isfinite(coeff))&&isfinite(eventVariance)&&eventVariance>=0;
singularFit=eventVariance<=max(1e-10,1e-8*mse);
summary=struct('converged',converged,'singular_fit',singularFit,'convergence_check','finite log-likelihood, coefficients and variance components','singularity_check','random-intercept variance relative to residual variance','event_variance',eventVariance,'residual_variance',mse,'max_vif',max(vif),'residual_mean',mean(res),'residual_sd',std(res),'residual_skewness',skewness(res),'residual_kurtosis',kurtosis(res),'random_intercept_sd',std(rv),'random_intercept_skewness',skewness(rv),'leave_one_earthquake_failures',sum(influence.status~="complete"),'diagnostic_note','Diagnostics summarize the fitted primary model; they do not establish causal validity.');
out=struct('residual',R,'random',random,'influence',influence,'summary',summary);
end

function directFigures(cfg,runDir)
P=readtable(fullfile(runDir,"validation","prediction_robustness_unified.csv"),"TextType","string");A=readtable(fullfile(runDir,"validation","p3_association_robustness.csv"),"TextType","string");R=readtable(fullfile(runDir,"figure_data","p3_residual_diagnostics.csv"));Q=readtable(fullfile(runDir,"figure_data","p3_random_intercepts.csv"));I=readtable(fullfile(runDir,"figure_data","p3_earthquake_influence.csv"),"TextType","string");
f=figure('Visible','off','Color','w','Position',[100 100 1500 760]);tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
families=["d0_model" "d0_representation" "restoration_milestone" "pga_summary" "small_positive_d0" "time_distribution"];
titles=["a  Initial-disruption model" "b  Initial-disruption representation" "c  Restoration milestone" "d  Regional PGA summary" "e  Small positive D0" "f  Restoration-time distribution"];
for kk=1:numel(families)
 nexttile; U=P(P.family==families(kk)&P.status=="complete",:);
 if families(kk)=="d0_model"
  bar(categorical(replace(U.specification,"_"," ")),U.primary_r2); ylabel('LOEO-CV R^2'); ylim([0 max(.3,max(U.primary_r2)*1.15)]);
 else
  bar(categorical(replace(U.specification,"_"," ")),U.delta_r2); yline(0,'k--'); ylabel('\Delta R^2_{LOEO-CV}');
 end
 title(titles(kk)); xtickangle(30);
end
savefig(f,fullfile(runDir,'figures','Supplementary_Figure_contextual_prediction_robustness.fig'));exportgraphics(f,fullfile(runDir,'figures','Supplementary_Figure_contextual_prediction_robustness.png'),'Resolution',cfg.figure_dpi);close(f);
f=figure('Visible','off','Color','w','Position',[100 100 1300 900]);tiledlayout(2,2,'Padding','compact','TileSpacing','compact');families=["pga_summary" "restoration_milestone" "d0_representation" "population_measure"];titles=["a  PGA specification" "b  Restoration milestone" "c  Initial-disruption specification" "d  Population measure"];for k=1:4,nexttile;D=A(A.family==families(k)&A.status=="complete"&~startsWith(A.term,"income" )&A.term~="(Intercept)",:);plotAssociationPanel(D);title(titles(k));end;savefig(f,fullfile(runDir,'figures','Supplementary_Figure_adjusted_association_robustness.fig'));exportgraphics(f,fullfile(runDir,'figures','Supplementary_Figure_adjusted_association_robustness.png'),'Resolution',cfg.figure_dpi);close(f);
f=figure('Visible','off','Color','w','Position',[100 100 1300 900]);tiledlayout(2,2,'Padding','compact','TileSpacing','compact');nexttile;scatter(R.fitted_log_t90,R.conditional_residual,10,'filled');yline(0,'k--');xlabel('Fitted ln(T_{90})');ylabel('Conditional residual');title('a  Residuals versus fitted');nexttile;plot(R.residual_qq_normal,R.residual_qq_observed,'.');refline(1,0);xlabel('Normal quantile');ylabel('Residual quantile');title('b  Residual Q-Q');nexttile;plot(Q.normal_quantile,Q.random_intercept_sorted,'.');refline(1,0);xlabel('Normal quantile');ylabel('Random-intercept quantile');title('c  Random-intercept Q-Q');nexttile;terms=unique(I.term,'stable');hold on;for j=1:numel(terms),D=I(I.term==terms(j)&I.status=="complete",:);plot(j+0*D.difference_log_t90,D.difference_log_t90,'.');end;yline(0,'k--');xlim([.5 numel(terms)+.5]);xticks(1:numel(terms));xticklabels(strrep(terms,'_',' '));xtickangle(25);ylabel('Change after omitting one earthquake');title('d  Earthquake-level influence');savefig(f,fullfile(runDir,'figures','Supplementary_Figure_mixed_model_diagnostics.fig'));exportgraphics(f,fullfile(runDir,'figures','Supplementary_Figure_mixed_model_diagnostics.png'),'Resolution',cfg.figure_dpi);close(f);
end

function plotAssociationPanel(D)
if isempty(D),axis off;return;end;terms=unique(D.term,'stable');hold on;for j=1:numel(terms),S=D(D.term==terms(j),:);x=1:height(S);errorbar(x,S.adjusted_t90_factor,S.adjusted_t90_factor-S.ci95_low,S.ci95_high-S.adjusted_t90_factor,'o');end;yline(1,'k--');set(gca,'YScale','log');ylabel('Adjusted T_{90} factor');xlabel('Specification');
end

function directValidate(cfg,runDir)
required=[fullfile(runDir,"models","primary_prediction_performance.csv"),fullfile(runDir,"validation","prediction_robustness_unified.csv"),fullfile(runDir,"models","p3_fixed_effect_primary.csv"),fullfile(runDir,"models","p3_r2.csv"),fullfile(runDir,"validation","p3_association_robustness.csv"),fullfile(runDir,"validation","p3_diagnostic_summary.json"),fullfile(runDir,"figures","Supplementary_Figure_contextual_prediction_robustness.fig"),fullfile(runDir,"figures","Supplementary_Figure_adjusted_association_robustness.fig"),fullfile(runDir,"figures","Supplementary_Figure_mixed_model_diagnostics.fig")];
present=cellfun(@isfile,cellstr(required)); status="complete"; if ~all(present), status="failed"; end
details=struct('required_files',{cellstr(required)},'present',present,'p3_r2_replicates_ok',false,'validated_at',char(datetime('now','Format',"yyyy-MM-dd'T'HH:mm:ss")));
if isfile(fullfile(runDir,"models","p3_r2.csv")), r2=readtable(fullfile(runDir,"models","p3_r2.csv")); details.p3_r2_replicates_ok=all(r2.successful_replicates>=cfg.bootstrap_replicates); if ~details.p3_r2_replicates_ok, status="failed"; end, end
details.status=char(status); writeJson(fullfile(runDir,"manifest","validation.json"),details); writeStageStatus(runDir,"validate",status,details);
if status~="complete", error("ContextualDirect:Validation","Validation failed; inspect manifest/validation.json."); end
end
function writeStageStatus(runDir,stage,status,details), payload=struct('stage',char(stage),'status',char(status),'timestamp',char(datetime('now','Format',"yyyy-MM-dd'T'HH:mm:ss")),'details',details); writeJson(fullfile(runDir,"manifest","stage_"+string(stage)+".json"),payload); end
function writeRunManifest(cfg,runDir,status)
sourcePath=fullfile(cfg.project_root,cfg.source_xlsx); runnerPath=string(mfilename("fullpath"))+".m"; configPath=fullfile(runDir,"manifest","resolved_config.json");
payload=struct('analysis_id',cfg.analysis_id,'package_version',cfg.package_version,'run_id',cfg.run_id,'status',char(status),'timestamp',char(datetime('now','Format',"yyyy-MM-dd'T'HH:mm:ss")),'source_xlsx',cfg.source_xlsx,'source_xlsx_sha256',fileHash(sourcePath),'runner_sha256',fileHash(runnerPath),'config_sha256',fileHash(configPath)); writeJson(fullfile(runDir,"manifest","run_metadata.json"),payload);
end
function hash=fileHash(path)
fid=fopen(path,'r'); if fid<0, error("ContextualDirect:HashRead","Unable to read %s",path); end
cleanup=onCleanup(@()fclose(fid)); bytes=fread(fid,Inf,'*uint8')'; digest=java.security.MessageDigest.getInstance('SHA-256'); digest.update(bytes); hash=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'))',1,[]));
end

function label=pgaSummaryContrastLabel(pgaField)
switch string(pgaField)
    case "avg_PGA", label="0.1 g service-region mean PGA";
    case "popWeightedPGA_g", label="0.1 g population-weighted PGA";
    case "popMeanPGA_g", label="0.1 g population-mean PGA";
    case "pgaMedian_g", label="0.1 g median PGA";
    case "pgaMax_g", label="0.1 g maximum PGA";
    otherwise, error("ContextualAnalysis:UnknownPGASummary","Unknown PGA summary: %s",string(pgaField));
end
end
function levels=incomeLevels(D,cfg), levels=sort(unique(string(D.development_level)));reference=string(cfg.reference_income_group);if ~ismember(reference,levels),error("ContextualDirect:IncomeReference","High income reference absent.");end;levels=[reference;levels(levels~=reference)];end
function m=metric(y,p),m=struct('r2',1-sum((y-p).^2)/sum((y-mean(y)).^2),'rmse',sqrt(mean((y-p).^2)),'mae',mean(abs(y-p)));end
function idx=sampleEventRows(eid),events=unique(eid,'stable');draw=events(randi(numel(events),numel(events),1));idx=[];for k=1:numel(draw),idx=[idx;find(eid==draw(k))];end,end
function b=linearFit(X,y),A=[ones(size(X,1),1),X];if rank(A)<size(A,2),b=lsqminnorm(A,y);else,b=A\y;end,end
function B=rcsBasis(x,k),k=unique(k);if numel(k)<4,error("ContextualDirect:Spline","Insufficient distinct spline knots.");end;last=k(end);penult=k(end-1);B=x;for j=1:numel(k)-2,h=max(x-k(j),0).^3-((last-k(j))/(last-penult))*max(x-penult,0).^3+((penult-k(j))/(last-penult))*max(x-last,0).^3;B(:,end+1)=h/(last-k(1))^2;end,end
function value=exactOutcomeDirect(T,outcome),c="C"+extractAfter(outcome,"T");b=outcome+"_bound_type";value=T.valid_system&isfinite(T.(outcome))&T.(outcome)>0&T.(c)==0&lower(string(T.(b)))=="exact";end
function theta=fitBeta(X,y),A=[ones(size(X,1),1),X];init=zeros(size(A,2)+1,1);init(1)=log(mean(y)/(1-mean(y)));init(end)=log(10);o=optimoptions('fminunc','Display','off','Algorithm','quasi-newton','MaxIterations',1000);[theta,~,flag]=fminunc(@(v)betaNll(v,A,y),init,o);if flag<=0,error("ContextualDirect:Beta","Beta model did not converge.");end,end
function v=betaNll(theta,A,y),mu=logistic(A*theta(1:end-1));phi=exp(theta(end));a=max(mu*phi,eps);b=max((1-mu)*phi,eps);v=-sum(gammaln(phi)-gammaln(a)-gammaln(b)+(a-1).*log(y)+(b-1).*log1p(-y));if ~isfinite(v),v=realmax;end,end
function y=logistic(x),y=1./(1+exp(-max(min(x,35),-35)));end
function theta=fitAft(X,y,dist),A=[ones(size(X,1),1),X];b=linearFit(X,y);init=[b;log(max(std(y-A*b),.05))];o=optimoptions('fminunc','Display','off','Algorithm','quasi-newton','MaxIterations',500);[theta,~,flag]=fminunc(@(v)aftNll(v,A,y,dist),init,o);if flag<=0,error("ContextualDirect:AFT","AFT model did not converge.");end,end
function v=aftNll(theta,A,y,dist),mu=A*theta(1:end-1);s=exp(theta(end));z=(y-mu)/s;if dist=="weibull",lp=-log(s)+z-exp(z);else,lp=-log(s)-z-2*softplus(-z);end;v=-sum(lp);end
function y=aftMedian(X,t,dist),mu=[ones(size(X,1),1),X]*t(1:end-1);s=exp(t(end));if dist=="weibull",y=mu+s*log(log(2));else,y=mu;end,end
function y=softplus(x),y=max(x,0)+log1p(exp(-abs(x)));end
function writeJson(path,value),fid=fopen(path,'w');cleanup=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(value,'PrettyPrint',true),'char');end









