function run_contextual_supplementary(projectRoot, inputPath, configPath, primaryDir, runDir, nReplicates)
if nargin < 1 || isempty(projectRoot) || strlength(string(projectRoot)) == 0
    moduleRoot = ctx_module_root();
    projectRoot = string(fileparts(fileparts(fileparts(moduleRoot))));
end
if nargin < 5, error("ContextualSupplementary:ArgumentsRequired", "Five path arguments are required."); end
if nargin < 6 || isempty(nReplicates), nReplicates = 2000; end
projectRoot=string(projectRoot); inputPath=string(inputPath); configPath=string(configPath); primaryDir=string(primaryDir); runDir=string(runDir);
primaryPath=fullfile(primaryDir,"models","primary_fitted_performance.csv"); incrementPath=fullfile(primaryDir,"models","primary_fitted_increment.csv");
fullEffectPath=fullfile(primaryDir,"models","t90_full_record_set_effects.csv"); alignedEffectPath=fullfile(primaryDir,"models","t90_D0_aligned_plus_D0_effects.csv");
required=[inputPath;configPath;primaryPath;incrementPath;fullEffectPath;alignedEffectPath]; if any(~isfile(required)), error("ContextualSupplementary:MissingInput", "Main analysis input is incomplete."); end
for d=["figure_data" "validation" "manifest" "figures" "logs"], if ~isfolder(fullfile(runDir,d)), mkdir(fullfile(runDir,d)); end, end
cfg=jsondecode(fileread(configPath)); cfg.bootstrap_replicates=nReplicates; seed=double(cfg.random_seed)+9010; manifest=struct(); manifest.analysis_id="contextual_supplementary"; manifest.run_id=char(runDir); manifest.input_path=char(inputPath); manifest.primary_dir=char(primaryDir); manifest.bootstrap_replicates=nReplicates; manifest.random_seed=seed; manifest.status="started";
manifest.input_sha256=fileHash(inputPath);manifest.config_sha256=fileHash(configPath);
manifest.runner_sha256=fileHash(mfilename("fullpath")+".m");
manifest.inference_runner_sha256=fileHash(fullfile(fileparts(mfilename("fullpath")),"run_contextual_supp_inference.m"));
manifest.fixed_factor_inference="CR2 Satterthwaite";manifest.mixed_factor_inference="Wald Satterthwaite";
writeJson(fullfile(runDir,"manifest","run_started.json"),manifest); S=load(inputPath,"T"); T=S.T;
[d0Specs,fullSpecs,pairSpecs]=specifications(cfg);
[d0Point,d0Sets]=pointD0(T,d0Specs,cfg);
[fullPoint,fullSets]=pointFull(T,fullSpecs,cfg);
[pairPoint,pairSets]=pointPair(T,pairSpecs,cfg);
[systemResult,d0FactorResult]=run_contextual_supp_inference(T,cfg);
checkpoint=fullfile(runDir,"validation","bootstrap_checkpoint.mat");
[boot,done]=bootstrapAll(d0Sets,d0Specs,fullSets,fullSpecs,pairSets,pairSpecs,cfg,nReplicates,seed,checkpoint,inputPath,configPath);
d0Result=attachIntervals(d0Point,boot.d0_r2); fullResult=attachIntervals(fullPoint,boot.full_r2);
pairResult=attachIntervals(pairPoint,boot.pair_delta);
[d0Result,fullResult,pairResult]=harmonizePrimary(d0Result,fullResult,pairResult,primaryPath,incrementPath,nReplicates);
validateFactorInference(systemResult,d0FactorResult,fullEffectPath,alignedEffectPath);
writetable(d0Result,fullfile(runDir,"figure_data","s1a_D0_fitted_R2_sensitivity.csv"));
writetable(fullResult,fullfile(runDir,"figure_data","s1b_T90_fitted_R2_sensitivity.csv"));
writetable(pairResult,fullfile(runDir,"figure_data","s1c_paired_fitted_delta_R2_sensitivity.csv"));
writetable(systemResult,fullfile(runDir,"figure_data","s3a_system_factor_sensitivity.csv"));
writetable(d0FactorResult,fullfile(runDir,"figure_data","s3b_D0_factor_sensitivity.csv"));
[loeoPrimary,loeoBaseline]=run_contextual_loeo_cv(T,cfg); loeoRobust=ctx_unified_prediction_robustness(T,cfg,loeoBaseline);
writetable(loeoPrimary,fullfile(runDir,"figure_data","s2a_primary_LOEO_performance.csv"));
writetable(loeoRobust,fullfile(runDir,"validation","loeo_sensitivity_status.csv"));
L=loeoRobust(loeoRobust.panel=="c",:);
if any(L.status~="complete"),error("ContextualSupplementary:LOEOSensitivity","A paired LOEO sensitivity failed; see validation/loeo_sensitivity_status.csv.");end
writetable(L,fullfile(runDir,"figure_data","s2b_paired_LOEO_delta_sensitivity.csv"));
writeDiagnostics(T,cfg,runDir); writeBootstrapStatus(boot,done,nReplicates,runDir);
validatePrimary(d0Result,fullResult,pairResult,primaryPath,incrementPath);
manifest.status="complete"; manifest.completed_at=char(datetime("now","Format","yyyy-MM-dd'T'HH:mm:ss")); writeJson(fullfile(runDir,"manifest","run_complete.json"),manifest);
fprintf("Contextual supplementary analysis complete: %s",runDir);
end
function [d0Specs,fullSpecs,pairSpecs] = specifications(cfg)
d0Specs = repmat(baseSpec(),0,1);
d0Specs(end+1)=spec("primary","Primary fractional logit","outcome model","D0","avg_PGA","raw","population","none","lognormal",nan,"fractional");
d0Specs(end+1)=spec("zero_one_inflated_beta","Zero-one-inflated beta","outcome model","D0","avg_PGA","raw","population","none","lognormal",nan,"zoi_beta");
for f=string(cfg.alternative_pga(:))'
    d0Specs(end+1)=spec(f,pgaLabel(f),"PGA summary","D0",f,"raw","population","none","lognormal",nan,"fractional"); %#ok<AGROW>
end
d0Specs(end+1)=spec("ln1p_pga_over_0p1","ln[1 + PGA / (0.1 g)]","PGA functional form","D0","avg_PGA","ln1p","population","none","lognormal",nan,"fractional");
d0Specs(end+1)=spec("population_density","Population density","population measure","D0","avg_PGA","raw","density","none","lognormal",nan,"fractional");

fullSpecs = repmat(baseSpec(),0,1);
fullSpecs(end+1)=spec("primary","Primary model","primary","T90","avg_PGA","raw","population","none","lognormal",nan,"linear");
for f=string(cfg.alternative_pga(:))'
    fullSpecs(end+1)=spec(f,pgaLabel(f),"PGA summary","T90",f,"raw","population","none","lognormal",nan,"linear"); %#ok<AGROW>
end
fullSpecs(end+1)=spec("ln1p_pga_over_0p1","ln[1 + PGA / (0.1 g)]","PGA functional form","T90","avg_PGA","ln1p","population","none","lognormal",nan,"linear");
fullSpecs(end+1)=spec("population_density","Population density","population measure","T90","avg_PGA","raw","density","none","lognormal",nan,"linear");

pairSpecs = repmat(baseSpec(),0,1);
pairSpecs(end+1)=spec("primary","log2(D0) (primary)","D0 representation","T90","avg_PGA","raw","population","log2_d0","lognormal",nan,"linear");
pairSpecs(end+1)=spec("raw_d0","Raw D0","D0 representation","T90","avg_PGA","raw","population","raw_d0","lognormal",nan,"linear");
pairSpecs(end+1)=spec("log2_d0_complete","log2(D0) + complete-disruption indicator","D0 representation","T90","avg_PGA","raw","population","log2_d0_complete","lognormal",nan,"linear");
pairSpecs(end+1)=spec("spline_log2_d0","Spline in log2(D0)","D0 representation","T90","avg_PGA","raw","population","spline_log2_d0","lognormal",nan,"linear");
pairSpecs(end+1)=spec("log2p1_d0_zero_inclusive","log2(1 + D0)","zero-inclusive D0","T90","avg_PGA","raw","population","log2p1_d0","lognormal",nan,"linear");
for out=["T80" "T95"]
    pairSpecs(end+1)=spec(out,out,"restoration milestone",out,"avg_PGA","raw","population","log2_d0","lognormal",nan,"linear"); %#ok<AGROW>
end
for f=string(cfg.alternative_pga(:))'
    pairSpecs(end+1)=spec(f,pgaLabel(f),"PGA summary","T90",f,"raw","population","log2_d0","lognormal",nan,"linear"); %#ok<AGROW>
end
pairSpecs(end+1)=spec("ln1p_pga_over_0p1","ln[1 + PGA / (0.1 g)]","PGA functional form","T90","avg_PGA","ln1p","population","log2_d0","lognormal",nan,"linear");
pairSpecs(end+1)=spec("population_density","Population density","population measure","T90","avg_PGA","raw","density","log2_d0","lognormal",nan,"linear");
pairSpecs(end+1)=spec("exclude_bottom_1pct","Exclude bottom 1% of positive D0","small-positive D0","T90","avg_PGA","raw","population","log2_d0","lognormal",0.01,"linear");
pairSpecs(end+1)=spec("exclude_bottom_5pct","Exclude bottom 5% of positive D0","small-positive D0","T90","avg_PGA","raw","population","log2_d0","lognormal",0.05,"linear");

end

function s=baseSpec()
s=struct("id","","label","","family","","outcome","T90","pga_field","avg_PGA", ...
    "pga_transform","raw","context_kind","population","d0_kind","none", ...
    "distribution","lognormal","trim_quantile",nan,"model_kind","linear");
end
function s=spec(id,label,family,outcome,pga,pgaTransform,context,d0,dist,trim,model)
s=baseSpec(); s.id=string(id);s.label=string(label);s.family=string(family);s.outcome=string(outcome);
s.pga_field=string(pga);s.pga_transform=string(pgaTransform);s.context_kind=string(context);
s.d0_kind=string(d0);s.distribution=string(dist);s.trim_quantile=trim;s.model_kind=string(model);
end
function label=pgaLabel(f)
switch string(f)
    case "popWeightedPGA_g",label="Population-weighted PGA";
    case "popMeanPGA_g",label="Mean PGA across populated cells";
    case "pgaMedian_g",label="Median PGA across populated cells";
    case "pgaMax_g",label="Maximum PGA across populated cells";
    otherwise,label=replace(string(f),"_"," ");
end
end

function [R,sets]=pointD0(T,specs,cfg)
n=numel(specs); values=nan(n,1); records=zeros(n,1); events=zeros(n,1); sets=cell(n,1);
for j=1:n
    D=selectD0(T,specs(j)); levels=incomeLevels(D,cfg); p=fitD0(D,specs(j),levels);
    values(j)=metric(D.outage0,p).r2;records(j)=height(D);events(j)=numel(unique(D.eid));sets{j}=struct("D",D,"levels",levels);
end
R=specTable(specs,records,events);R.estimate=values;
end
function [R,sets]=pointFull(T,specs,cfg)
n=numel(specs); values=nan(n,1);records=zeros(n,1);events=zeros(n,1);sets=cell(n,1);
for j=1:n
    D=selectTime(T,specs(j),false);levels=incomeLevels(D,cfg);fit=fitTime(D,specs(j),false,levels);
    values(j)=metric(log(D.(specs(j).outcome)),fit.prediction).r2;records(j)=height(D);events(j)=numel(unique(D.eid));sets{j}=struct("D",D,"levels",levels);
end
R=specTable(specs,records,events);R.estimate=values;
end
function [R,sets]=pointPair(T,specs,cfg)
n=numel(specs);delta=nan(n,1);partial=nan(n,1);factor=nan(n,1);records=zeros(n,1);events=zeros(n,1);sets=cell(n,1);
for j=1:n
    D=selectTime(T,specs(j),true);levels=incomeLevels(D,cfg);a=fitTime(D,specs(j),false,levels);b=fitTime(D,specs(j),true,levels);
    y=log(D.(specs(j).outcome));m0=metric(y,a.prediction);m1=metric(y,b.prediction);
    delta(j)=m1.r2-m0.r2;partial(j)=(m0.sse-m1.sse)/m0.sse;factor(j)=d0Factor(D,specs(j),b,levels);
    records(j)=height(D);events(j)=numel(unique(D.eid));sets{j}=struct("D",D,"levels",levels);
end
R=specTable(specs,records,events);R.estimate=delta;R.partial_r2=partial;R.d0_factor=factor;
end


function R=specTable(specs,records,events)
R=table(string({specs.id})',string({specs.label})',string({specs.family})',records,events, ...
    'VariableNames',{'specification','label','family','record_count','earthquake_count'});
end

function [boot,done]=bootstrapAll(d0Sets,d0Specs,fullSets,fullSpecs,pairSets,pairSpecs,cfg,n,seed,checkpoint,inputPath,configPath)
names=["d0_r2","full_r2","pair_delta"];
sets={d0Sets,fullSets,pairSets}; specs={d0Specs,fullSpecs,pairSpecs};
identity=struct("schema",2,"input_sha256",fileHash(inputPath), ...
    "config_sha256",fileHash(configPath),"runner_sha256",fileHash(mfilename("fullpath")+".m"),"replicates",n,"seed",seed);
boot=struct(); topupAttempts=struct(); failureRows=cell(0,5); done=false(n,1);
for k=1:numel(names)
    boot.(names(k))=nan(n,numel(specs{k}));
    topupAttempts.(names(k))=zeros(1,numel(specs{k}));
end
if isfile(checkpoint)
    C=load(checkpoint);
    if ~isfield(C,"identity") || ~isequal(C.identity,identity)
        error("ContextualSupplementary:CheckpointIdentity","Checkpoint input, configuration or schema differs; use a new run directory.");
    end
    boot=C.boot;done=C.done;topupAttempts=C.topupAttempts;failureRows=C.failureRows;
end
blockSize=double(cfg.bootstrap_checkpoint_interval);
if n>=100 && any(~done) && isempty(gcp("nocreate")),parpool("threads");end
for first=1:blockSize:n
    last=min(first+blockSize-1,n);ids=(first:last)';ids=ids(~done(ids));
    if isempty(ids),continue;end
    values=cell(numel(ids),1); failures=cell(numel(ids),1);
    if n>=100
        parfor k=1:numel(ids)
            [values{k},failures{k}]=oneBootstrap(ids(k),seed,sets,specs,names);
        end
    else
        for k=1:numel(ids)
            [values{k},failures{k}]=oneBootstrap(ids(k),seed,sets,specs,names);
        end
    end
    for k=1:numel(ids)
        for q=1:numel(names),boot.(names(q))(ids(k),:)=values{k}.(names(q));end
        failureRows=[failureRows;failures{k}]; %#ok<AGROW>
        done(ids(k))=true;
    end
    save(checkpoint,"boot","done","identity","topupAttempts","failureRows","-v7.3");
end
% Complete each statistic independently; aligned comparisons share one draw.
maxTopupAttempts=max(1000,n);
for k=1:numel(names)
    for j=1:numel(specs{k})
        while any(~isfinite(boot.(names(k))(:,j)))
            if topupAttempts.(names(k))(j)>=maxTopupAttempts
                error("ContextualSupplementary:BootstrapIncomplete", ...
                    "Could not obtain %d valid resamples for %s/%s. Checkpoint and failures retained.",n,names(k),specs{k}(j).id);
            end
            topupAttempts.(names(k))(j)=topupAttempts.(names(k))(j)+1;
            attempt=topupAttempts.(names(k))(j);
            rng(seed+1000000*k+10000*j+attempt,"twister");
            try
                value=metricResample(names(k),sets{k}{j},specs{k}(j));
                row=find(~isfinite(boot.(names(k))(:,j)),1);
                boot.(names(k))(row,j)=value;
            catch ME
                failureRows(end+1,:)={names(k),specs{k}(j).id,n+attempt,string(ME.identifier),string(ME.message)}; %#ok<AGROW>
            end
            save(checkpoint,"boot","done","identity","topupAttempts","failureRows","-v7.3");
        end
    end
end
F=cell2table(failureRows,'VariableNames',{'quantity','specification','attempt','error_id','error_message'});
writetable(F,fullfile(fileparts(checkpoint),"bootstrap_failures.csv"));
end

function [v,failures]=oneBootstrap(b,seed,sets,specs,names)
rng(seed+b,"twister");v=struct();failures=cell(0,5);
for k=1:numel(names)
    v.(names(k))=nan(1,numel(specs{k}));
    for j=1:numel(specs{k})
        try
            v.(names(k))(j)=metricResample(names(k),sets{k}{j},specs{k}(j));
        catch ME
            failures(end+1,:)={names(k),specs{k}(j).id,b,string(ME.identifier),string(ME.message)}; %#ok<AGROW>
        end
    end
end
end

function value=metricResample(kind,set,s)
D=set.D(sampleEventRows(set.D.eid),:);
if kind=="d0_r2"
    prediction=fitD0(D,s,set.levels);value=metric(D.outage0,prediction).r2;
elseif kind=="full_r2"
    fit=fitTime(D,s,false,set.levels);value=metric(log(D.(s.outcome)),fit.prediction).r2;
else
    a=fitTime(D,s,false,set.levels);b=fitTime(D,s,true,set.levels);
    y=log(D.(s.outcome));m0=metric(y,a.prediction);m1=metric(y,b.prediction);
    value=m1.r2-m0.r2;
end
if ~isfinite(value),error("ContextualSupplementary:NonFiniteMetric","Non-finite %s for %s.",kind,s.id);end
end

function R=attachIntervals(R,B)
if any(~isfinite(B),"all")
    error("ContextualSupplementary:BootstrapIncomplete","Every statistic requires the requested number of valid resamples.");
end
R.ci95_low=prctile(B,2.5,1)';R.ci95_high=prctile(B,97.5,1)';
R.bootstrap_successes=repmat(size(B,1),height(R),1);
end

function [d0,full,pair]=harmonizePrimary(d0,full,pair,performancePath,incrementPath,n)
P=readtable(performancePath,"TextType","string");I=readtable(incrementPath);
statusPath=fullfile(fileparts(fileparts(performancePath)),"validation","bootstrap_status.csv");
S=readtable(statusPath);
if any(S.successful_replicates~=n)
    error("ContextualSupplementary:PrimaryBootstrapCount","Main and supplementary bootstrap counts must agree.");
end
row=d0.specification=="primary";src=P.model=="system_and_context_full_D0_record_set";d0.ci95_low(row)=P.r2_ci95_low(src);d0.ci95_high(row)=P.r2_ci95_high(src);
row=full.specification=="primary";src=P.model=="system_and_context_full_T90_record_set";full.ci95_low(row)=P.r2_ci95_low(src);full.ci95_high(row)=P.r2_ci95_high(src);
row=pair.specification=="primary";pair.ci95_low(row)=I.delta_r2_ci95_low;pair.ci95_high(row)=I.delta_r2_ci95_high;
end

function validateFactorInference(system,d0Factor,fullEffectPath,alignedEffectPath)
F=readtable(fullEffectPath,"TextType","string");A=readtable(alignedEffectPath,"TextType","string");
row=system.specification=="primary";w=F.effect_id=="water_vs_electric";g=F.effect_id=="gas_vs_electric";
actual=[system.water_factor(row),system.water_ci95_low(row),system.water_ci95_high(row),system.water_p_value(row), ...
    system.gas_factor(row),system.gas_ci95_low(row),system.gas_ci95_high(row),system.gas_p_value(row)];
expected=[F.estimate(w),F.ci95_low(w),F.ci95_high(w),F.p_value(w),F.estimate(g),F.ci95_low(g),F.ci95_high(g),F.p_value(g)];
assert(all(abs(actual-expected)<1e-8),"Primary system inference differs between main and supplementary fits.");
row=d0Factor.specification=="primary";src=A.effect_id=="D0_doubling";
actual=[d0Factor.estimate(row),d0Factor.ci95_low(row),d0Factor.ci95_high(row),d0Factor.p_value(row)];
expected=[A.estimate(src),A.ci95_low(src),A.ci95_high(src),A.p_value(src)];
assert(all(abs(actual-expected)<1e-8),"Primary D0 inference differs between main and supplementary fits.");
end

function D=selectD0(T,s)
D=T(T.d0_model_eligible,:);D=filterContext(D,s);
end
function D=selectTime(T,s,paired)
out=s.outcome;
if paired
    if s.d0_kind=="log2p1_d0",D=T(T.t90_zero_inclusive_eligible,:);
    elseif out=="T90",D=T(T.t90_primary_eligible,:);
    else,D=T(exactOutcome(T,out)&T.valid_context&T.valid_d0&T.outage0>0,:);end
else
    if out=="T90",D=T(T.t90_context_eligible,:);
    else,D=T(exactOutcome(T,out)&T.valid_context,:);end
end
D=filterContext(D,s);
if paired&&isfinite(s.trim_quantile)
    threshold=quantile(D.outage0,s.trim_quantile);D=D(D.outage0>threshold,:);
end
end
function D=filterContext(D,s)
ok=isfinite(D.(s.pga_field))&D.(s.pga_field)>=0;
if s.context_kind=="density",ok=ok&isfinite(D.popDensity_km2)&D.popDensity_km2>0;else,ok=ok&isfinite(D.population)&D.population>0;end
D=D(ok,:);
end
function tf=exactOutcome(T,out)
c="C"+extractAfter(out,"T");b=out+"_bound_type";tf=T.valid_system&isfinite(T.(out))&T.(out)>0&T.(c)==0&lower(string(T.(b)))=="exact";
end

function p=fitD0(D,s,levels)
X=baseDesign(D,s,levels);
if s.model_kind=="fractional"
    lastwarn("");m=fitglm(X,D.outage0,"linear","Distribution","binomial","Link","logit","Options",statset("MaxIter",5000));[~,wid]=lastwarn;if strlength(string(wid))>0,error("ContextualSupplementary:FractionalWarning",wid);end;p=predict(m,X);
else
    mu=mean(X,1);sd=std(X,0,1);sd(sd==0|~isfinite(sd))=1;Z=(X-mu)./sd;state=2*ones(height(D),1);state(D.outage0==0)=1;state(D.outage0==1)=3;
    B=mnrfit(Z,state,'model','nominal');P=mnrval(B,Z,'model','nominal');part=D.outage0>0&D.outage0<1;theta=fitBeta(Z(part,:),D.outage0(part));mub=logistic([ones(height(D),1),Z]*theta(1:end-1));p=P(:,3)+P(:,2).*mub;
end
end
function X=baseDesign(D,s,levels)
if s.pga_transform=="ln1p",pga=log1p(D.(s.pga_field)/0.1);else,pga=D.(s.pga_field)/0.1;end
if s.context_kind=="density",ctx=log2(D.popDensity_km2);else,ctx=log2(D.population);end
X=[pga,(D.year-2000)/10,ctx,double(D.sys=="W"),double(D.sys=="G")];
for j=2:numel(levels),X(:,end+1)=double(string(D.development_level)==levels(j));end
end

function fit=fitTime(D,s,includeD0,levels)
[X,names,knots]=timeDesign(D,s,includeD0,levels,[]);
y=log(D.(s.outcome));
if s.distribution=="lognormal"
    beta=linearFit(X,y);pred=[ones(height(D),1),X]*beta;
else
    theta=fitAft(X,y,s.distribution);beta=theta(1:end-1);pred=aftMedian(X,theta,s.distribution);
end
fit=struct("beta",beta,"prediction",pred,"names",names,"knots",knots,"spec",s);
end
function [X,names,knots]=timeDesign(D,s,includeD0,levels,knots)
X=baseDesign(D,s,levels);names=["pga";"year";"context";"water_vs_electric";"gas_vs_electric"];
for j=2:numel(levels),names(end+1,1)="income_"+matlab.lang.makeValidName(levels(j));end
if ~includeD0,knots=[];return;end
switch s.d0_kind
    case "log2_d0",X(:,end+1)=log2(D.outage0);names(end+1,1)="d0_main";
    case "raw_d0",X(:,end+1)=D.outage0;names(end+1,1)="d0_main";
    case "log2_d0_complete",X(:,end+1)=log2(D.outage0);X(:,end+1)=double(D.outage0==1);names(end+1:end+2,1)=["d0_main";"d0_complete"];
    case "log2p1_d0",X(:,end+1)=log2(1+D.outage0);names(end+1,1)="d0_main";
    case "spline_log2_d0"
        z=log2(D.outage0);if isempty(knots),knots=quantile(z,[.05 .35 .65 .95]);end;B=rcsBasis(z,knots);X=[X,B];for k=1:size(B,2),names(end+1,1)="d0_spline_"+k;end
    otherwise,error("ContextualSupplementary:D0Kind","Unknown D0 representation %s",s.d0_kind);
end
end


function f=d0Factor(D,s,fit,levels)
lo=D;hi=D;lo.outage0(:)=0.10;hi.outage0(:)=0.20;
[Xlo,~,~]=timeDesign(lo,s,true,levels,fit.knots);[Xhi,~,~]=timeDesign(hi,s,true,levels,fit.knots);
f=exp(mean((Xhi-Xlo)*fit.beta(2:end)));
end



function writeDiagnostics(T,cfg,runDir)
sD=spec("primary","Primary","primary","D0","avg_PGA","raw","population","none","lognormal",nan,"fractional");D0=selectD0(T,sD);p=fitD0(D0,sD,incomeLevels(D0,cfg));
[~,~,bin]=histcounts(p,quantile(p,linspace(0,1,11)));bin(bin==0)=1;G=findgroups(bin);cal=table((1:max(G))',splitapply(@numel,p,G),splitapply(@mean,p,G),splitapply(@mean,D0.outage0,G),'VariableNames',{'bin','record_count','mean_fitted_D0','mean_observed_D0'});writetable(cal,fullfile(runDir,"figure_data","s4b_D0_calibration.csv"));
sF=spec("primary","Primary","primary","T90","avg_PGA","raw","population","none","lognormal",nan,"linear");F=selectTime(T,sF,false);ff=fitTime(F,sF,false,incomeLevels(F,cfg));writeTimeDiagnostic(F,ff,"full_T90",runDir);
sA=spec("primary","Primary","primary","T90","avg_PGA","raw","population","log2_d0","lognormal",nan,"linear");A=selectTime(T,sA,true);fa=fitTime(A,sA,true,incomeLevels(A,cfg));writeTimeDiagnostic(A,fa,"aligned_plus_D0",runDir);
end
function writeTimeDiagnostic(D,fit,label,runDir)
y=log(D.(fit.spec.outcome));res=y-fit.prediction;n=numel(res);q=norminv(((1:n)'-.5)/n);
writetable(table(D.record_id,D.eid,fit.prediction,res,'VariableNames',{'record_id','eid','fitted_log_time','residual'}),fullfile(runDir,"figure_data","s4_"+label+"_residuals.csv"));
writetable(table(q,sort(res),'VariableNames',{'normal_quantile','residual_quantile'}),fullfile(runDir,"figure_data","s4_"+label+"_qq.csv"));
end

function writeBootstrapStatus(boot,done,n,runDir)
names=fieldnames(boot);rows=cell(0,4);
for k=1:numel(names),B=boot.(names{k});for j=1:size(B,2),rows(end+1,:)={string(names{k}),j,sum(isfinite(B(:,j))),n};end,end
R=cell2table(rows,'VariableNames',{'quantity','specification_index','successful_replicates','requested_replicates'});writetable(R,fullfile(runDir,"validation","bootstrap_success_by_specification.csv"));
writeJson(fullfile(runDir,"manifest","validation.json"),struct("status","complete","completed_replicates",sum(done),"requested_replicates",n,"minimum_successful_specification",min(R.successful_replicates)));
end
function validatePrimary(d0,full,pair,primaryPath,incrementPath)
P=readtable(primaryPath,"TextType","string");I=readtable(incrementPath);
assert(abs(d0.estimate(d0.specification=="primary")-P.fitted_r2(P.model=="system_and_context_full_D0_record_set"))<1e-10);
assert(abs(full.estimate(full.specification=="primary")-P.fitted_r2(P.model=="system_and_context_full_T90_record_set"))<1e-10);
assert(abs(pair.estimate(pair.specification=="primary")-I.delta_fitted_r2)<1e-10);
end

function levels=incomeLevels(D,cfg)
levels=sort(unique(string(D.development_level)));ref=string(cfg.reference_income_group);if ~ismember(ref,levels),error("ContextualSupplementary:IncomeReference","Reference income group is absent.");end;levels=[ref;levels(levels~=ref)];
end
function m=metric(y,p)
r=y-p;m=struct("sse",sum(r.^2),"sst",sum((y-mean(y)).^2),"r2",1-sum(r.^2)/sum((y-mean(y)).^2));
end
function idx=sampleEventRows(eid)
events=unique(eid,"stable");draw=events(randi(numel(events),numel(events),1));idx=zeros(0,1);for k=1:numel(draw),idx=[idx;find(eid==draw(k))];end
end
function beta=linearFit(X,y)
A=[ones(size(X,1),1),X];if rank(A)<size(A,2),beta=lsqminnorm(A,y);else,beta=A\y;end
end
function B=rcsBasis(x,k)
k=unique(k);if numel(k)<4,error("ContextualSupplementary:Spline","Insufficient spline knots.");end;last=k(end);penult=k(end-1);B=x;for j=1:numel(k)-2,h=max(x-k(j),0).^3-((last-k(j))/(last-penult))*max(x-penult,0).^3+((penult-k(j))/(last-penult))*max(x-last,0).^3;B(:,end+1)=h/(last-k(1))^2;end
end
function theta=fitBeta(X,y)
A=[ones(size(X,1),1),X];init=zeros(size(A,2)+1,1);init(1)=log(mean(y)/(1-mean(y)));init(end)=log(10);o=optimoptions('fminunc','Display','off','Algorithm','quasi-newton','MaxIterations',1000);[theta,~,flag]=fminunc(@(v)betaNll(v,A,y),init,o);if flag<=0,error("ContextualSupplementary:Beta","Beta model did not converge.");end
end
function v=betaNll(theta,A,y)
mu=logistic(A*theta(1:end-1));phi=exp(theta(end));a=max(mu*phi,eps);b=max((1-mu)*phi,eps);v=-sum(gammaln(phi)-gammaln(a)-gammaln(b)+(a-1).*log(y)+(b-1).*log1p(-y));if ~isfinite(v),v=realmax;end
end
function y=logistic(x),y=1./(1+exp(-max(min(x,35),-35)));end
function theta=fitAft(X,y,dist)
A=[ones(size(X,1),1),X];b=linearFit(X,y);init=[b;log(max(std(y-A*b),.05))];o=optimoptions('fminunc','Display','off','Algorithm','quasi-newton','MaxIterations',500);[theta,~,flag]=fminunc(@(v)aftNll(v,A,y,dist),init,o);if flag<=0,error("ContextualSupplementary:AFT","AFT model did not converge.");end
end
function v=aftNll(theta,A,y,dist)
mu=A*theta(1:end-1);s=exp(theta(end));z=(y-mu)/s;if dist=="weibull",lp=-log(s)+z-exp(z);else,lp=-log(s)-z-2*softplus(-z);end;v=-sum(lp);
end
function y=aftMedian(X,t,dist)
mu=[ones(size(X,1),1),X]*t(1:end-1);s=exp(t(end));if dist=="weibull",y=mu+s*log(log(2));else,y=mu;end
end
function y=softplus(x),y=max(x,0)+log1p(exp(-abs(x)));end
function writeJson(path,value)
fid=fopen(path,"w");if fid<0,error("ContextualSupplementary:JsonWrite","Cannot write %s",path);end;c=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(value,"PrettyPrint",true),"char");
end
function hash=fileHash(path)
fid=fopen(path,"r");if fid<0,error("ContextualSupplementary:HashRead","Cannot read %s",path);end
cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,"*uint8")';
digest=java.security.MessageDigest.getInstance("SHA-256");digest.update(bytes);
hash=lower(reshape(dec2hex(typecast(digest.digest(),"uint8"))',1,[]));
end
