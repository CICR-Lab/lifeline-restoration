function [s13,s14,s15,s16]=compute_analytical_supplementary(cfg)

%% Supplementary Figure S13: repair and reconnection limits
s13=build_repair_reconnection_limits(cfg);

%% Supplementary Figures S14-S15: loss scaling and indicator maps
s15=build_depth_indicator_grid(cfg);
s14=build_depth_loss_fit(s15);

%% Supplementary Figure S16: generalized analytical sensitivities
s16=build_generalized_sensitivities(cfg);
end

function s13=build_repair_reconnection_limits(cfg)
theta=linspace(0,10,2001)'; identity=struct('type','identity');
s13.theta=theta;
s13.responseA=analytical_response('parallel',theta,3,struct( ...
    'distName','exponential','meanDuration',1,'durationCV',1, ...
    'releaseSpec',struct('type','power','q',.75,'h',.70,'s',20)));
s13.depthsB=0:3;
s13.closedB=analytical_response('frontier',theta,s13.depthsB, ...
    struct('meanDuration',1,'releaseSpec',identity));
model.L=3; model.lambda=ones(1,4); model.accessType='single';
model.m=ones(1,4); model.pi=ones(1,4)/4;
model.release=struct('type','identity','q',1,'h',.70,'s',20);
s13.numericalB=solve_generalized_meanfield(model,cfg.solver);
s13.Xnumerical=interp1(s13.numericalB.theta,s13.numericalB.X,theta,'pchip');
closedFull=analytical_response('frontier',s13.numericalB.theta,s13.depthsB, ...
    struct('meanDuration',1,'releaseSpec',identity));
s13.maxError=max(abs(closedFull.taskFraction-s13.numericalB.X),[],'all');
s13.depthsC=1:4;
s13.responseC=analytical_response('parallel',theta,s13.depthsC,struct( ...
    'distName','exponential','meanDuration',1,'durationCV',1, ...
    'releaseSpec',identity));
s13.depthsD=0:3;
s13.frontierRate=analytical_layer_rate(theta,max(s13.depthsD));
s13.baseCompletion=duration_cdf(theta,'exponential',1,1);
end

function grid=build_depth_indicator_grid(cfg)
Lvals=cfg.sensitivity.depthValues(:)';
dvals=cfg.sensitivity.depthNormGrid(:)';
nL=numel(Lvals); nD=numel(dvals); nM=numel(cfg.sensitivity.indicatorNames);
limitNames=cfg.sensitivity.repairLimits;
nR=numel(limitNames);
V=nan(nL,nD,nM,nR);
T50=nan(nL,nD,nR); LossArea=T50; EffectiveDepth=T50;
Pi=cell(nL,nD);
theta=cfg.theta.grid;
for r=1:nR
    limitName=limitNames{r};
    model=model_for_limit(cfg,limitName);
    for i=1:nL
        L=Lvals(i);
        response=analytical_response(limitName,theta,0:L,model);
        layerRestoration=response.restoredFraction;
        for j=1:nD
            pi=service_distribution_for_mean_depth(L,dvals(j));
            met=closed_form_profile_metrics(theta,layerRestoration,pi,limitName);
            V(i,j,:,r)=[met.tau80_tau50 met.tau90_tau80 met.tau95_tau90 ...
                met.kappaInfinite met.eta90Infinite];
            T50(i,j,r)=met.tau50;
            LossArea(i,j,r)=met.areaInfinite;
            EffectiveDepth(i,j,r)=met.effectiveDepth;
            Pi{i,j}=pi;
        end
    end
end
[T,complete]=load_empirical_indicator_ranges( ...
    cfg.data.empiricalRanges,cfg.sensitivity.indicatorNames);
C=nan(size(V));
for r=1:nR
    for m=1:nM
        range=get_empirical_range(T,cfg.sensitivity.indicatorNames{m});
        C(:,:,m,r)=classify_empirical_range(V(:,:,m,r),range);
    end
end
grid.L=Lvals; grid.depthNorm=dvals;
grid.limitNames=limitNames;
grid.limitLabels=cfg.sensitivity.repairLimitLabels;
grid.valuesByLimit=V; grid.classesByLimit=C;
grid.t50ByLimit=T50; grid.lossAreaByLimit=LossArea;
grid.effectiveDepthByLimit=EffectiveDepth; grid.pi=Pi;
grid.values=V(:,:,:,1); grid.classes=C(:,:,:,1);
grid.empiricalRanges=T; grid.empiricalComplete=complete;
end

function model=model_for_limit(cfg,limitName)
model=struct('meanDuration',1,'releaseSpec',struct('type','identity'));
if strcmpi(limitName,'parallel')
    if ~strcmpi(cfg.sensitivity.parallel.durationFamily,'exponential') || ...
            cfg.sensitivity.parallel.meanDuration~=1
        error(['The exact parallel S5 loss indicators currently require ' ...
            'unit-mean exponential task durations.']);
    end
    model.distName=cfg.sensitivity.parallel.durationFamily;
    model.meanDuration=cfg.sensitivity.parallel.meanDuration;
    model.durationCV=cfg.sensitivity.parallel.durationCV;
end
end

function s16=build_generalized_sensitivities(cfg)
families={'dependency','rates','release_duration'};
s16=cell(3,2);
for i=1:3
    s16{i,1}=build_sensitivity_family(cfg,families{i},'frontier');
    s16{i,2}=build_sensitivity_family(cfg,families{i},'parallel');
end
end

function result=build_depth_loss_fit(grid)
[iL,iD,iR]=ndgrid(1:numel(grid.L),1:numel(grid.depthNorm),1:numel(grid.limitNames));
n=numel(iL); RepairLimit=strings(n,1); MaximumDepth=zeros(n,1);
NormalizedServiceWeightedDepth=zeros(n,1); EffectiveDepth=zeros(n,1);
Tau50=zeros(n,1); LossRest=zeros(n,1); Kappa=zeros(n,1); Pi=strings(n,1);
for row=1:n
    a=iL(row); b=iD(row); c=iR(row);
    RepairLimit(row)=string(grid.limitNames{c}); MaximumDepth(row)=grid.L(a);
    NormalizedServiceWeightedDepth(row)=grid.depthNorm(b);
    EffectiveDepth(row)=grid.effectiveDepthByLimit(a,b,c);
    Tau50(row)=grid.t50ByLimit(a,b,c); LossRest(row)=grid.lossAreaByLimit(a,b,c);
    Kappa(row)=grid.valuesByLimit(a,b,4,c);
    Pi(row)="["+strjoin(compose('%.8g',grid.pi{a,b}),',')+"]";
end
Dmax=ones(n,1); DmaxTau50=Tau50;
caseTable=table(RepairLimit,MaximumDepth,NormalizedServiceWeightedDepth, ...
    EffectiveDepth,Dmax,Tau50,DmaxTau50,LossRest,Kappa,Pi);
scope=["Pooled";string(grid.limitLabels(:))];
K=zeros(numel(scope),1); RSquared=K; RMSE=K; MAE=K; N=K;
masks=cell(numel(scope),1); masks{1}=true(n,1);
for r=1:numel(grid.limitNames), masks{r+1}=RepairLimit==string(grid.limitNames{r}); end
for r=1:numel(scope)
    f=fit_origin(DmaxTau50(masks{r}),LossRest(masks{r}));
    K(r)=f.k; RSquared(r)=f.rSquared; RMSE(r)=f.rmse; MAE(r)=f.mae; N(r)=f.n;
end
fitSummary=table(scope,K,RSquared,RMSE,MAE,N,'VariableNames', ...
    {'Scope','K','RSquared','RMSE','MAE','N'});
caseTable.PooledPredictedLoss=K(1)*DmaxTau50;
caseTable.PooledResidual=LossRest-caseTable.PooledPredictedLoss;
result=struct('repairLimits',{grid.limitNames},'limitLabels',{grid.limitLabels}, ...
    'caseTable',caseTable,'fitSummary',fitSummary);
end

function fit=fit_origin(x,y)
k=sum(x.*y)/sum(x.^2); residual=y-k*x;
fit=struct('k',k,'rSquared',1-sum(residual.^2)/sum((y-mean(y)).^2), ...
    'rmse',sqrt(mean(residual.^2)),'mae',mean(abs(residual)),'n',numel(x));
end
