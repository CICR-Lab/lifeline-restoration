function mapped = map_all_factor_inputs(unitPoints,cfg)

u=double(unitPoints);
if size(u,2)~=cfg.nFactors || any(~isfinite(u(:))) || ...
        any(u(:)<0 | u(:)>1)
    error('RepairTask:InvalidAllFactorUnitPoints', ...
        'unitPoints must be a finite n-by-%d matrix on [0,1].', ...
        cfg.nFactors);
end
u=min(max(u,eps('double')),1-eps('double'));
n=size(u,1);

nTasks=round(log_map(u(:,1),cfg.nTasksBounds));
maximumDepth=discrete_numeric(u(:,2),cfg.maximumDepth);
[taskProfileCode,taskLayerProfile]=categorical_map( ...
    u(:,3),cfg.taskLayerProfiles);
taskCountConcentration=log_map( ...
    u(:,4),cfg.taskCountConcentrationBounds);
[dependencyCode,dependencyLogic]=categorical_map( ...
    u(:,5),cfg.dependencyLogic);
parentConcentration=log_map(u(:,6),cfg.parentConcentrationBounds);
[serviceProfileCode,serviceLayerProfile]=categorical_map( ...
    u(:,7),cfg.serviceLayerProfiles);
durationDepthRatio=log_map(u(:,8),cfg.durationDepthRatioBounds);
cvDurationWithinLayer=linear_map(u(:,9),cfg.cvDurationBounds);
cvWeightWithinLayer=linear_map(u(:,10),cfg.cvWeightBounds);
rhoDurationWeight=linear_map(u(:,11),cfg.rhoDurationWeightBounds);
crewFraction=log_map(u(:,12),cfg.crewFractionBounds);
qF=linear_map(u(:,13),cfg.qFBounds);

mapped=table(nTasks,maximumDepth,taskProfileCode,taskLayerProfile, ...
    taskCountConcentration,dependencyCode,dependencyLogic, ...
    parentConcentration,serviceProfileCode,serviceLayerProfile, ...
    durationDepthRatio,cvDurationWithinLayer,cvWeightWithinLayer, ...
    rhoDurationWeight,crewFraction,qF);
if height(mapped)~=n
    error('RepairTask:AllFactorMappingSize', ...
        'Mapped input table has an unexpected number of rows.');
end
end

function x=log_map(u,bounds)
x=bounds(1).*(bounds(2)./bounds(1)).^u;
end

function x=linear_map(u,bounds)
x=bounds(1)+(bounds(2)-bounds(1)).*u;
end

function x=discrete_numeric(u,support)
k=numel(support);
index=min(floor(u.*k)+1,k);
x=reshape(support(index(:)),[],1);
end

function [code,label]=categorical_map(u,support)
k=numel(support);
code=min(floor(u.*k)+1,k);
code=code(:);
label=reshape(support(code),[],1);
end
