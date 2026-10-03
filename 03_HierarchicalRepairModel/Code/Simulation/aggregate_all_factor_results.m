function summary = aggregate_all_factor_results(cfg,design,empirical)

schema=all_factor_result_schema();
nCases=height(design.cases);
nSeeds=cfg.nSeeds;
if mod(nSeeds,2)~=0
    error('AllFactorGSA:OddSeedCount','The seed count must split into two equal batches.');
end
batchSize=nSeeds/2;
if strcmp(cfg.mode,'full') && batchSize~=25
    error('AllFactorGSA:ProductionSeedContract', ...
        'Production requires two independent batches of 25 seeds.');
end

dummy=zeros(1,schema.nColumns);
metricColumns={'tau80_tau50','tau90_tau80','tau95_tau90','kappa','eta90', ...
    'tau80_tau50_step','tau90_tau80_step','tau95_tau90_step', ...
    'kappa_step','eta90_step'};
for j=1:numel(metricColumns), dummy(schema.index.(metricColumns{j}))=1; end
dummy(schema.index.eta90)=0.5; dummy(schema.index.eta90_step)=0.5;
outputTemplate=summarize_all_factor_outputs(dummy,schema.names,empirical,1e-6);
nOutcome=size(outputTemplate.values,2);
sobolOutcomeIndex=find(outputTemplate.outcomeInfo.isTransformed & ...
    ~outputTemplate.outcomeInfo.isMembership);
if numel(sobolOutcomeIndex)~=10
    error('AllFactorGSA:SobolOutcomeContract', ...
        'Sobol analysis requires five primary and five exact-step transforms.');
end

caseBatchSum=zeros(nCases,2,nOutcome);
caseBatchCount=zeros(nCases,2);
caseSumSq=zeros(nCases,nOutcome);
caseCount=zeros(nCases,1);
seen=false(nCases,nSeeds);

nABExpected=2*cfg.nBase*cfg.nSeeds;
abValues=nan(nABExpected,nOutcome);
abRaw=nan(nABExpected,10);
abDiagnosticNames={ ...
    'n_tasks_target','n_tasks_realized', ...
    'crew_fraction_target','crew_fraction_realized', ...
    'duration_depth_ratio','duration_depth_ratio_realized', ...
    'cv_duration_within_layer_target','cv_duration_within_layer_realized', ...
    'cv_weight_target','cv_weight_within_layer_realized', ...
    'rho_dw_target','rho_dw_realized','q_f_target','chi_frontier', ...
    'crew_utilization','maximum_depth','min_layer_count', ...
    'task_layer_entropy_normalized','effective_service_task_count', ...
    'parent_concentration', ...
    'parent_outdegree_cv_realized','parent_use_effective_fraction', ...
    'max_restoration_jump'};
diagnosticColumns=zeros(1,numel(abDiagnosticNames));
for d=1:numel(abDiagnosticNames)
    diagnosticColumns(d)=schema.index.(abDiagnosticNames{d});
end
abDiagnostics=nan(nABExpected,numel(abDiagnosticNames));
nAB=0;
nRowsRead=0;

files=dir(fullfile(cfg.paths.checkpointDir,'*.mat'));
files=files(~contains({files.name},'.partial'));
if isempty(files)
    error('AllFactorGSA:NoCheckpoints', ...
        'No completed MAT checkpoints were found in %s.',cfg.paths.checkpointDir);
end
[~,order]=sort({files.name}); files=files(order);
for f=1:numel(files)
    pathText=fullfile(files(f).folder,files(f).name);
    loaded=load(pathText,'resultMatrix','resultNames','checkpointMeta');
    validate_checkpoint_contract(loaded,schema,cfg,pathText);
    m=loaded.resultMatrix;
    nRowsRead=nRowsRead+size(m,1);
    status=m(:,schema.index.status_code);
    if any(status~=0)
        bad=find(status~=0,1);
        error('AllFactorGSA:FailedTrajectory', ...
            'Checkpoint %s contains failed run %.0f.', ...
            files(f).name,m(bad,schema.index.run_id));
    end

    caseID=round(m(:,schema.index.case_id));
    seedID=round(m(:,schema.index.seed_id));
    batchID=round(m(:,schema.index.seed_batch_id));
    linearRun=sub2ind([nCases,nSeeds],caseID,seedID);
    if any(seen(linearRun)) || numel(unique(linearRun))~=numel(linearRun)
        error('AllFactorGSA:DuplicateTrajectory', ...
            'Checkpoint %s duplicates a case-by-seed trajectory.',files(f).name);
    end
    seen(linearRun)=true;

    prepared=summarize_all_factor_outputs(m,schema.names,empirical,1e-6);
    values=prepared.values;
    linearBatch=sub2ind([nCases,2],caseID,batchID);
    caseBatchCount=caseBatchCount+reshape(accumarray(linearBatch,1, ...
        [nCases*2,1]),nCases,2);
    caseCount=caseCount+accumarray(caseID,1,[nCases,1]);
    for o=1:nOutcome
        contribution=accumarray(linearBatch,values(:,o),[nCases*2,1]);
        caseBatchSum(:,:,o)=caseBatchSum(:,:,o)+reshape(contribution,nCases,2);
        caseSumSq(:,o)=caseSumSq(:,o)+ ...
            accumarray(caseID,values(:,o).^2,[nCases,1]);
    end

    isAB=m(:,schema.index.design_type_code)==1 & ...
        (m(:,schema.index.matrix_id)==1 | m(:,schema.index.matrix_id)==2);
    count=nnz(isAB);
    if count>0
        destination=nAB+(1:count);
        if destination(end)>nABExpected
            error('AllFactorGSA:ABOverflow','More A/B trajectories than expected.');
        end
        abValues(destination,:)=values(isAB,:);
        abRaw(destination,:)=[prepared.primaryRaw(isAB,:),prepared.stepRaw(isAB,:)];
        abDiagnostics(destination,:)=m(isAB,diagnosticColumns);
        nAB=nAB+count;
    end
end

if ~all(seen,'all') || nRowsRead~=cfg.expectedTrajectories
    error('AllFactorGSA:IncompleteCheckpointSet', ...
        'Expected %d unique trajectories, read %d rows and observed %d missing.', ...
        cfg.expectedTrajectories,nRowsRead,nnz(~seen));
end
if any(caseBatchCount(:)~=batchSize) || any(caseCount~=nSeeds)
    error('AllFactorGSA:SeedBatchMismatch', ...
        'Every case must contain two complete equal seed batches.');
end
if nAB~=nABExpected
    error('AllFactorGSA:ABCountMismatch', ...
        'Expected %d A/B trajectories but aggregated %d.',nABExpected,nAB);
end

caseBatchMean=caseBatchSum./caseBatchCount;
caseTotalSum=squeeze(sum(caseBatchSum,2));
caseMean=caseTotalSum./caseCount;
caseVariance=(caseSumSq-(caseTotalSum.^2)./caseCount)./(caseCount-1);
caseVariance=max(caseVariance,0);

[factorInfo,groupInfo]=all_factor_metadata();
[factorData,groupData]=assemble_sobol_data(design, ...
    caseBatchMean(:,:,sobolOutcomeIndex),caseVariance(:,sobolOutcomeIndex), ...
    outputTemplate.outcomeInfo(sobolOutcomeIndex,:),factorInfo,groupInfo);
allowDiagnosticFallback=~strcmp(cfg.mode,'full');
factorData.allowNonpositiveVarianceFallback=allowDiagnosticFallback;
groupData.allowNonpositiveVarianceFallback=allowDiagnosticFallback;
factorEstimate=compute_sobol_indices(factorData);
groupEstimate=compute_sobol_indices(groupData);

factorBlock.estimate=factorEstimate;
groupBlock.estimate=groupEstimate;
if cfg.nBootstrap>0 && numel(unique(factorData.scrambleID))>=2
    factorBlock.bootstrap=bootstrap_sobol_scrambles(factorData, ...
        cfg.nBootstrap,cfg.bootstrapSeed);
    groupBlock.bootstrap=bootstrap_sobol_scrambles(groupData, ...
        cfg.nBootstrap,cfg.bootstrapSeed+1);
end
baseCounts=nested_base_counts(cfg.nBasePerScramble);
factorBlock.convergence=compute_sobol_convergence(factorData,baseCounts);
groupBlock.convergence=compute_sobol_convergence(groupData,baseCounts);

consistency=global_consistency_summary(abValues,abRaw, ...
    outputTemplate.index,outputTemplate.metricNames,empirical);
responseSummary.taskService=summarize_matched_response( ...
    design.taskServiceCases,caseMean,outputTemplate.outcomeInfo, ...
    cfg,'task_service',cfg.nBootstrap,cfg.bootstrapSeed+100);
responseSummary.capacityQ=summarize_matched_response( ...
    design.capacityQCases,caseMean,outputTemplate.outcomeInfo, ...
    cfg,'capacity_q',cfg.nBootstrap,cfg.bootstrapSeed+101);

summary.contractVersion='finite-network-summary-v1';
summary.modelVersion=cfg.modelVersion;
summary.mode=cfg.mode;
summary.experimentID=cfg.experimentID;
summary.empirical=empirical;
summary.createdUTC=char(datetime('now','TimeZone','UTC', ...
    'Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));
summary.factorInfo=factorInfo;
summary.groupInfo=groupInfo;
summary.outcomeInfo=outputTemplate.outcomeInfo;
summary.sobol.factor=factorBlock;
summary.sobol.group=groupBlock;
summary.consistency=consistency;
summary.responseSummary=responseSummary;
summary.diagnostics=struct('mode',cfg.mode,'nCases',nCases, ...
    'nSeeds',nSeeds,'nResponseBackgrounds',cfg.nResponseBackgrounds, ...
    'nResponseScrambles',cfg.nResponseScrambles, ...
    'nBackgroundsPerResponseScramble',cfg.nBackgroundsPerResponseScramble, ...
    'nResponseCells',cfg.nResponseCells, ...
    'seedBatchSize',batchSize,'nTrajectories',nRowsRead, ...
    'nABTrajectories',nAB,'logitEpsilon',1e-6, ...
    'quickVarianceFallbackAllowed',allowDiagnosticFallback, ...
    'factorVarianceFallbackUsed',any(factorEstimate.usedUncorrectedFallback), ...
    'groupVarianceFallbackUsed',any(groupEstimate.usedUncorrectedFallback), ...
    'bootstrapUnit','complete_scramble_block', ...
    'responseBootstrapUnit','complete_response_scramble_block', ...
    'nBootstrap',cfg.nBootstrap, ...
    'responseModulesExcludedFromSobol',true, ...
    'responseModulesExcludedFromGlobalConsistency',true);
summary.diagnostics.realizedParentChoice= ...
    realized_parent_choice_diagnostics(abDiagnostics,abDiagnosticNames);
summary.diagnostics.realizedInput= ...
    realized_input_diagnostics(abDiagnostics,abDiagnosticNames);
end

function validate_checkpoint_contract(x,schema,cfg,pathText)
    required={'resultMatrix','resultNames','checkpointMeta'};
for j=1:numel(required)
    if ~isfield(x,required{j})
        error('AllFactorGSA:CheckpointContract', ...
            '%s omits variable %s.',pathText,required{j});
    end
end
if iscell(x.resultNames)
    checkpointNames=x.resultNames(:);
else
    checkpointNames=cellstr(x.resultNames(:));
end
if ~isequal(checkpointNames,schema.names(:)) || ...
        size(x.resultMatrix,2)~=schema.nColumns
    error('AllFactorGSA:CheckpointSchemaMismatch', ...
        '%s does not use the locked all-factor result schema.',pathText);
end
meta=x.checkpointMeta;
if ~checkpoint_matches_experiment(meta,cfg)
    error('AllFactorGSA:CheckpointProvenance', ...
        '%s belongs to another model or experiment.',pathText);
end
end

function [factorData,groupData]=assemble_sobol_data(design,batchMean, ...
    caseVariance,outcomeInfo,factorInfo,groupInfo)
sobol=design.sobolCases;
[aRows,pairScramble,pairBase]=ordered_role_rows(sobol,'A',0,0);
[bRows,~,~]=ordered_role_rows(sobol,'B',0,0);
nPair=numel(aRows); nOutcome=size(batchMean,3);
A=batchMean(aRows,:,:); B=batchMean(bRows,:,:);
factorHybrid=zeros(nPair,height(factorInfo),2,nOutcome);
for k=1:height(factorInfo)
    rows=ordered_role_rows(sobol,'AB_factor',k,0);
    factorHybrid(:,k,:,:)=reshape(batchMean(rows,:,:),nPair,1,2,nOutcome);
end
groupHybrid=zeros(nPair,height(groupInfo),2,nOutcome);
for k=1:height(groupInfo)
    rows=ordered_role_rows(sobol,'AB_group',0,k);
    groupHybrid(:,k,:,:)=reshape(batchMean(rows,:,:),nPair,1,2,nOutcome);
end
rowStochasticVariance=0.5*(caseVariance(aRows,:)+caseVariance(bRows,:));
stochasticVariance=mean(rowStochasticVariance,1,'omitnan');
factorData=struct('A',A,'B',B,'hybrid',factorHybrid, ...
    'targetNames',{factorInfo.factorName'}, ...
    'outcomeNames',{outcomeInfo.outcomeName'}, ...
    'scrambleID',pairScramble,'baseID',pairBase, ...
    'rowStochasticVariance',rowStochasticVariance, ...
    'stochasticVariance',stochasticVariance);
groupData=factorData;
groupData.hybrid=groupHybrid;
groupData.targetNames=groupInfo.groupName';
end

function [caseRows,scramble,base]=ordered_role_rows(sobol,role,factorID,groupID)
mask=strcmp(sobol.matrixRole,role);
if factorID>0, mask=mask & sobol.hybridFactorID==factorID; end
if groupID>0, mask=mask & sobol.hybridGroupID==groupID; end
selected=sobol(mask,:);
selected=sortrows(selected,{'scrambleID','baseIDWithinScramble'});
caseRows=selected.caseID;
scramble=selected.scrambleID;
base=selected.baseIDWithinScramble;
if size(unique([scramble,base],'rows'),1)~=height(selected)
    error('AllFactorGSA:DuplicateSobolMatrixRow', ...
        'Role %s does not contain one row per scramble-by-base pair.',role);
end
end

function counts=nested_base_counts(maxBase)
candidate=2.^(6:floor(log2(maxBase)));
candidate=candidate(candidate<=maxBase);
if isempty(candidate) || candidate(end)~=maxBase, candidate=[candidate,maxBase]; end
counts=unique(candidate);
end

function consistency=global_consistency_summary(values,raw,index,metricNames,empirical)
n=size(values,1);
readout=[repmat({'primary'},5,1);repmat({'exact_step'},5,1)];
metricName=[metricNames(:);metricNames(:)];
belowFraction=nan(10,1); marginalMembership=nan(10,1); aboveFraction=nan(10,1);
for j=1:5
    for step=0:1
        r=j+5*step; x=raw(:,r);
        belowFraction(r)=mean(x<empirical.lower(j),'omitnan');
        aboveFraction(r)=mean(x>empirical.upper(j),'omitnan');
        if step
            member=index.(sprintf('step_membership_%s',metricNames{j}));
        else
            member=index.(sprintf('membership_%s',metricNames{j}));
        end
        marginalMembership(r)=mean(values(:,member),'omitnan');
    end
end
nRuns=repmat(n,10,1);
withinFraction=marginalMembership;
if any(abs(belowFraction+withinFraction+aboveFraction-1)>1e-12)
    error('AllFactorGSA:ConsistencyPartition', ...
        'Below/within/above empirical-range fractions must sum to one.');
end
globalTable=table(readout,metricName,nRuns,belowFraction,withinFraction, ...
    aboveFraction,marginalMembership);
allFiveIntersection=table({'primary';'exact_step'},n*[1;1], ...
    [mean(values(:,index.all_five_marginal_interval_intersection),'omitnan'); ...
     mean(values(:,index.step_all_five_marginal_interval_intersection),'omitnan')], ...
    'VariableNames',{'readout','nRuns','allFiveIntersection'});

rows=20; readoutD=cell(rows,1); metricD=cell(rows,1); scale=cell(rows,1);
meanValue=nan(rows,1); p05=meanValue; medianValue=meanValue; p95=meanValue;
row=0;
for step=0:1
    for transformed=0:1
        for j=1:5
            row=row+1;
            if step, readoutD{row}='exact_step'; else, readoutD{row}='primary'; end
            metricD{row}=metricNames{j};
            if transformed
                scale{row}='transformed';
                if step
                    column=index.([transformed_name(metricNames{j}),'_step']);
                else
                    column=index.(transformed_name(metricNames{j}));
                end
                x=values(:,column);
            else
                scale{row}='raw';
                x=raw(:,j+5*step);
            end
            meanValue(row)=mean(x,'omitnan');
            p05(row)=quantile_local(x,0.05);
            medianValue(row)=quantile_local(x,0.50);
            p95(row)=quantile_local(x,0.95);
        end
    end
end
distribution=table(readoutD,metricD,scale,repmat(n,rows,1), ...
    meanValue,p05,medianValue,p95,'VariableNames', ...
    {'readout','metricName','scale','nRuns','mean','p05','median','p95'});
transformedDistribution=distribution(strcmp(distribution.scale,'transformed'),:);
rawDistribution=distribution(strcmp(distribution.scale,'raw'),:);

empiricalP05=repmat(empirical.lower(:),2,1);
empiricalMedian=repmat(empirical.median(:),2,1);
empiricalP95=repmat(empirical.upper(:),2,1);
simulationP05=nan(10,1); simulationMedian=nan(10,1); simulationP95=nan(10,1);
comparisonRole=[repmat({'primary readout'},5,1); ...
    repmat({'exact-step readout sensitivity; primary empirical bounds unchanged'},5,1)];
for r=1:10
    hit=strcmp(rawDistribution.readout,readout{r}) & ...
        strcmp(rawDistribution.metricName,metricName{r});
    simulationP05(r)=rawDistribution.p05(hit);
    simulationMedian(r)=rawDistribution.median(hit);
    simulationP95(r)=rawDistribution.p95(hit);
end
modelDataAudit=table(readout,metricName,nRuns,empiricalP05, ...
    empiricalMedian,empiricalP95,simulationP05,simulationMedian, ...
    simulationP95,belowFraction,withinFraction,aboveFraction,comparisonRole);

disagreement=nan(5,1);
for j=1:5
    primaryColumn=index.(sprintf('membership_%s',metricNames{j}));
    stepColumn=index.(sprintf('step_membership_%s',metricNames{j}));
    disagreement(j)=mean(values(:,primaryColumn)~=values(:,stepColumn),'omitnan');
end
readoutSensitivity=table(metricNames(:),n*ones(5,1),disagreement, ...
    'VariableNames',{'metricName','nRuns','classificationDisagreement'});
consistency.global=globalTable;
consistency.allFiveIntersection=allFiveIntersection;
consistency.distribution=rawDistribution;
consistency.transformedDistribution=transformedDistribution;
consistency.readoutSensitivity=readoutSensitivity;
consistency.modelDataAudit=modelDataAudit;
consistency.empiricalReference=table(metricNames(:),empirical.lower(:), ...
    empirical.median(:),empirical.upper(:), ...
    repmat({char(empirical.intervalLabel)},5,1), ...
    'VariableNames',{'metricName','empiricalP05','empiricalMedian', ...
    'empiricalP95','intervalLabel'});
consistency.meanMarginalConsistency=table({'primary';'exact_step'}, ...
    [mean(values(:,index.mean_marginal_consistency),'omitnan'); ...
     mean(values(:,index.step_mean_marginal_consistency),'omitnan')], ...
    'VariableNames',{'readout','meanMarginalConsistency'});
end

function name=transformed_name(metric)
if strcmp(metric,'eta90'), name='logit_eta90'; else, name=['log_',metric]; end
end

function diagnostics=realized_parent_choice_diagnostics(x,names)
if size(x,2)~=numel(names) || any(~isfinite(x),'all')
    error('AllFactorGSA:RealizedDiagnosticContract', ...
        'Realized parent-choice diagnostics are incomplete or nonfinite.');
end
nMetric=numel(names); nRuns=size(x,1);
metricName=names(:); meanValue=mean(x,1)'; p05=nan(nMetric,1);
medianValue=p05; p95=p05;
for j=1:nMetric
    p05(j)=quantile_local(x(:,j),0.05);
    medianValue(j)=quantile_local(x(:,j),0.50);
    p95(j)=quantile_local(x(:,j),0.95);
end
diagnostics.distribution=table(metricName,repmat(nRuns,nMetric,1), ...
    meanValue,p05,medianValue,p95,'VariableNames', ...
    {'metricName','nRuns','mean','p05','median','p95'});

alpha=x(:,strcmp(names,'parent_concentration'));
nTask=x(:,strcmp(names,'n_tasks_target'));
outdegree=x(:,strcmp(names,'parent_outdegree_cv_realized'));
effective=x(:,strcmp(names,'parent_use_effective_fraction'));
predictor={'log(parent_concentration)';'log(n_tasks_target)'; ...
    'log(parent_concentration)';'log(n_tasks_target)'};
realizedMetric={'parent_outdegree_cv_realized';'parent_outdegree_cv_realized'; ...
    'parent_use_effective_fraction';'parent_use_effective_fraction'};
u={log(alpha),log(nTask),log(alpha),log(nTask)};
v={outdegree,outdegree,effective,effective};
pearsonCorrelation=nan(4,1);
for j=1:4
    c=corrcoef(u{j},v{j}); pearsonCorrelation(j)=c(1,2);
end
diagnostics.association=table(predictor,realizedMetric, ...
    repmat(nRuns,4,1),pearsonCorrelation);
diagnostics.note=[ ...
    'Descriptive design-by-seed diagnostics only; interactions of alpha_p ' ...
    'with task count and layer width are retained in Sobol total-order effects.'];
end

function tab=realized_input_diagnostics(x,names)
diagnosticName={ ...
    'repair_task_count';'relative_repair_capacity';'duration_depth_ratio'; ...
    'within_layer_duration_cv';'within_layer_service_weight_cv'; ...
    'duration_weight_association';'frontier_ready_start_fraction'; ...
    'crew_utilization';'maximum_dependency_depth';'minimum_layer_count'; ...
    'normalized_task_layer_entropy';'effective_service_task_count'; ...
    'parent_choice_concentration'; ...
    'parent_outdegree_cv';'effective_parent_use_fraction'; ...
    'maximum_restoration_jump'};
targetColumn={ ...
    'n_tasks_target';'crew_fraction_target';'duration_depth_ratio'; ...
    'cv_duration_within_layer_target';'cv_weight_target';'rho_dw_target'; ...
    'q_f_target';'';'maximum_depth';'';'';'';'parent_concentration';'';'';''};
realizedColumn={ ...
    'n_tasks_realized';'crew_fraction_realized';'duration_depth_ratio_realized'; ...
    'cv_duration_within_layer_realized';'cv_weight_within_layer_realized'; ...
    'rho_dw_realized';'chi_frontier';'crew_utilization';''; ...
    'min_layer_count';'task_layer_entropy_normalized'; ...
    'effective_service_task_count';''; ...
    'parent_outdegree_cv_realized';'parent_use_effective_fraction'; ...
    'max_restoration_jump'};
role={ ...
    'target-realized pair';'target-realized pair';'target-realized pair'; ...
    'target-realized pair';'target-realized pair'; ...
    'target-realized association proxy'; ...
    'target-realized mediator pair';'realized mediator';'design input'; ...
    'realized diagnostic';'realized diagnostic';'realized diagnostic'; ...
    'design input'; ...
    'realized mediator';'realized mediator';'realized diagnostic'};
nRow=numel(diagnosticName); n=zeros(nRow,1);
targetMin=nan(nRow,1); targetP05=targetMin; targetMedian=targetMin;
targetP95=targetMin; targetMax=targetMin;
realizedMin=targetMin; realizedP05=targetMin; realizedMedian=targetMin;
realizedP95=targetMin; realizedMax=targetMin;
targetRealizedCorrelation=targetMin; meanSignedError=targetMin;
medianAbsoluteError=targetMin; rootMeanSquaredError=targetMin;
for j=1:nRow
    target=[]; realized=[];
    if ~isempty(targetColumn{j}), target=x(:,column_index(names,targetColumn{j})); end
    if ~isempty(realizedColumn{j}), realized=x(:,column_index(names,realizedColumn{j})); end
    if ~isempty(target)
        [targetMin(j),targetP05(j),targetMedian(j),targetP95(j),targetMax(j)]= ...
            five_number(target);
    end
    if ~isempty(realized)
        [realizedMin(j),realizedP05(j),realizedMedian(j),realizedP95(j), ...
            realizedMax(j)]=five_number(realized);
    end
    if ~isempty(target) && ~isempty(realized)
        ok=isfinite(target)&isfinite(realized); n(j)=nnz(ok);
        c=corrcoef(target(ok),realized(ok)); targetRealizedCorrelation(j)=c(1,2);
        if strcmp(role{j},'target-realized pair')
            errorValue=realized(ok)-target(ok);
            meanSignedError(j)=mean(errorValue);
            medianAbsoluteError(j)=quantile_local(abs(errorValue),0.5);
            rootMeanSquaredError(j)=sqrt(mean(errorValue.^2));
        end
    elseif ~isempty(target)
        n(j)=nnz(isfinite(target));
    else
        n(j)=nnz(isfinite(realized));
    end
end
tab=table(diagnosticName,targetColumn,realizedColumn,role,n, ...
    targetMin,targetP05,targetMedian,targetP95,targetMax, ...
    realizedMin,realizedP05,realizedMedian,realizedP95,realizedMax, ...
    targetRealizedCorrelation,meanSignedError,medianAbsoluteError, ...
    rootMeanSquaredError);
end

function index=column_index(names,wanted)
index=find(strcmp(names,wanted),1);
if isempty(index)
    error('AllFactorGSA:MissingRealizedDiagnostic','Missing diagnostic %s.',wanted);
end
end

function [minimum,p05,medianValue,p95,maximum]=five_number(x)
x=x(isfinite(x));
if isempty(x)
    minimum=nan; p05=nan; medianValue=nan; p95=nan; maximum=nan;
    return;
end
minimum=min(x); p05=quantile_local(x,0.05);
medianValue=quantile_local(x,0.50); p95=quantile_local(x,0.95);
maximum=max(x);
end
