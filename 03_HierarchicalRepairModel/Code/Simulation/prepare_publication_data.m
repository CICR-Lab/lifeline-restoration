function pub = prepare_publication_data(summary)

primaryReadout='primary';
responseOutcome='mean_marginal_consistency';

factorInfo = sortrows(summary.factorInfo,'factorID');
groupInfo = sortrows(summary.groupInfo,'groupID');
if isfield(summary,'diagnostics'), diagnostics=summary.diagnostics;
else, diagnostics=struct; end
if isfield(summary,'mode') && ~isfield(diagnostics,'mode')
    diagnostics.mode=summary.mode;
end
validate_full_inference(summary,diagnostics);

[outcomeNames,outcomeLabels,metricNames] = primary_outcomes(summary);

pub.variance = make_variance_table(summary,outcomeNames,outcomeLabels);
pub.factorSobol = make_sobol_table(summary.sobol.factor,factorInfo, ...
    outcomeNames,outcomeLabels,'factor');
pub.groupSobol = make_sobol_table(summary.sobol.group,groupInfo, ...
    outcomeNames,outcomeLabels,'group');
pub.factorSobol.readout=repmat("primary",height(pub.factorSobol),1);
pub.groupSobol.readout=repmat("primary",height(pub.groupSobol),1);
pub.consistency = make_consistency_table(summary,primaryReadout,metricNames);
taskServiceAll = normalize_response_all( ...
    summary.responseSummary.taskService.cellTable,'task_service');
capacityQAll = normalize_response_all( ...
    summary.responseSummary.capacityQ.cellTable,'capacity_q');
pub.taskService = select_response(taskServiceAll,responseOutcome,'task_service');
pub.capacityQ = select_response(capacityQAll,responseOutcome,'capacity_q');
pub.modelDataAudit = required_model_data_audit(summary.consistency);
pub.tableS8 = build_table_s8(factorInfo);
end

function validate_full_inference(summary,diagnostics)
if ~isfield(diagnostics,'mode') || ...
        ~strcmpi(string(diagnostics.mode),'full'), return; end
blocks={summary.sobol.factor,summary.sobol.group};
for j=1:2
    if ~isfield(blocks{j},'bootstrap') || ...
            ~isfield(blocks{j}.bootstrap,'nBootstrap') || ...
            double(blocks{j}.bootstrap.nBootstrap)~=2000 || ...
            ~isfield(blocks{j}.bootstrap,'bootstrapUnit') || ...
            ~strcmp(string(blocks{j}.bootstrap.bootstrapUnit),'complete_scramble_block')
        error('AllFactorPublication:ProductionSobolBootstrap', ...
            ['Production factor and group Sobol results require 2,000 ' ...
             'complete-scramble bootstrap replicates.']);
    end
end
end

function tab=required_model_data_audit(consistency)
if ~isfield(consistency,'modelDataAudit') || ...
        ~istable(consistency.modelDataAudit)
    error('AllFactorPublication:MissingModelDataAudit', ...
        ['summary.consistency.modelDataAudit is required for the ' ...
         'distributional-calibration audit.']);
end
tab=consistency.modelDataAudit;
required={'readout','metricName','nRuns','empiricalP05','empiricalMedian', ...
    'empiricalP95','simulationP05','simulationMedian','simulationP95', ...
    'belowFraction','withinFraction','aboveFraction','comparisonRole'};
if ~all(ismember(required,tab.Properties.VariableNames)) || height(tab)~=10
    error('AllFactorPublication:ModelDataAuditSchema', ...
        'modelDataAudit must contain the frozen 10-row distribution-audit contract.');
end
numericFields={'empiricalP05','empiricalMedian','empiricalP95', ...
    'simulationP05','simulationMedian','simulationP95','belowFraction', ...
    'withinFraction','aboveFraction'};
for j=1:numel(numericFields)
    if any(~isfinite(tab.(numericFields{j})))
        error('AllFactorPublication:ModelDataAuditValues', ...
            'modelDataAudit contains nonfinite %s.',numericFields{j});
    end
end
if any(~isfinite(tab.nRuns) | tab.nRuns<1 | mod(tab.nRuns,1)~=0)
    error('AllFactorPublication:ModelDataAuditRuns', ...
        'modelDataAudit nRuns must be positive integers.');
end
if any(tab.empiricalP05>tab.empiricalMedian) || ...
        any(tab.empiricalMedian>tab.empiricalP95) || ...
        any(tab.simulationP05>tab.simulationMedian) || ...
        any(tab.simulationMedian>tab.simulationP95) || ...
        any(tab.belowFraction<0 | tab.withinFraction<0 | tab.aboveFraction<0) || ...
        any(tab.belowFraction>1 | tab.withinFraction>1 | tab.aboveFraction>1) || ...
        any(abs(tab.belowFraction+tab.withinFraction+tab.aboveFraction-1)>1e-8)
    error('AllFactorPublication:ModelDataAuditFractions', ...
        'Below, within and above fractions must sum to one in every audit row.');
end
end

function [names,labels,metrics] = primary_outcomes(summary)
info = summary.outcomeInfo;
variables = string(info.Properties.VariableNames);
mask = true(height(info),1);
if any(variables=="isTransformed"), mask = mask & logical(info.isTransformed); end
if any(variables=="isMembership"), mask = mask & ~logical(info.isMembership); end
if any(variables=="isSecondary"), mask = mask & ~logical(info.isSecondary); end
preferred = ["log_tau80_tau50","log_tau90_tau80", ...
    "log_tau95_tau90","log_kappa","logit_eta90"];
available = string(info.outcomeName(mask));
if all(ismember(preferred,available))
    row = zeros(5,1);
    for j=1:5, row(j)=find(string(info.outcomeName)==preferred(j),1); end
else
    row = find(mask);
    if numel(row)~=5
        error('AllFactorPublication:PrimaryOutcomeCount', ...
            'Exactly five primary transformed outcomes are required.');
    end
end
names = string(info.outcomeName(row));
labels = string(info.outcomeLabel(row));
metrics = string(info.metricName(row));
end

function tab = make_variance_table(summary,outcomeNames,outcomeLabels)
estimate = summary.sobol.factor.estimate;
idx = outcome_index(estimate,outcomeNames);
tab = table(outcomeNames(:),outcomeLabels(:), ...
    estimate.parameterVariance(idx)',estimate.stochasticVariance(idx)', ...
    estimate.parameterShare(idx)',estimate.stochasticShare(idx)', ...
    'VariableNames',{'outcomeName','outcomeLabel','designVariance', ...
    'stochasticVariance','designShare','stochasticShare'});
if isfield(summary.sobol.factor,'bootstrap')
    b=summary.sobol.factor.bootstrap;
    tab.designShareP025 = row_values(b,'parameterShareP025',idx);
    tab.designShareP975 = row_values(b,'parameterShareP975',idx);
    tab.stochasticShareP025 = row_values(b,'stochasticShareP025',idx);
    tab.stochasticShareP975 = row_values(b,'stochasticShareP975',idx);
end
end

function tab = make_sobol_table(block,metadata,outcomeNames,outcomeLabels,kind)
estimate=block.estimate;
idx=outcome_index(estimate,outcomeNames);
if strcmp(kind,'factor')
    idName='factorID'; nameName='factorName'; labelName='factorLabel';
else
    idName='groupID'; nameName='groupName'; labelName='groupLabel';
end
nTarget=height(metadata); nOutcome=numel(idx); n=nTarget*nOutcome;
targetID=zeros(n,1); targetName=strings(n,1); targetLabel=strings(n,1);
outName=strings(n,1); outLabel=strings(n,1);
S1=nan(n,1); ST=nan(n,1); gap=nan(n,1); row=0;
for i=1:nTarget
    for o=1:nOutcome
        row=row+1;
        targetID(row)=metadata.(idName)(i);
        targetName(row)=string(metadata.(nameName)(i));
        targetLabel(row)=string(metadata.(labelName)(i));
        outName(row)=outcomeNames(o); outLabel(row)=outcomeLabels(o);
        S1(row)=estimate.firstOrderRaw(i,idx(o));
        ST(row)=estimate.totalOrderRaw(i,idx(o));
        gap(row)=estimate.interactionGapRaw(i,idx(o));
    end
end
tab=table(targetID,targetName,targetLabel,outName,outLabel,S1,ST,gap, ...
    max(min(S1,1),0),max(min(ST,1),0), ...
    'VariableNames',{idName,nameName,labelName,'outcomeName','outcomeLabel', ...
    'firstOrderRaw','totalOrderRaw','interactionGapRaw', ...
    'firstOrderDisplay','totalOrderDisplay'});
if isfield(block,'bootstrap')
    b=block.bootstrap; row=0;
    p1=nan(n,1); p2=p1; p3=p1; p4=p1;
    for i=1:nTarget
        for o=1:nOutcome
            row=row+1;
            p1(row)=matrix_value(b,'firstOrderP025',i,idx(o));
            p2(row)=matrix_value(b,'firstOrderP975',i,idx(o));
            p3(row)=matrix_value(b,'totalOrderP025',i,idx(o));
            p4(row)=matrix_value(b,'totalOrderP975',i,idx(o));
        end
    end
    tab.firstOrderP025=p1; tab.firstOrderP975=p2;
    tab.totalOrderP025=p3; tab.totalOrderP975=p4;
    if isfield(b,'nBootstrap'), tab.nBootstrap=repmat(double(b.nBootstrap),n,1);
    else, tab.nBootstrap=nan(n,1); end
    if isfield(b,'bootstrapSeed'), tab.bootstrapSeed=repmat(double(b.bootstrapSeed),n,1);
    else, tab.bootstrapSeed=nan(n,1); end
    if isfield(b,'bootstrapUnit'), tab.bootstrapUnit=repmat(string(b.bootstrapUnit),n,1);
    else, tab.bootstrapUnit=repmat("not_reported",n,1); end
else
    tab.firstOrderP025=nan(n,1); tab.firstOrderP975=nan(n,1);
    tab.totalOrderP025=nan(n,1); tab.totalOrderP975=nan(n,1);
    tab.nBootstrap=zeros(n,1); tab.bootstrapSeed=nan(n,1);
    tab.bootstrapUnit=repmat("not_available",n,1);
end
tab.estimator=repmat("direct noise-corrected pick-freeze",n,1);
tab.surrogateUsed=false(n,1);
tab.varianceDenominator=repmat("design-driven conditional-mean variance",n,1);
if strcmp(kind,'factor')
    tab.groupID=zeros(n,1); tab.groupName=strings(n,1); tab.groupLabel=strings(n,1);
    for r=1:n
        i=find(metadata.factorID==tab.factorID(r),1);
        tab.groupID(r)=metadata.groupID(i);
        tab.groupName(r)=string(metadata.groupName(i));
        tab.groupLabel(r)=string(metadata.groupLabel(i));
    end
end
end

function tab = make_consistency_table(summary,primaryReadout,metricNames)
g=summary.consistency.global;
required={'readout','metricName','marginalMembership'};
if ~all(ismember(required,g.Properties.VariableNames))
    error('AllFactorPublication:ConsistencySchema', ...
        'consistency.global lacks required columns.');
end
readouts=string(g.readout);
primary=strcmpi(readouts,primaryReadout);
if ~any(primary)
    primary=~contains(lower(readouts),'step');
end
if ~any(primary), primary=readouts==readouts(1); end
metricLabels=raw_metric_labels(metricNames);
tab=table((1:5)',metricNames(:),metricLabels(:),nan(5,1),nan(5,1), ...
    'VariableNames',{'displayOrder','metricName','metricLabel','estimate','nRuns'});
for j=1:5
    hit=find(primary & strcmpi(string(g.metricName),metricNames(j)),1);
    if isempty(hit)
        error('AllFactorPublication:MissingConsistencyMetric', ...
            'Missing primary marginal consistency for %s.',metricNames(j));
    end
    tab.estimate(j)=g.marginalMembership(hit);
    if ismember('nRuns',g.Properties.VariableNames), tab.nRuns(j)=g.nRuns(hit); end
end
tab.role=repmat("primary marginal interval",5,1);
tab(end+1,:)={6,"mean_five","Mean across five indicators", ...
    mean(tab.estimate),tab.nRuns(1),"mean of primary marginal intervals"};

allFive=extract_all_five(summary,primary,readouts);
tab(end+1,:)={7,"all_five_intersection", ...
    "All-five marginal-interval intersection",allFive,tab.nRuns(1), ...
    "secondary simultaneous intersection"};
end

function labels=raw_metric_labels(metrics)
metrics=string(metrics); labels=replace(metrics,'_','/');
for j=1:numel(metrics)
    switch lower(metrics(j))
        case 'tau80_tau50', labels(j)='tau80/tau50';
        case 'tau90_tau80', labels(j)='tau90/tau80';
        case 'tau95_tau90', labels(j)='tau95/tau90';
        case 'kappa', labels(j)='kappa';
        case 'eta90', labels(j)='eta90';
    end
end
end

function value=extract_all_five(summary,primary,readouts)
value=nan;
if isfield(summary.consistency,'allFiveIntersection')
    a=summary.consistency.allFiveIntersection;
    if istable(a)
        if ismember('readout',a.Properties.VariableNames)
            hit=find(strcmpi(string(a.readout),readouts(find(primary,1))),1);
            if isempty(hit), hit=1; end
        else
            hit=1;
        end
        candidates={'allFiveIntersection','estimate','value'};
        for j=1:numel(candidates)
            if ismember(candidates{j},a.Properties.VariableNames)
                value=a.(candidates{j})(hit); return;
            end
        end
    elseif isstruct(a)
        candidates={'allFiveIntersection','estimate','value'};
        for j=1:numel(candidates)
            if isfield(a,candidates{j})
                v=a.(candidates{j});
                if isnumeric(v) && ~isempty(v), value=v(1); return; end
            end
        end
    elseif isnumeric(a) && isscalar(a)
        value=a;
    end
end
g=summary.consistency.global;
if isnan(value) && ismember('allFiveIntersection',g.Properties.VariableNames)
    hit=find(primary,1); value=g.allFiveIntersection(hit);
end
end

function tab=normalize_response_all(tab,moduleName)
if ~istable(tab)
    error('AllFactorPublication:ResponseSchema','%s cellTable must be a table.',moduleName);
end
required={'responseCellID','responseAxis1ID','responseAxis2ID','outcomeName','nBackgrounds', ...
    'nResponseScrambles','nBackgroundsPerScramble','nSeeds', ...
    'mean','meanP025','meanP975','median','p05','p10','p90','p95', ...
    'nBootstrap','bootstrapUnit'};
required=[required,{'responseAxis1Value','responseAxis2Value', ...
    'responseAxis1Label','responseAxis2Label'}];
if ~all(ismember(required,tab.Properties.VariableNames))
    error('AllFactorPublication:ResponseSchema', ...
        '%s cellTable lacks required columns.',moduleName);
end
if any(~strcmp(string(tab.bootstrapUnit),'complete_response_scramble_block'))
    error('AllFactorPublication:ResponseBootstrapUnit', ...
        '%s cell means must use complete response-scramble-block resampling.',moduleName);
end
tab.module=repmat(string(moduleName),height(tab),1);
end

function tab=select_response(tab,responseOutcome,moduleName)
hit=strcmpi(string(tab.outcomeName),string(responseOutcome));
if ~any(hit)
    alternatives=contains(lower(string(tab.outcomeName)),'mean') & ...
        contains(lower(string(tab.outcomeName)),'consistency');
    hit=alternatives;
end
if ~any(hit)
    error('AllFactorPublication:MissingResponseOutcome', ...
        '%s lacks mean marginal consistency response rows.',moduleName);
end
tab=tab(hit,:);
end

function tab=build_table_s8(f)
keep={'factorID','factorName','factorLabel','groupID','groupName','groupLabel', ...
    'designSupport','samplingDistribution','designRole','factorDefinition'};
keep=keep(ismember(keep,f.Properties.VariableNames));
tab=f(:,keep);
end


function idx=outcome_index(estimate,names)
if ~isfield(estimate,'outcomeNames')
    error('AllFactorPublication:OutcomeNames','Sobol estimate lacks outcomeNames.');
end
idx=names_to_index(string(estimate.outcomeNames),names);
end

function idx=names_to_index(source,wanted)
idx=zeros(1,numel(wanted));
for j=1:numel(wanted)
    hit=find(strcmpi(source,wanted(j)),1);
    if isempty(hit), error('AllFactorPublication:OutcomeMapping', ...
            'Cannot map outcome %s.',wanted(j)); end
    idx(j)=hit;
end
end

function values=row_values(s,field,idx)
if isfield(s,field), x=s.(field); values=reshape(x(idx),[],1);
else, values=nan(numel(idx),1); end
end

function value=matrix_value(s,field,i,j)
if isfield(s,field), x=s.(field); value=x(i,j); else, value=nan; end
end
