function outputs = summarize_all_factor_outputs(resultMatrix,resultNames,empirical,logitEpsilon)

metricNames={'tau80_tau50','tau90_tau80','tau95_tau90','kappa','eta90'};
metricLabels={'log(\tau_{80}/\tau_{50})', ...
    'log(\tau_{90}/\tau_{80})','log(\tau_{95}/\tau_{90})', ...
    'log(\kappa)','logit(\eta_{90})'};
raw=extract_columns(resultMatrix,resultNames,metricNames);
stepNames=strcat(metricNames,'_step');
stepRaw=extract_columns(resultMatrix,resultNames,stepNames);
validate_metrics(raw,'primary');
validate_metrics(stepRaw,'exact-step');

primaryTransformed=transform_metrics(raw,logitEpsilon);
stepTransformed=transform_metrics(stepRaw,logitEpsilon);
lower=empirical.lower;
upper=empirical.upper;
primaryMembership=raw>=lower & raw<=upper;
stepMembership=stepRaw>=lower & stepRaw<=upper;
primaryMeanMarginal=mean(primaryMembership,2);
stepMeanMarginal=mean(stepMembership,2);
primaryAllFive=all(primaryMembership,2);
stepAllFive=all(stepMembership,2);

values=[primaryTransformed,stepTransformed, ...
    primaryMembership,primaryMeanMarginal,primaryAllFive, ...
    stepMembership,stepMeanMarginal,stepAllFive];
primaryOutcomeNames={'log_tau80_tau50','log_tau90_tau80', ...
    'log_tau95_tau90','log_kappa','logit_eta90'};
stepOutcomeNames=strcat(primaryOutcomeNames,'_step');
outcomeName=[primaryOutcomeNames,stepOutcomeNames, ...
    strcat('membership_',metricNames), ...
    {'mean_marginal_consistency'}, ...
    {'all_five_marginal_interval_intersection'}, ...
    strcat('step_membership_',metricNames), ...
    {'step_mean_marginal_consistency'}, ...
    {'step_all_five_marginal_interval_intersection'}]';
outcomeLabel=[metricLabels,strcat(metricLabels,' [exact step]'), ...
    strcat('Within empirical range: ',metricNames), ...
    {'Mean marginal empirical-range consistency'}, ...
    {'All-five marginal-interval intersection'}, ...
    strcat('Exact-step within empirical range: ',metricNames), ...
    {'Exact-step mean marginal consistency'}, ...
    {'Exact-step all-five marginal-interval intersection'}]';
metricColumn=[metricNames,metricNames,metricNames,{'mean_five'}, ...
    {'all_five'},metricNames,{'mean_five'},{'all_five'}]';
readout=[repmat({'primary'},1,5),repmat({'exact_step'},1,5), ...
    repmat({'primary'},1,7),repmat({'exact_step'},1,7)]';
isTransformed=[true(1,10),false(1,14)]';
isMembership=[false(1,10),true(1,14)]';
isSecondary=[false(1,5),true(1,5),false(1,6),true(1,8)]';
outcomeID=(1:numel(outcomeName))';
outcomeInfo=table(outcomeID,outcomeName,outcomeLabel,metricColumn,readout, ...
    isTransformed,isMembership,isSecondary, ...
    'VariableNames',{'outcomeID','outcomeName','outcomeLabel','metricName', ...
    'readout','isTransformed','isMembership','isSecondary'});
index=struct();
for j=1:numel(outcomeName)
    index.(outcomeName{j})=j;
end

outputs.values=values;
outputs.outcomeInfo=outcomeInfo;
outputs.index=index;
outputs.metricNames=metricNames;
outputs.primaryRaw=raw;
outputs.stepRaw=stepRaw;
end

function x=extract_columns(matrix,names,wanted)
if ~iscell(names), names=cellstr(names); end
names=names(:);
index=zeros(1,numel(wanted));
for j=1:numel(wanted)
    q=find(strcmp(names,wanted{j}),1);
    if isempty(q)
        error('AllFactorGSA:MissingResultColumn', ...
            'Required result column %s is absent.',wanted{j});
    end
    index(j)=q;
end
x=double(matrix(:,index));
end

function validate_metrics(x,label)
if any(~isfinite(x),'all') || any(x(:,1:4)<=0,'all') || ...
        any(x(:,5)<0 | x(:,5)>1)
    error('AllFactorGSA:InvalidIndicators', ...
        'The %s indicators are outside their admissible domains.',label);
end
end

function y=transform_metrics(x,epsilon)
y=zeros(size(x));
y(:,1:4)=log(x(:,1:4));
eta=min(max(x(:,5),epsilon),1-epsilon);
y(:,5)=log(eta./(1-eta));
end
