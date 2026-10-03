function convergence = compute_sobol_convergence(data,baseCounts)
arguments
    data struct
    baseCounts (1,:) double {mustBeInteger,mustBePositive}
end
if ~isfield(data,'baseID') || numel(data.baseID)~=size(data.A,1)
    error('AllFactorGSA:MissingBaseID', ...
        'data.baseID is required for nested-prefix convergence.');
end
baseCounts=unique(baseCounts(:)');
maxBase=max(data.baseID);
if any(baseCounts>maxBase)
    error('AllFactorGSA:InvalidBaseCount', ...
        'A requested convergence prefix exceeds the available base design.');
end
template=compute_sobol_indices(data);
nTarget=template.nTargets; nOutcome=template.nOutcomes;
first=nan(numel(baseCounts),nTarget,nOutcome);
total=nan(numel(baseCounts),nTarget,nOutcome);
parameterShare=nan(numel(baseCounts),nOutcome);
stochasticShare=nan(numel(baseCounts),nOutcome);
parameterVariance=nan(numel(baseCounts),nOutcome);
valid=false(numel(baseCounts),1);
failureIdentifier=repmat({''},numel(baseCounts),1);
for j=1:numel(baseCounts)
    index=data.baseID(:)<=baseCounts(j);
    subset=data;
    subset.A=data.A(index,:,:);
    subset.B=data.B(index,:,:);
    subset.hybrid=data.hybrid(index,:,:,:);
    subset.baseID=data.baseID(index);
    if isfield(data,'scrambleID'), subset.scrambleID=data.scrambleID(index); end
    if isfield(data,'rowStochasticVariance')
        subset.rowStochasticVariance=data.rowStochasticVariance(index,:);
        subset.stochasticVariance=mean(subset.rowStochasticVariance,1,'omitnan');
    end
    try
        e=compute_sobol_indices(subset);
    catch ME
        if strcmp(ME.identifier,'AllFactorGSA:NonpositiveParameterVariance')
            failureIdentifier{j}=ME.identifier;
            continue;
        end
        rethrow(ME);
    end
    valid(j)=true;
    first(j,:,:)=e.firstOrderRaw;
    total(j,:,:)=e.totalOrderRaw;
    parameterShare(j,:)=e.parameterShare;
    stochasticShare(j,:)=e.stochasticShare;
    parameterVariance(j,:)=e.parameterVariance;
end
convergence.baseCounts=baseCounts;
convergence.firstOrder=first;
convergence.totalOrder=total;
convergence.parameterShare=parameterShare;
convergence.stochasticShare=stochasticShare;
convergence.parameterVariance=parameterVariance;
convergence.valid=valid;
convergence.failureIdentifier=failureIdentifier;
convergence.targetNames=template.targetNames;
convergence.outcomeNames=template.outcomeNames;
end
