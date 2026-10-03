function [net,realized]=af_generate_repair_task_network(caseRow,seedID,cfg)

if ~strcmp(cfg.durationFamily,'gamma') || ...
        ~strcmp(cfg.gammaSampler,'inverse_cdf_gammaincinv') || ...
        ~strcmp(cfg.serviceReleaseRule,'task_weight_on_reconnection') || ...
        ~strcmp(cfg.durationLayerNormalization, ...
        'exact_layer_mean_affine_floor')
    error('RepairTask:UnimplementedModelVariant', ...
        ['The primary code implements inverse-CDF gamma sampling, exact ' ...
         'layer-mean duration normalization and task-weight release only.']);
end

L=caseRow.maximumDepth;
N=round(caseRow.nTasks);
if N < cfg.minTasksPerLayer * (L+1)
    error('RepairTask:InsufficientTasks', ...
        'N=%d cannot support %d layers with the required minimum.', N, L+1);
end

blockID = af_make_random_block_id( ...
    caseRow.randomBlockID,seedID,cfg);
layerStream = make_component_stream(cfg, blockID, 1);
edgeStream = make_component_stream(cfg, blockID, 2);
durationStream = make_component_stream(cfg, blockID, 3);
weightStream = make_component_stream(cfg, blockID, 4);
priorityStream = make_component_stream(cfg, blockID, 5);

taskProfile=caseRow.taskLayerProfile{1};
serviceProfile=caseRow.serviceLayerProfile{1};
dependencyLogic=caseRow.dependencyLogic{1};

taskProb = layer_profile_weights(taskProfile, L);
alpha=caseRow.taskCountConcentration.*taskProb;
randomShares = gamma_from_stream(layerStream, alpha, size(alpha));
randomShares = max(randomShares,realmin);
randomShares = randomShares ./ sum(randomShares);
layerCounts = allocate_integer_counts(N, randomShares, cfg.minTasksPerLayer);

layer = zeros(N,1);
layerTasks = cell(L+1,1);
cursor = 0;
for ell = 0:L
    ids = (cursor+1):(cursor+layerCounts(ell+1));
    layer(ids) = ell;
    layerTasks{ell+1} = ids(:);
    cursor = cursor + layerCounts(ell+1);
end

layerDurationMean=caseRow.durationDepthRatio.^((0:L)'./L);
overallDurationMean = sum(layerCounts .* layerDurationMean) ./ N;
layerDurationMean = layerDurationMean ./ overallDurationMean;
duration = zeros(N,1);
durationLayerCV = nan(L+1,1);
for ell = 0:L
    ids = layerTasks{ell+1};
    rawDuration = gamma_unit_mean(durationStream, ...
        caseRow.cvDurationWithinLayer,numel(ids));
    rawDuration = rawDuration ./ mean(rawDuration);

    relativeMinimum = cfg.minimumTaskDuration ./ ...
        layerDurationMean(ell+1);
    floorFraction = max(100.*eps('double'), 1.01.*relativeMinimum);
    if floorFraction >= 1
        error('RepairTask:InfeasibleDurationFloor', ...
            ['The minimum task duration is incompatible with the ' ...
             'prescribed layer mean.']);
    end
    rawDuration = floorFraction + ...
        (1-floorFraction) .* rawDuration;
    duration(ids) = layerDurationMean(ell+1) .* rawDuration;
    durationLayerCV(ell+1)=coefficient_of_variation(duration(ids));
end

piLayer = layer_profile_weights(serviceProfile, L);
weight = zeros(N,1);
layerCV = nan(L+1,1);
layerRho = nan(L+1,1);
for ell = 0:L
    ids = layerTasks{ell+1};
    rawWeight = gamma_unit_mean(weightStream, ...
        caseRow.cvWeightWithinLayer,numel(ids));
    rawWeight = impose_rank_association(duration(ids), rawWeight, ...
        caseRow.rhoDurationWeight,weightStream);
    rawWeight = rawWeight ./ sum(rawWeight) .* piLayer(ell+1);
    weight(ids) = rawWeight;
    layerCV(ell+1) = coefficient_of_variation(rawWeight);
    layerRho(ell+1) = spearman_pair(duration(ids), rawWeight);
end
weight = weight ./ sum(weight);

parentCount = zeros(N,1,'uint8');
children = cell(N,1);
for ell = 1:L
    parentIDs = layerTasks{ell};
    childIDs = layerTasks{ell+1};
    prop = gamma_from_stream(edgeStream, ...
        caseRow.parentConcentration,[numel(parentIDs),1]);
    prop = max(prop,realmin);
    prop = prop ./ sum(prop);
    switch dependencyLogic
        case 'single'
            inDegree = 1;
        case {'and2','or2'}
            inDegree = 2;
        otherwise
            error('RepairTask:UnknownDependencyLogic', ...
                'Unknown dependency logic: %s', dependencyLogic);
    end
    if numel(parentIDs) < 2
        error('RepairTask:InsufficientParents', ...
            'Layer %d has fewer than two distinct candidate parents.', ell-1);
    end
    for q = 1:numel(childIDs)
        candidateLocal=weighted_without_replacement(edgeStream,prop,2);
        chosen=parentIDs(candidateLocal(1:inDegree));
        child = childIDs(q);
        parentCount(child) = uint8(inDegree);
        for h = 1:inDegree
            p = chosen(h);
            children{p}(end+1,1) = child;
        end
    end
end

net.N = N;
net.layer = layer;
net.parentCount = parentCount;
net.children = children;
net.dependencyLogic = dependencyLogic;
net.duration = duration;
net.weight = weight;
net.randomPriority = rand(priorityStream,N,1);
net.hybridChoice = rand(priorityStream,N,1);

validate_generated_network(net,L,cfg);

realized.N = N;
realized.minLayerCount = min(layerCounts);
realized.effectiveServiceTaskCount = 1 ./ sum(weight.^2);
realized.durationDepthRatio = mean(duration(layer==L)) ./ ...
    mean(duration(layer==0));
outdegree = cellfun(@numel,children);
outdegreeCV = nan(L,1);
parentUseFraction = nan(L,1);
parentUseWeight = zeros(L,1);
for ell=0:(L-1)
    layerOutdegree = outdegree(layer==ell);
    outdegreeCV(ell+1)=coefficient_of_variation(layerOutdegree);
    parentUseWeight(ell+1) = sum(layerOutdegree);
    if parentUseWeight(ell+1) > 0
        parentProbability = layerOutdegree(layerOutdegree>0) ./ ...
            parentUseWeight(ell+1);
        effectiveParentCount = exp(-sum(parentProbability .* ...
            log(parentProbability)));
        parentUseFraction(ell+1) = effectiveParentCount ./ ...
            numel(layerOutdegree);
    end
end
realized.parentOutdegreeCV = weighted_nanmean(outdegreeCV, ...
    layerCounts(1:L));
realized.parentUseEffectiveFraction = weighted_nanmean( ...
    parentUseFraction,parentUseWeight);
realized.cvDurationWithinLayer = weighted_nanmean(durationLayerCV,layerCounts);
realized.cvWeightWithinLayer = weighted_nanmean(layerCV, layerCounts);
realized.rhoDurationWeight = weighted_nanmean(layerRho, layerCounts);
taskLayerShare=layerCounts./sum(layerCounts);
realized.taskLayerEntropyNormalized=-sum(taskLayerShare.* ...
    log(taskLayerShare))./log(L+1);
if ~isfinite(realized.taskLayerEntropyNormalized) || ...
        realized.taskLayerEntropyNormalized < -cfg.numericTolerance || ...
        realized.taskLayerEntropyNormalized > 1+cfg.numericTolerance
    error('RepairTask:TaskLayerEntropy', ...
        'Normalized task-layer entropy is outside [0,1].');
end

expectedRatio=caseRow.durationDepthRatio;
if abs(realized.durationDepthRatio-expectedRatio) > ...
        1000*eps(max(1,abs(expectedRatio)))
    error('RepairTask:DurationDepthControl', ...
        'Realized duration-depth ratio does not match its target.');
end
end

function counts = allocate_integer_counts(N, shares, minimumCount)
nLayer = numel(shares);
remaining = N - minimumCount*nLayer;
if remaining < 0
    error('RepairTask:IntegerAllocation', 'N is below the layer minimum.');
end
target = remaining .* shares(:);
extra = floor(target);
left = remaining - sum(extra);
[~, order] = sort(target-extra, 'descend');
extra(order(1:left)) = extra(order(1:left)) + 1;
counts = minimumCount + extra;
end

function x = gamma_unit_mean(stream, cv, n)
if cv <= 1e-12
    x = ones(n,1);
else
    shape = 1 ./ (cv.^2);
    scale = cv.^2;
    x = gamma_from_stream(stream, shape, [n, 1]) .* scale;
    x = max(x, realmin);
end
end

function x = gamma_from_stream(stream, shape, outputSize)
if ~isscalar(shape) && ~isequal(size(shape),outputSize)
    error('RepairTask:GammaShapeSize', ...
        'Gamma shape must be scalar or match the requested output size.');
end
if any(~isfinite(shape(:)) | shape(:)<=0)
    error('RepairTask:InvalidGammaShape', ...
        'Gamma shape parameters must be finite and strictly positive.');
end
u = rand(stream,outputSize);
u = min(max(u,realmin('double')),1-eps('double'));
x = gammaincinv(u,shape,'lower');
end

function reordered = impose_rank_association(duration, values, targetRho, stream)
n = numel(values);
if n < 3 || abs(targetRho) < 1e-12 || std(duration) <= eps
    reordered = values;
    return;
end
rD = tiedrank(duration);
pD = (rD-0.5)./n;
zD = sqrt(2).*erfinv(2.*pD-1);
rW = tiedrank(values);
pW = (rW-0.5)./n;
zW = sqrt(2).*erfinv(2.*pW-1);
zW = zW(randperm(stream,n));
score = targetRho.*zD + sqrt(max(0,1-targetRho.^2)).*zW;
[~, scoreOrder] = sort(score, 'ascend');
sortedValues = sort(values, 'ascend');
reordered = zeros(n,1);
reordered(scoreOrder) = sortedValues;
end

function chosen = weighted_without_replacement(stream, weights, k)
weights = weights(:);
chosen = zeros(k,1);
for j = 1:k
    total = sum(weights);
    u = rand(stream,1,1) .* total;
    idx = find(cumsum(weights) > u, 1, 'first');
    chosen(j) = idx;
    weights(idx) = 0;
end
end

function cv = coefficient_of_variation(x)
m = mean(x);
if m <= 0 || numel(x) < 2
    cv = 0;
else
    cv = std(x,0) ./ m;
end
end

function rho = spearman_pair(x, y)
if numel(x) < 3 || std(x)<=eps || std(y)<=eps
    rho = NaN;
    return;
end
rx = tiedrank(x);
ry = tiedrank(y);
c = corrcoef(rx,ry);
rho = c(1,2);
end

function y = weighted_nanmean(x, w)
ok = isfinite(x) & isfinite(w) & w>0;
if any(ok)
    y = sum(x(ok).*w(ok)) ./ sum(w(ok));
else
    y = NaN;
end
end

function validate_generated_network(net,L,cfg)
assert(net.N==numel(net.duration));
assert(all(net.duration>0 & isfinite(net.duration)));
assert(all(net.duration>=cfg.minimumTaskDuration));
assert(abs(mean(net.duration)-1)<=1000*eps('double'));
assert(all(net.weight>0 & isfinite(net.weight)));
assert(abs(sum(net.weight)-1)<=100*cfg.numericTolerance);
assert(all(net.layer>=0 & net.layer<=L & mod(net.layer,1)==0));

incoming=zeros(net.N,1);
for parent=1:net.N
    child=double(net.children{parent}(:));
    assert(all(child>=1 & child<=net.N));
    assert(numel(unique(child))==numel(child));
    assert(all(net.layer(child)==net.layer(parent)+1));
    incoming(child)=incoming(child)+1;
end
assert(all(incoming(net.layer==0)==0));
assert(all(incoming(net.layer>0)>=1));
assert(all(double(net.parentCount)==incoming));
end
