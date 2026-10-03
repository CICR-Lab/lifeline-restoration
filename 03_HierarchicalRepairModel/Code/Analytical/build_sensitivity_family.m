function result = build_sensitivity_family(cfg,familyName,limitName)
if nargin < 3 || isempty(limitName)
    limitName = 'frontier';
end
familyName = normalize_family_name(familyName);
limitName = lower(char(limitName));

profiles = service_profiles(3);
basePi = profiles(3).pi;
[names,models] = make_models(cfg,familyName,limitName,basePi);
n = numel(names);
theta = sensitivity_time_grid(cfg,familyName);

values = nan(5,n);
classes = nan(5,n);
curves = cell(n,1);
fits = cell(n,1);
[T,~] = load_empirical_indicator_ranges(cfg.data.empiricalRanges, ...
    cfg.sensitivity.indicatorNames);

for k = 1:n
    sol = analytical_sensitivity_response(theta,models{k});
    u = sol.theta/sol.metrics.tau50;
    keep = u <= cfg.fit.uMax;
    fits{k} = fit_weibull_cdf(u(keep),sol.R(keep),cfg.fit);
    curves{k} = struct('u',u,'R',sol.R,'theta',sol.theta, ...
        'model',models{k});
    values(:,k) = indicator_vector(sol.metrics)';
end

for m = 1:5
    range = get_empirical_range(T,cfg.sensitivity.indicatorNames{m});
    classes(m,:) = classify_empirical_range(values(m,:),range);
end

result.family = familyName;
result.limitName = limitName;
result.names = names;
result.models = models;
result.values = values;
result.classes = classes;
result.curves = curves;
result.fits = fits;
end

function theta = sensitivity_time_grid(cfg,familyName)
theta = cfg.theta.grid;
if strcmp(familyName,'release_duration') && ...
        cfg.sensitivity.releaseDurationTimeMax > theta(end)
    step = theta(2)-theta(1);
    n = ceil(cfg.sensitivity.releaseDurationTimeMax/step)+1;
    theta = linspace(0,cfg.sensitivity.releaseDurationTimeMax,n)';
end
end

function [names,models] = make_models(cfg,familyName,limitName,basePi)
identity = struct('type','identity','q',1,'h',0.70,'s',20);
switch familyName
    case 'dependency'
        names = cfg.sensitivity.accessNames;
        n = numel(names);
        models = cell(n,1);
        for k = 1:n
            models{k} = base_model(cfg,limitName,basePi, ...
                cfg.sensitivity.accessTypes{k},cfg.sensitivity.accessM(k), ...
                ones(1,4),'exponential',1,identity);
        end

    case 'rates'
        names = cfg.sensitivity.rateNames;
        n = numel(names);
        models = cell(n,1);
        for k = 1:n
            rates = cfg.sensitivity.rateProfiles(k,:);
            rates = rates/mean(rates);
            models{k} = base_model(cfg,limitName,basePi,'single',1, ...
                rates,'exponential',1,identity);
        end

    case 'release_duration'
        names = cfg.sensitivity.releaseDurationNames;
        n = numel(names);
        models = cell(n,1);
        for k = 1:n
            models{k} = base_model(cfg,limitName,basePi,'single',1, ...
                ones(1,4),cfg.sensitivity.releaseDurationFamilies{k}, ...
                cfg.sensitivity.releaseDurationCV(k), ...
                cfg.sensitivity.releaseDurationSpecs{k});
        end

    otherwise
        error('Unknown sensitivity family: %s',familyName);
end
end

function model = base_model(cfg,limitName,pi,dependencyType,m,lambda, ...
        durationFamily,durationCV,releaseSpec)
model.L = 3;
model.pi = pi;
model.limitName = limitName;
model.dependencyType = dependencyType;
model.m = m;
model.lambda = lambda;
model.durationFamily = durationFamily;
model.durationCV = durationCV;
model.meanDuration = 1;
model.releaseSpec = releaseSpec;
model.solverCfg = cfg.solver;
end

function familyName = normalize_family_name(familyName)
familyName = lower(char(familyName));
switch familyName
    case 'enabling'
        familyName = 'dependency';
    case 'release'
        familyName = 'release_duration';
end
end
