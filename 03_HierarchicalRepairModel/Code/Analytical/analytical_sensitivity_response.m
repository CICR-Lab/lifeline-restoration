function sol = analytical_sensitivity_response(t,model)
t = t(:);
L = model.L;
pi = model.pi(:);
pi = pi/sum(pi);
lambda = model.lambda(:);

if numel(pi) ~= L+1 || numel(lambda) ~= L+1
    error('model.pi and model.lambda must contain L+1 entries.');
end

limitName = lower(char(model.limitName));
dependencyType = lower(char(model.dependencyType));
durationFamily = lower(char(model.durationFamily));

if ismember(durationFamily,{'gamma','lognormal'}) && ...
        ~strcmp(dependencyType,'single')
    error('Non-exponential S6 scenarios require single prerequisites.');
end

switch durationFamily
    case 'exponential'
        task = zeros(numel(t),L+1);
        for l = 0:L
            task(:,l+1) = 1-exp(-lambda(l+1)*t);
        end
    case 'gamma'
        if any(abs(lambda-lambda(1))>1e-12)
            error('Gamma-duration S6 scenarios require a common mean rate.');
        end
        alpha = 1/model.durationCV^2;
        scale = model.meanDuration/alpha;
        baseCdf = gammainc(t/scale,alpha,'lower');
        task = repmat(baseCdf,1,L+1);
    case 'lognormal'
        baseCdf = duration_cdf(t,'lognormal',model.meanDuration, ...
            model.durationCV);
        task = repmat(baseCdf,1,L+1);
    otherwise
        error('Unsupported S6 duration family: %s',durationFamily);
end

connected = zeros(size(task));
switch limitName
    case 'frontier'
        if strcmp(durationFamily,'exponential')
            odeModel = struct('L',L,'lambda',lambda,'pi',pi, ...
                'accessType',dependencyType,'m',model.m*ones(1,L+1), ...
                'release',struct('type','identity'));
            solverCfg = model.solverCfg;
            solverCfg.tMax = t(end);
            solverCfg.nEval = numel(t);
            odeSol = solve_generalized_meanfield(odeModel,solverCfg);
            connected = odeSol.X;
            task = connected;
        elseif strcmp(durationFamily,'gamma')
            alpha = 1/model.durationCV^2;
            scale = model.meanDuration/alpha;
            for l = 0:L
                connected(:,l+1) = gammainc(t/scale,(l+1)*alpha,'lower');
            end
            task = connected;
        else
            refinement = 4;
            fineT = linspace(t(1),t(end),refinement*(numel(t)-1)+1)';
            fineBaseCdf = duration_cdf(fineT,'lognormal', ...
                model.meanDuration,model.durationCV);
            fineConnected = iid_sum_cdf_fft(fineBaseCdf,L+1);
            connected = interp1(fineT,fineConnected,t,'linear');
            task = connected;
        end

    case 'parallel'
        connected(:,1) = task(:,1);
        for l = 2:(L+1)
            previous = connected(:,l-1);
            switch dependencyType
                case 'single'
                    upstream = previous;
                case 'mandatory'
                    upstream = previous.^model.m;
                case 'alternative'
                    upstream = 1-(1-previous).^model.m;
                otherwise
                    error('Unknown dependency type: %s',dependencyType);
            end
            connected(:,l) = task(:,l).*upstream;
        end

    otherwise
        error('Unknown repair limit: %s',limitName);
end

restored = zeros(size(connected));
for l = 0:L
    restored(:,l+1) = service_release(connected(:,l+1), ...
        model.releaseSpec,l+1);
end
R = restored*pi;
R = min(max(cummax(R),0),1);

sol.theta = t;
sol.X = task;
sol.Y = connected;
sol.H = restored;
sol.R = R;
sol.metrics = restoration_curve_metrics(t,R,'linear');
end
