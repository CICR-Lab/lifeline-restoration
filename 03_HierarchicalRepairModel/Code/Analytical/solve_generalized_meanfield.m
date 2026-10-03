function sol = solve_generalized_meanfield(model, solverCfg)
if ~isstruct(model) || ~isscalar(model)
    error('model must be a scalar structure.');
end
if ~isfield(model,'L') || ~isscalar(model.L) || ~isfinite(model.L) || ...
        model.L < 0 || model.L ~= floor(model.L)
    error('model.L must be a non-negative integer.');
end
L = model.L;
requiredFields = {'lambda','accessType','release','pi'};
for k = 1:numel(requiredFields)
    if ~isfield(model,requiredFields{k})
        error('Missing model field: %s.',requiredFields{k});
    end
end

lambda = model.lambda(:);
if numel(lambda) ~= L+1 || any(~isfinite(lambda)) || any(lambda <= 0)
    error('model.lambda must contain L+1 positive finite entries.');
end
if ~isfield(model,'m') || isempty(model.m)
    model.m = ones(1,L+1);
end
if ~isnumeric(model.m) || any(~isfinite(model.m(:))) || ...
        any(model.m(:) < 1) || any(model.m(:) ~= floor(model.m(:)))
    error('model.m must contain positive integers.');
end

pi = model.pi(:);
if numel(pi) ~= L+1 || any(~isfinite(pi)) || any(pi < 0) || sum(pi) <= 0
    error('model.pi must contain L+1 non-negative finite weights with positive sum.');
end
pi = pi / sum(pi);
model.lambda = lambda;
model.pi = pi;

solverFields = {'relTol','absTol','tMax','nEval'};
for k = 1:numel(solverFields)
    if ~isfield(solverCfg,solverFields{k})
        error('Missing solverCfg field: %s.',solverFields{k});
    end
end
if ~isscalar(solverCfg.relTol) || ~isfinite(solverCfg.relTol) || solverCfg.relTol <= 0 || ...
        ~isscalar(solverCfg.absTol) || ~isfinite(solverCfg.absTol) || solverCfg.absTol <= 0
    error('Solver tolerances must be positive finite scalars.');
end
if ~isscalar(solverCfg.tMax) || ~isfinite(solverCfg.tMax) || solverCfg.tMax <= 0
    error('solverCfg.tMax must be a positive finite scalar.');
end
if ~isscalar(solverCfg.nEval) || ~isfinite(solverCfg.nEval) || ...
        solverCfg.nEval < 2 || solverCfg.nEval ~= floor(solverCfg.nEval)
    error('solverCfg.nEval must be an integer greater than or equal to 2.');
end

opts = odeset('RelTol',solverCfg.relTol,'AbsTol',solverCfg.absTol, ...
    'NonNegative',1:(L+1));
tEval = linspace(0,solverCfg.tMax,solverCfg.nEval)';
[theta,X] = ode45(@(t,x) generalized_meanfield_rhs(t,x,model), ...
    tEval, zeros(L+1,1), opts);
X = min(max(X,0),1);
H = zeros(size(X));
for l = 0:L
    H(:,l+1) = service_release(X(:,l+1),model.release,l+1);
end
R = H * pi;
R = min(max(cummax(R),0),1);
sol.theta = theta;
sol.X = X;
sol.H = H;
sol.R = R;
sol.metrics = restoration_curve_metrics(theta,R,'linear');
end
