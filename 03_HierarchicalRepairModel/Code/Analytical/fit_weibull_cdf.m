function fit = fit_weibull_cdf(u, R, fitCfg)
if nargin < 3 || isempty(fitCfg)
    fitCfg.uMax = 4; fitCfg.nStarts = 25; fitCfg.maxIter = 2500;
end
u = u(:); R = R(:);
keep = isfinite(u) & isfinite(R) & u >= 0 & u <= fitCfg.uMax;
u = u(keep); R = min(max(R(keep),0),1);
if numel(u) < 10
    fit = struct('a',NaN,'b',NaN,'mae',NaN,'rmse',NaN,'Rhat',nan(size(u)));
    return;
end
if exist("lsqcurvefit", "file") ~= 2
    error("HierarchicalRepairModel:MissingOptimizer", ...
        "Weibull curve fitting requires lsqcurvefit from MATLAB Optimization Toolbox.");
end
startB = logspace(log10(0.7),log10(6),fitCfg.nStarts);
bestConverged = struct('p',[NaN NaN],'SSE',Inf,'exitflag',-1);
bestAny = bestConverged;
opts = optimoptions("lsqcurvefit","Display","off", ...
    "MaxIterations",fitCfg.maxIter, ...
    "MaxFunctionEvaluations",4*fitCfg.maxIter, ...
    "FunctionTolerance",1e-10, ...
    "StepTolerance",1e-10, ...
    "OptimalityTolerance",1e-8);
lb = [1e-12 1e-12];
ub = [Inf Inf];
for k = 1:numel(startB)
    b0 = startB(k);
    a0 = 1/(log(2)^(1/b0));
    try
        [p,~,residual,exitflag] = lsqcurvefit( ...
            @(p,x) weibull_eval(x,p(1),p(2)), ...
            [a0 b0],u,R,lb,ub,opts);
        if any(~isfinite(residual))
            SSE = Inf;
        else
            SSE = sum(residual.^2);
        end
        if isfinite(SSE) && SSE < bestAny.SSE
            bestAny = struct('p',p,'SSE',SSE,'exitflag',exitflag);
        end
        if exitflag > 0 && isfinite(SSE) && SSE < bestConverged.SSE
            bestConverged = struct('p',p,'SSE',SSE,'exitflag',exitflag);
        end
    catch ME
        if k == numel(startB) && ~isfinite(bestAny.SSE)
            rethrow(ME);
        end
    end
end
if isfinite(bestConverged.SSE)
    best = bestConverged;
else
    best = bestAny;
end
if ~isfinite(best.SSE)
    fit = struct('a',NaN,'b',NaN,'mae',NaN,'rmse',NaN, ...
        'u',u,'R',R,'Rhat',nan(size(u)));
    return;
end
a=best.p(1); b=best.p(2); Rhat=weibull_eval(u,a,b);
fit.a=a; fit.b=b; fit.mae=mean(abs(R-Rhat)); fit.rmse=sqrt(mean((R-Rhat).^2));
fit.u=u; fit.R=R; fit.Rhat=Rhat;
end

function y = weibull_eval(u,a,b)
y = 1-exp(-(max(u,0)./a).^b);
y = min(max(y,0),1);
end
