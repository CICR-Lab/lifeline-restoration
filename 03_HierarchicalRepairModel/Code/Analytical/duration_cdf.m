function F = duration_cdf(t,distName,meanDuration,cv)

if ~isnumeric(t) || ~isreal(t)
    error('t must be a real numeric array.');
end
if ~isscalar(meanDuration) || ~isfinite(meanDuration) || meanDuration<=0
    error('meanDuration must be a positive finite scalar.');
end
if ~isscalar(cv) || ~isfinite(cv) || cv<0
    error('cv must be a non-negative finite scalar.');
end

t=max(t,0);
name=lower(string(distName));

switch name
    case "exponential"
        lambda=1/meanDuration;
        F=1-exp(-lambda*t);

    case "deterministic"
        F=double(t>=meanDuration);

    case "gamma"
        if cv<=1e-12
            F=double(t>=meanDuration);
        else
            shape=1/(cv^2);
            scale=meanDuration/shape;
            F=gammainc(t/scale,shape,'lower');
        end

    case "lognormal"
        if cv<=1e-12
            F=double(t>=meanDuration);
        else
            sigma2=log(1+cv^2);
            sigma=sqrt(sigma2);
            mu=log(meanDuration)-sigma2/2;
            z=(log(max(t,realmin))-mu)/(sigma*sqrt(2));
            F=0.5*(1+erf(z));
            F(t<=0)=0;
        end

    case "weibull"
        if cv<=1e-12
            F=double(t>=meanDuration);
        else
            k=weibull_shape_from_cv(cv);
            scale=meanDuration/gamma(1+1/k);
            F=1-exp(-(t/scale).^k);
        end

    otherwise
        error('Unknown duration distribution: %s',distName);
end

F=min(max(F,0),1);
end

function k = weibull_shape_from_cv(targetCV)
obj=@(logk) (sqrt(gamma(1+2/exp(logk))/ ...
    gamma(1+1/exp(logk))^2-1)-targetCV).^2;
logk=fminsearch(obj,log(1.5),optimset('Display','off'));
k=exp(logk);
end
