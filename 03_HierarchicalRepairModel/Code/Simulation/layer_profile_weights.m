function p=layer_profile_weights(profileName,L)
x=((0:L)'+0.5)./(L+1);
switch lower(char(profileName))
    case 'shallow'
        p=beta_kernel(x,1.5,4);
    case 'uniform'
        p=ones(L+1,1);
    case 'intermediate'
        p=beta_kernel(x,3,3);
    case 'deep'
        p=beta_kernel(x,4,1.5);
    case 'deep_tail'
        p=beta_kernel(x,1.5,4);
        p=p/sum(p);
        p=.9*p;
        p(end)=p(end)+.1;
    otherwise
        error('RepairTask:UnknownLayerProfile', ...
            'Unknown layer profile: %s',profileName);
end
p=max(p,realmin); p=p/sum(p);
end

function y=beta_kernel(x,a,b)
y=exp((a-1).*log(x)+(b-1).*log(1-x)- ...
    (gammaln(a)+gammaln(b)-gammaln(a+b)));
end
