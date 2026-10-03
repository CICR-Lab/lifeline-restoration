function pi = service_distribution_for_mean_depth(L, depthNorm)

if L < 1 || L ~= floor(L)
    error('L must be a positive integer.');
end
if ~isscalar(depthNorm) || depthNorm < 0 || depthNorm > 1
    error('depthNorm must lie in [0,1].');
end
l = 0:L;
if depthNorm <= 1e-12
    pi = zeros(1,L+1); pi(1) = 1; return;
elseif depthNorm >= 1-1e-12
    pi = zeros(1,L+1); pi(end) = 1; return;
end

target = depthNorm * L;
lo = -60; hi = 60;
for iter = 1:120
    mid = 0.5*(lo+hi);
    z = mid*l;
    z = z - max(z);
    p = exp(z); p = p/sum(p);
    mu = sum(l.*p);
    if mu < target
        lo = mid;
    else
        hi = mid;
    end
end
alpha = 0.5*(lo+hi);
z = alpha*l; z = z-max(z);
pi = exp(z); pi = pi/sum(pi);
end
