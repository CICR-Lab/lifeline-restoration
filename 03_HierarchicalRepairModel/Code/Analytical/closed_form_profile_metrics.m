function met = closed_form_profile_metrics(theta, X, pi, limitName)
if nargin < 4 || isempty(limitName)
    limitName = 'frontier';
end
pi = pi(:)';
L = numel(pi)-1;
if size(X,2) ~= L+1
    error('X must contain L+1 analytical layer-response columns.');
end
R = X*pi(:);
R = min(max(cummax(R),0),1);
met = restoration_curve_metrics(theta,R,'linear');
switch lower(char(limitName))
    case 'frontier'
        layerArea = 1:(L+1);
    case 'parallel'
        layerArea = arrayfun(@(n) sum(1./(1:n)),1:(L+1));
    otherwise
        error('Unknown repair limit: %s',limitName);
end
met.areaInfinite = sum(pi .* layerArea);
met.kappaInfinite = met.areaInfinite/met.tau50;
met.eta90Infinite = profile_tail_area(pi,met.tau90,limitName)/met.areaInfinite;
met.effectiveDepth = sum((0:L).*pi);
met.R = R;
end
