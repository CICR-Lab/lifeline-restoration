function Y = service_release(X, spec, layerIndex)
if nargin < 3
    layerIndex = 1;
end
if ~isfield(spec,'type')
    spec.type = 'identity';
end
kind = lower(char(spec.type));

switch kind
    case 'identity'
        Y = X;
    case 'power'
        q = pick_parameter(spec, 'q', layerIndex, 1.0);
        Y = X.^q;
    case 'threshold'
        h = pick_parameter(spec, 'h', layerIndex, 0.70);
        s = pick_parameter(spec, 's', layerIndex, 20);
        lo = 1 ./ (1 + exp(s*h));
        hi = 1 ./ (1 + exp(-s*(1-h)));
        raw = 1 ./ (1 + exp(-s*(X-h)));
        Y = (raw - lo) ./ (hi - lo);
    otherwise
        error('Unknown service-release type: %s', kind);
end
Y = min(max(Y,0),1);
end

function value = pick_parameter(spec, fieldName, idx, defaultValue)
if ~isfield(spec, fieldName) || isempty(spec.(fieldName))
    value = defaultValue;
else
    v = spec.(fieldName);
    if isscalar(v)
        value = v;
    else
        value = v(min(idx,numel(v)));
    end
end
end
