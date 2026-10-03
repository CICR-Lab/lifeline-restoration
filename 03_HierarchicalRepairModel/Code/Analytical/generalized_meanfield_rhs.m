function dx = generalized_meanfield_rhs(~, x, model)
L = model.L;
dx = zeros(L+1,1);
lambda = model.lambda(:);
dx(1) = lambda(1) * max(0, 1-x(1));
for l = 2:(L+1)
    prev = min(max(x(l-1),0),1);
    m = model.m(min(l,numel(model.m)));
    switch lower(char(model.accessType))
        case 'single'
            a = prev;
        case 'mandatory'
            a = prev.^m;
        case 'alternative'
            a = 1 - (1-prev).^m;
        otherwise
            error('Unknown accessType: %s', model.accessType);
    end
    dx(l) = lambda(l) * max(0, a-x(l));
end
end
