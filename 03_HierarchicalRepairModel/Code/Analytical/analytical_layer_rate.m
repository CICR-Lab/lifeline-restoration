function V = analytical_layer_rate(theta, L)
theta = theta(:);
V = zeros(numel(theta), L+1);
for l = 0:L
    V(:,l+1) = exp(-theta) .* theta.^l / factorial(l);
end
end
