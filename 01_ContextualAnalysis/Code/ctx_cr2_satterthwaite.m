function [result, details] = ctx_cr2_satterthwaite(X, y, cluster, contrasts, labels)
% CR2 cluster-robust inference for Gaussian linear-model coefficients.
% X must include the intercept. Each column of contrasts defines one test.

arguments
    X double
    y double
    cluster
    contrasts double = eye(size(X, 2))
    labels string = "contrast_" + string(1:size(contrasts, 2))
end

y = y(:);
cluster = cluster(:);
labels = labels(:);
[n, p] = size(X);
if numel(y) ~= n || numel(cluster) ~= n
    error("ContextualCR2:RowMismatch", "X, y and cluster must have the same number of rows.");
end
if size(contrasts, 1) ~= p || size(contrasts, 2) ~= numel(labels)
    error("ContextualCR2:ContrastShape", "Contrasts or labels have incompatible dimensions.");
end
if any(~isfinite(X), "all") || any(~isfinite(y))
    error("ContextualCR2:NonFinite", "X and y must be finite.");
end
if rank(X) ~= p
    error("ContextualCR2:RankDeficient", "The Gaussian-model design matrix is rank deficient.");
end

xtxInverse = inv(X' * X);
beta = xtxInverse * (X' * y);
residual = y - X * beta;
hat = X * xtxInverse * X';
residualMaker = eye(n) - hat;
residualMaker = (residualMaker + residualMaker') / 2;

[group, groupLabels] = findgroups(cluster);
gCount = max(group);
adjustments = cell(gCount, 1);
minimumEigenvalue = inf;
maximumAdjustmentEigenvalue = 0;
meat = zeros(p, p);
for g = 1:gCount
    idx = group == g;
    Xg = X(idx, :);
    eg = residual(idx);
    block = eye(sum(idx)) - Xg * xtxInverse * Xg';
    block = (block + block') / 2;
    [vectors, values] = eig(block, "vector");
    scale = max(1, max(abs(values)));
    tolerance = max(size(block)) * eps(scale) * 100;
    minimumEigenvalue = min(minimumEigenvalue, min(values));
    if any(values <= tolerance)
        error("ContextualCR2:SingularAdjustment", ...
            "CR2 adjustment is singular for cluster %s (minimum eigenvalue %.3g).", ...
            string(groupLabels(g)), min(values));
    end
    adjustment = vectors * diag(1 ./ sqrt(values)) * vectors';
    adjustment = (adjustment + adjustment') / 2;
    adjustments{g} = adjustment;
    maximumAdjustmentEigenvalue = max(maximumAdjustmentEigenvalue, max(eig(adjustment)));
    adjustedResidual = adjustment * eg;
    meat = meat + Xg' * (adjustedResidual * adjustedResidual') * Xg;
end
covariance = xtxInverse * meat * xtxInverse;
covariance = (covariance + covariance') / 2;

q = size(contrasts, 2);
estimate = nan(q, 1);
standardError = nan(q, 1);
df = nan(q, 1);
ciLow = nan(q, 1);
ciHigh = nan(q, 1);
pValue = nan(q, 1);
for j = 1:q
    contrast = contrasts(:, j);
    estimate(j) = contrast' * beta;
    variance = contrast' * covariance * contrast;
    if ~(isfinite(variance) && variance > 0)
        error("ContextualCR2:InvalidVariance", "Contrast %s has invalid CR2 variance.", labels(j));
    end
    standardError(j) = sqrt(variance);

    embedded = zeros(n, gCount);
    for g = 1:gCount
        idx = group == g;
        Xg = X(idx, :);
        embedded(idx, g) = adjustments{g} * Xg * xtxInverse * contrast;
    end
    projected = residualMaker * embedded;
    gram = projected' * projected;
    traceB = sum(diag(gram));
    traceB2 = sum(gram .^ 2, "all");
    df(j) = traceB ^ 2 / traceB2;
    if ~(isfinite(df(j)) && df(j) > 0)
        error("ContextualCR2:InvalidDf", "Contrast %s has invalid Satterthwaite degrees of freedom.", labels(j));
    end
    critical = tinv(0.975, df(j));
    ciLow(j) = estimate(j) - critical * standardError(j);
    ciHigh(j) = estimate(j) + critical * standardError(j);
    pValue(j) = 2 * tcdf(-abs(estimate(j) / standardError(j)), df(j));
end

result = table(labels, estimate, standardError, df, ciLow, ciHigh, pValue, ...
    'VariableNames', {'term', 'estimate', 'se_cr2', 'satterthwaite_df', ...
    'ci95_low', 'ci95_high', 'p_value'});
details = struct();
details.n = n;
details.p = p;
details.cluster_count = gCount;
details.beta = beta;
details.covariance_cr2 = covariance;
details.minimum_adjustment_block_eigenvalue = minimumEigenvalue;
details.maximum_adjustment_eigenvalue = maximumAdjustmentEigenvalue;
details.covariance_minimum_eigenvalue = min(eig(covariance));
details.working_model = "identity";
details.small_sample_correction = "CR2";
details.degrees_of_freedom = "Satterthwaite";
end
