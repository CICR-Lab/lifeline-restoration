function validation = validate_ctx_cr2_satterthwaite()
% Deterministic checks for the contextual CR2 implementation.

y = [1.2; 2.1; 2.7; 4.4; 5.0; 7.3; 8.1];
n = numel(y);
[single, ~] = ctx_cr2_satterthwaite(ones(n, 1), y, (1:n)', 1, "intercept");
expectedEstimate = mean(y);
expectedSe = std(y, 0) / sqrt(n);
expectedDf = n - 1;

rng(20260902, "twister");
cluster = repelem((1:18)', 3);
x = randn(numel(cluster), 1);
z = randn(numel(cluster), 1);
clusterEffect = repelem(randn(18, 1), 3);
outcome = 0.7 + 0.35 * x - 0.2 * z + clusterEffect + 0.3 * randn(numel(cluster), 1);
X = [ones(numel(cluster), 1), x, z];
C = [0; 1; 0];
[base, baseDetails] = ctx_cr2_satterthwaite(X, outcome, cluster, C, "x");

scaledX = [ones(numel(cluster), 1), 7 * x, z];
scaledC = [0; 7; 0];
scaled = ctx_cr2_satterthwaite(scaledX, outcome, cluster, scaledC, "x");

permutation = randperm(numel(cluster));
permuted = ctx_cr2_satterthwaite(X(permutation, :), outcome(permutation), ...
    cluster(permutation), C, "x");

tolerance = 1e-10;
checks = struct();
checks.singleton_estimate = abs(single.estimate - expectedEstimate) < tolerance;
checks.singleton_se = abs(single.se_cr2 - expectedSe) < tolerance;
checks.singleton_df = abs(single.satterthwaite_df - expectedDf) < tolerance;
checks.scaling_estimate = abs(base.estimate - scaled.estimate) < tolerance;
checks.scaling_se = abs(base.se_cr2 - scaled.se_cr2) < tolerance;
checks.scaling_df = abs(base.satterthwaite_df - scaled.satterthwaite_df) < tolerance;
checks.permutation_estimate = abs(base.estimate - permuted.estimate) < tolerance;
checks.permutation_se = abs(base.se_cr2 - permuted.se_cr2) < tolerance;
checks.permutation_df = abs(base.satterthwaite_df - permuted.satterthwaite_df) < tolerance;
checks.covariance_symmetric = norm(baseDetails.covariance_cr2 - baseDetails.covariance_cr2', "fro") < tolerance;
checks.covariance_psd = baseDetails.covariance_minimum_eigenvalue > -tolerance;

names = string(fieldnames(checks));
passed = cellfun(@(name) checks.(name), cellstr(names));
if ~all(passed)
    error("ContextualCR2:ValidationFailure", "CR2 validation failed: %s", strjoin(names(~passed), ", "));
end
validation = struct();
validation.status = "complete";
validation.checks = checks;
validation.singleton_expected_se = expectedSe;
validation.singleton_observed_se = single.se_cr2;
validation.singleton_expected_df = expectedDf;
validation.singleton_observed_df = single.satterthwaite_df;
validation.synthetic_cluster_count = baseDetails.cluster_count;
end
