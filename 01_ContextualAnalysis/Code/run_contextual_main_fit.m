function run_contextual_main_fit(projectRoot, inputPath, runDir, configPath, nReplicates, referenceRunDir)
% Contextual full fits with CR2 coefficient inference and separate R2 bootstraps.

if nargin < 1 || isempty(projectRoot) || strlength(string(projectRoot)) == 0
    moduleRoot = ctx_module_root();
    projectRoot = string(fileparts(fileparts(fileparts(moduleRoot))));
else
    moduleRoot = ctx_module_root();
end
if nargin < 2 || isempty(inputPath) || strlength(string(inputPath)) == 0, inputPath = fullfile(moduleRoot, "outputs", "intermediate", "analysis_input", "input", "canonical_analysis_table.mat"); end
if nargin < 3 || isempty(runDir) || strlength(string(runDir)) == 0, runDir = fullfile(moduleRoot, "outputs", "analysis"); end
if nargin < 4 || isempty(configPath) || strlength(string(configPath)) == 0, configPath = fullfile(moduleRoot, "contextual_config.json"); end
if nargin < 5 || isempty(nReplicates), nReplicates = 2000; end
if nargin < 6 || isempty(referenceRunDir), referenceRunDir = ""; end
projectRoot = string(projectRoot); inputPath = string(inputPath); runDir = string(runDir); configPath = string(configPath); referenceRunDir = string(referenceRunDir);
required = [inputPath, configPath];
if any(~isfile(required)), error("ContextualMainFit:MissingInput", "Required analysis input or configuration is missing."); end
if isfolder(runDir), error("ContextualMainFit:ExistingRun", "Analysis output already exists: %s", runDir); end
for d = ["models", "validation", "figure_data", "supplementary", "manifest", "logs", "figures"], mkdir(fullfile(runDir, d)); end
cfg = jsondecode(fileread(configPath));
seed = double(cfg.random_seed) + 9020;
manifest = struct();
manifest.analysis_id = "contextual_full_fit";
manifest.run_id = char(string(runDir));
manifest.upstream_run_id = "none";
manifest.input_path = char(inputPath);
manifest.input_sha256 = fileHash(inputPath);
manifest.config_sha256 = fileHash(configPath);
manifest.runner_sha256 = fileHash(mfilename("fullpath") + ".m");
manifest.cr2_runner_sha256 = fileHash(fullfile(fileparts(mfilename("fullpath")), "ctx_cr2_satterthwaite.m"));
manifest.random_seed = seed;
manifest.bootstrap_successful_replicates = nReplicates;
manifest.started_at = char(datetime("now", "Format", "yyyy-MM-dd'T'HH:mm:ss"));
manifest.status = "started";
writeJson(fullfile(runDir, "manifest", "run_started.json"), manifest);

cr2Validation = validate_ctx_cr2_satterthwaite();
writeJson(fullfile(runDir, "validation", "cr2_implementation_validation.json"), cr2Validation);

S = load(inputPath, "T");
T = S.T;
d0 = T(T.d0_model_eligible, :);
fullT90 = T(T.t90_context_eligible, :);
aligned = T(T.t90_primary_eligible, :);
assertRecordSet(d0, 985, 187, "full D0 record set");
assertRecordSet(fullT90, 951, 187, "full T90 record set");
assertRecordSet(aligned, 613, 120, "D0-aligned T90 record set");

levelsD0 = incomeLevels(d0, cfg);
levelsFull = incomeLevels(fullT90, cfg);
levelsAligned = incomeLevels(aligned, cfg);
[Xd0, namesD0] = fixedDesign(d0, false, levelsD0);
[Xfull, namesFull] = fixedDesign(fullT90, false, levelsFull);
[Xaligned0, namesAligned0] = fixedDesign(aligned, false, levelsAligned);
[Xaligned1, namesAligned1] = fixedDesign(aligned, true, levelsAligned);

d0Fit = fitFractional(Xd0, d0.outage0);
fullFit = fitLinear(Xfull, log(fullT90.T90));
aligned0Fit = fitLinear(Xaligned0, log(aligned.T90));
aligned1Fit = fitLinear(Xaligned1, log(aligned.T90));

metrics = [metricRow("D0", "system_and_context_full_D0_record_set", d0.outage0, ...
        d0Fit.prediction, height(d0), numel(unique(d0.eid)), "D0 original scale"); ...
    metricRow("T90", "system_and_context_full_T90_record_set", log(fullT90.T90), ...
        fullFit.prediction, height(fullT90), numel(unique(fullT90.eid)), "ln(T90)"); ...
    metricRow("T90", "system_and_context_D0_aligned_T90_record_set", log(aligned.T90), ...
        aligned0Fit.prediction, height(aligned), numel(unique(aligned.eid)), "ln(T90)"); ...
    metricRow("T90", "system_and_context_plus_D0_D0_aligned_T90_record_set", log(aligned.T90), ...
        aligned1Fit.prediction, height(aligned), numel(unique(aligned.eid)), "ln(T90)")];
increment = pairedIncrement(log(aligned.T90), aligned0Fit.prediction, aligned1Fit.prediction);
d0Effects = d0EffectTable(d0Fit.model, Xd0, namesD0, levelsD0);
[fullEffects, fullCr2] = linearEffectTable(fullFit, Xfull, log(fullT90.T90), fullT90.eid, ...
    namesFull, levelsFull, "full_T90_record_set");
[alignedEffects, alignedCr2] = linearEffectTable(aligned1Fit, Xaligned1, log(aligned.T90), ...
    aligned.eid, namesAligned1, levelsAligned, "D0_aligned_T90_record_set");
writeCr2Details(fullCr2, fullfile(runDir, "validation", "cr2_full_T90_covariance.csv"));
writeCr2Details(alignedCr2, fullfile(runDir, "validation", "cr2_aligned_plus_D0_covariance.csv"));

[d0Boot, d0Status, d0Failures] = runMetricBootstrap("D0", d0, levelsD0, nReplicates, ...
    seed + 10000, fullfile(runDir, "validation", "bootstrap_D0_checkpoint.mat"), ...
    cfg.bootstrap_checkpoint_interval);
[fullBoot, fullStatus, fullFailures] = runMetricBootstrap("full_T90", fullT90, levelsFull, ...
    nReplicates, seed + 20000, fullfile(runDir, "validation", "bootstrap_full_T90_checkpoint.mat"), ...
    cfg.bootstrap_checkpoint_interval);
[alignedBoot, alignedStatus, alignedFailures] = runMetricBootstrap("aligned_T90", aligned, ...
    levelsAligned, nReplicates, seed + 30000, ...
    fullfile(runDir, "validation", "bootstrap_aligned_T90_checkpoint.mat"), ...
    cfg.bootstrap_checkpoint_interval);
bootMetrics = combineBootstrapMetrics(d0Boot, fullBoot, alignedBoot);
metrics = attachMetricIntervals(metrics, bootMetrics);
increment = attachIncrementIntervals(increment, bootMetrics);

writetable(metrics, fullfile(runDir, "models", "primary_fitted_performance.csv"));
writetable(increment, fullfile(runDir, "models", "primary_fitted_increment.csv"));
writetable(d0Effects, fullfile(runDir, "models", "d0_full_record_set_effects.csv"));
writetable(fullEffects, fullfile(runDir, "models", "t90_full_record_set_effects.csv"));
writetable(alignedEffects, fullfile(runDir, "models", "t90_D0_aligned_plus_D0_effects.csv"));
writetable(bootMetrics, fullfile(runDir, "validation", "bootstrap_metric_replicates.csv"));
bootstrapStatus = [d0Status; fullStatus; alignedStatus];
writetable(bootstrapStatus, fullfile(runDir, "validation", "bootstrap_status.csv"));
writetable([d0Failures; fullFailures; alignedFailures], ...
    fullfile(runDir, "validation", "bootstrap_failures.csv"));

writeFittedValues(d0, d0Fit.prediction, fullT90, fullFit.prediction, aligned, ...
    aligned0Fit.prediction, aligned1Fit.prediction, runDir);
writeFigureData(metrics, increment, alignedEffects, d0, fullT90, runDir);
loeoPerformancePath = "";
loeoIncrementPath = "";
loeoRobustnessPath = "";
writeSupplementaryTables(metrics, increment, d0Effects, fullEffects, alignedEffects, ...
    loeoPerformancePath, loeoIncrementPath, loeoRobustnessPath, runDir);

qa = struct();
qa.status = "complete";
qa.d0_records = height(d0);
qa.d0_earthquakes = numel(unique(d0.eid));
qa.full_t90_records = height(fullT90);
qa.full_t90_earthquakes = numel(unique(fullT90.eid));
qa.aligned_t90_records = height(aligned);
qa.aligned_t90_earthquakes = numel(unique(aligned.eid));
qa.aligned_design_nested = isequal(namesAligned1(1:end-1), namesAligned0);
qa.sse_plus_not_greater = increment.sse_plus_D0 <= increment.sse_system_and_context + 1e-10;
qa.bootstrap_all_blocks_complete = all(bootstrapStatus.successful_replicates == nReplicates);
qa.cr2_validation_complete = cr2Validation.status == "complete";
qa.reference_comparison_requested = strlength(referenceRunDir) > 0;
qa.reference_comparison_status = "not_requested";
qa.reference_comparison_passed = [];
if qa.reference_comparison_requested
    qa.reference_comparison_passed = comparePreviousEstimates(referenceRunDir, metrics, increment, ...
        fullEffects, alignedEffects);
    if qa.reference_comparison_passed, qa.reference_comparison_status = "passed";
    else, qa.reference_comparison_status = "failed"; end
end
qa.completed_at = char(datetime("now", "Format", "yyyy-MM-dd'T'HH:mm:ss"));
if ~qa.aligned_design_nested || ~qa.sse_plus_not_greater || ~qa.bootstrap_all_blocks_complete || ...
        ~qa.cr2_validation_complete || (qa.reference_comparison_requested && ~qa.reference_comparison_passed)
    qa.status = "failed";
end
writeJson(fullfile(runDir, "manifest", "validation.json"), qa);
if qa.status ~= "complete"
    error("ContextualMainFit:Validation", "Numerical validation failed.");
end

manifest.status = "complete";
manifest.completed_at = qa.completed_at;
manifest.bootstrap_attempts = sum(bootstrapStatus.attempted_replicates);
manifest.bootstrap_failures = sum(bootstrapStatus.failed_attempts);
writeJson(fullfile(runDir, "manifest", "run_complete.json"), manifest);
fprintf("Contextual full-fit inference update complete: %s\n", runDir);
end

function assertRecordSet(D, expectedN, expectedE, label)
if height(D) ~= expectedN || numel(unique(D.eid)) ~= expectedE
    error("ContextualMainFit:RecordSetMismatch", ...
        "%s has %d records and %d earthquakes; expected %d and %d.", ...
        label, height(D), numel(unique(D.eid)), expectedN, expectedE);
end
end

function levels = incomeLevels(D, cfg)
levels = sort(unique(string(D.development_level)));
reference = string(cfg.reference_income_group);
if ~ismember(reference, levels)
    error("ContextualMainFit:IncomeReference", "Reference income group is absent.");
end
levels = [reference; levels(levels ~= reference)];
end

function [X, names] = fixedDesign(D, includeD0, income)
X = [D.avg_PGA / 0.1, (D.year - 2000) / 10, log2(D.population), ...
    double(D.sys == "W"), double(D.sys == "G")];
names = ["pga_0p1g"; "earthquake_year_10y"; "population_doubling"; ...
    "water_vs_electric"; "gas_vs_electric"];
for j = 2:numel(income)
    X(:, end + 1) = double(string(D.development_level) == income(j)); %#ok<AGROW>
    names(end + 1, 1) = "income_" + matlab.lang.makeValidName(income(j)); %#ok<AGROW>
end
if includeD0
    X(:, end + 1) = log2(D.outage0);
    names(end + 1, 1) = "D0_doubling";
end
end

function fit = fitFractional(X, y)
lastwarn("");
model = fitglm(X, y, "linear", "Distribution", "binomial", "Link", "logit", ...
    "Options", statset("MaxIter", 5000));
[warningMessage, warningId] = lastwarn;
if strlength(string(warningId)) > 0
    error("ContextualMainFit:FractionalFitWarning", "%s: %s", warningId, warningMessage);
end
prediction = predict(model, X);
if any(~isfinite(prediction))
    error("ContextualMainFit:FractionalPrediction", "Fractional-logit predictions are non-finite.");
end
fit = struct("model", model, "prediction", prediction);
end

function fit = fitLinear(X, y)
A = [ones(size(X, 1), 1), X];
if rank(A) < size(A, 2), beta = lsqminnorm(A, y); else, beta = A \ y; end
fit = struct("beta", beta, "prediction", A * beta, "rank", rank(A));
end

function T = metricRow(outcome, model, y, prediction, n, e, scale)
m = metric(y, prediction);
T = table(string(outcome), string(model), n, e, string(scale), m.r2, m.rmse, m.mae, ...
    'VariableNames', {'outcome', 'model', 'record_count', 'earthquake_count', ...
    'outcome_scale', 'fitted_r2', 'rmse', 'mae'});
end

function m = metric(y, prediction)
residual = y - prediction;
m = struct("sse", sum(residual .^ 2), "sst", sum((y - mean(y)) .^ 2), ...
    "r2", 1 - sum(residual .^ 2) / sum((y - mean(y)) .^ 2), ...
    "rmse", sqrt(mean(residual .^ 2)), "mae", mean(abs(residual)));
end

function T = pairedIncrement(y, contextPrediction, plusPrediction)
c = metric(y, contextPrediction); p = metric(y, plusPrediction);
T = table(c.r2, p.r2, p.r2 - c.r2, (c.sse - p.sse) / c.sse, c.sse, p.sse, ...
    p.rmse - c.rmse, p.mae - c.mae, ...
    'VariableNames', {'system_and_context_fitted_r2', 'plus_D0_fitted_r2', ...
    'delta_fitted_r2', 'partial_r2_D0', 'sse_system_and_context', ...
    'sse_plus_D0', 'delta_rmse', 'delta_mae'});
end

function T = d0EffectTable(model, X, names, income)
ids = ["water_vs_electric"; "gas_vs_electric"; "pga_0p1g"; ...
    "earthquake_year_10y"; "population_doubling"];
for j = 2:numel(income)
    ids(end + 1, 1) = "income_" + matlab.lang.makeValidName(income(j)); %#ok<AGROW>
end
values = nan(numel(ids), 1);
for j = 1:numel(ids), values(j) = d0MarginalChange(model, X, names, ids(j)); end
T = table(ids, values, repmat("expected D0 difference", numel(ids), 1), ...
    'VariableNames', {'effect_id', 'estimate', 'effect_scale'});
end

function value = d0MarginalChange(model, X, names, id)
X0 = X; X1 = X;
if id == "water_vs_electric" || id == "gas_vs_electric"
    X0(:, names == "water_vs_electric" | names == "gas_vs_electric") = 0;
    X1 = X0; X1(:, names == id) = 1;
elseif startsWith(id, "income_")
    X0(:, startsWith(names, "income_")) = 0;
    X1 = X0; X1(:, names == id) = 1;
else
    X1(:, names == id) = X1(:, names == id) + 1;
end
value = mean(predict(model, X1) - predict(model, X0));
end

function [T, details] = linearEffectTable(fit, X, y, eid, names, income, recordSet)
A = [ones(size(X, 1), 1), X];
contrasts = [zeros(1, numel(names)); eye(numel(names))];
[inference, details] = ctx_cr2_satterthwaite(A, y, eid, contrasts, names);
if max(abs(fit.beta(2:end) - inference.estimate)) > 1e-10
    error("ContextualMainFit:CoefficientMismatch", "CR2 and fitted coefficients differ.");
end
[groups, labels] = effectLabels(names, income);
T = table(names, groups, labels, inference.estimate, inference.se_cr2, ...
    inference.satterthwaite_df, exp(inference.estimate), exp(inference.ci95_low), ...
    exp(inference.ci95_high), inference.p_value, repmat(string(recordSet), numel(names), 1), ...
    repmat("adjusted T90 factor", numel(names), 1), repmat("CR2 Satterthwaite", numel(names), 1), ...
    'VariableNames', {'effect_id', 'group', 'reported_contrast', 'estimate_log_T90', ...
    'se_cr2', 'satterthwaite_df', 'estimate', 'ci95_low', 'ci95_high', 'p_value', ...
    'record_set', 'effect_scale', 'ci_method'});
end

function [groups, labels] = effectLabels(names, income)
groups = strings(numel(names), 1); labels = strings(numel(names), 1);
for j = 1:numel(names)
    id = names(j);
    if id == "water_vs_electric", groups(j) = "Lifeline system"; labels(j) = "Water supply versus electric power";
    elseif id == "gas_vs_electric", groups(j) = "Lifeline system"; labels(j) = "Natural gas versus electric power";
    elseif id == "pga_0p1g", groups(j) = "Service-region context"; labels(j) = "0.1 g service-region mean PGA";
    elseif id == "earthquake_year_10y", groups(j) = "Service-region context"; labels(j) = "Ten-year earthquake year";
    elseif id == "population_doubling", groups(j) = "Population"; labels(j) = "Twofold matched service-region population";
    elseif id == "D0_doubling", groups(j) = "Initial disruption"; labels(j) = "Twofold positive D0";
    elseif startsWith(id, "income_")
        groups(j) = "Service-region income group";
        index = find("income_" + matlab.lang.makeValidName(income) == id, 1);
        labels(j) = income(index) + " versus high income";
    end
end
end

function [values, status, failures] = runMetricBootstrap(kind, D, levels, n, seed, checkpoint, interval)
switch kind
    case {"D0", "full_T90"}, width = 3;
    case "aligned_T90", width = 8;
    otherwise, error("ContextualMainFit:BootstrapKind", "Unknown bootstrap block %s.", kind);
end
if isfile(checkpoint)
    C = load(checkpoint, "values", "successCount", "attemptCount", "failureAttempts", "failureIds", "failureMessages");
    values = C.values; successCount = C.successCount; attemptCount = C.attemptCount;
    failureAttempts = C.failureAttempts; failureIds = string(C.failureIds); failureMessages = string(C.failureMessages);
else
    values = nan(n, width); successCount = 0; attemptCount = 0;
    failureAttempts = zeros(0, 1); failureIds = strings(0, 1); failureMessages = strings(0, 1);
end
if size(values, 1) ~= n || size(values, 2) ~= width
    error("ContextualMainFit:CheckpointShape", "Bootstrap checkpoint shape does not match.");
end
maximumAttempts = max(n * 10, n + 100);
while successCount < n && attemptCount < maximumAttempts
    attemptCount = attemptCount + 1;
    rng(seed + attemptCount, "twister");
    try
        B = D(sampleEventRows(D.eid), :);
        row = bootstrapMetricRow(kind, B, levels);
        successCount = successCount + 1;
        values(successCount, :) = row;
    catch ME
        failureAttempts(end + 1, 1) = attemptCount; %#ok<AGROW>
        failureIds(end + 1, 1) = string(ME.identifier); %#ok<AGROW>
        failureMessages(end + 1, 1) = string(ME.message); %#ok<AGROW>
    end
    if mod(attemptCount, interval) == 0 || successCount == n
        save(checkpoint, "values", "successCount", "attemptCount", "failureAttempts", ...
            "failureIds", "failureMessages", "-v7.3");
        fprintf("%s bootstrap: %d/%d successful (%d attempts).\n", kind, successCount, n, attemptCount);
    end
end
if successCount ~= n
    error("ContextualMainFit:BootstrapIncomplete", ...
        "%s bootstrap obtained %d of %d successful replicates after %d attempts.", ...
        kind, successCount, n, attemptCount);
end
status = table(string(kind), n, successCount, attemptCount, numel(failureAttempts), ...
    'VariableNames', {'analysis_block', 'requested_successful_replicates', ...
    'successful_replicates', 'attempted_replicates', 'failed_attempts'});
failures = table(repmat(string(kind), numel(failureAttempts), 1), failureAttempts, failureIds, failureMessages, ...
    'VariableNames', {'analysis_block', 'attempt', 'error_id', 'error_message'});
end

function row = bootstrapMetricRow(kind, D, levels)
switch kind
    case "D0"
        [X, ~] = fixedDesign(D, false, levels); fit = fitFractional(X, D.outage0);
        m = metric(D.outage0, fit.prediction); row = [m.r2 m.rmse m.mae];
    case "full_T90"
        [X, ~] = fixedDesign(D, false, levels); fit = fitLinear(X, log(D.T90));
        m = metric(log(D.T90), fit.prediction); row = [m.r2 m.rmse m.mae];
    case "aligned_T90"
        [X0, ~] = fixedDesign(D, false, levels); [X1, ~] = fixedDesign(D, true, levels);
        f0 = fitLinear(X0, log(D.T90)); f1 = fitLinear(X1, log(D.T90));
        m0 = metric(log(D.T90), f0.prediction); m1 = metric(log(D.T90), f1.prediction);
        row = [m0.r2 m1.r2 m1.r2-m0.r2 (m0.sse-m1.sse)/m0.sse ...
            m0.rmse m1.rmse m0.mae m1.mae];
end
end

function T = combineBootstrapMetrics(d0, full, aligned)
n = size(d0, 1);
T = table(d0(:,1), d0(:,2), d0(:,3), full(:,1), full(:,2), full(:,3), ...
    aligned(:,1), aligned(:,2), aligned(:,3), aligned(:,4), aligned(:,5), aligned(:,6), ...
    aligned(:,7), aligned(:,8), (1:n)', ...
    'VariableNames', {'d0_r2','d0_rmse','d0_mae','full_t90_r2','full_t90_rmse', ...
    'full_t90_mae','aligned_context_r2','aligned_plus_r2','delta_r2','partial_r2', ...
    'aligned_context_rmse','aligned_plus_rmse','aligned_context_mae','aligned_plus_mae','replicate'});
end

function idx = sampleEventRows(eid)
events = unique(eid, "stable"); draw = events(randi(numel(events), numel(events), 1));
idx = zeros(0, 1);
for k = 1:numel(draw), idx = [idx; find(eid == draw(k))]; %#ok<AGROW>
end
end

function T = attachMetricIntervals(T, B)
prefix = ["d0"; "full_t90"; "aligned_context"; "aligned_plus"];
for j = 1:height(T)
    for metricName = ["r2", "rmse", "mae"]
        q = prctile(B.(prefix(j) + "_" + metricName), [2.5 97.5]);
        T.(metricName + "_ci95_low")(j) = q(1);
        T.(metricName + "_ci95_high")(j) = q(2);
    end
end
end

function T = attachIncrementIntervals(T, B)
for source = ["delta_r2", "partial_r2"]
    q = prctile(B.(source), [2.5 97.5]);
    T.(source + "_ci95_low") = q(1); T.(source + "_ci95_high") = q(2);
end
end

function writeFittedValues(d0, d0p, fullT90, fullp, aligned, aligned0, aligned1, runDir)
A = d0(:, ["record_id", "eid", "outage0"]); A.fitted_D0 = d0p;
writetable(A, fullfile(runDir, "models", "d0_full_record_set_fitted_values.csv"));
B = fullT90(:, ["record_id", "eid", "T90"]); B.fitted_log_T90 = fullp; B.fitted_T90 = exp(fullp);
writetable(B, fullfile(runDir, "models", "t90_full_record_set_fitted_values.csv"));
C = aligned(:, ["record_id", "eid", "T90", "outage0"]);
C.fitted_system_and_context_log_T90 = aligned0; C.fitted_plus_D0_log_T90 = aligned1;
writetable(C, fullfile(runDir, "models", "t90_D0_aligned_fitted_values.csv"));
end

function writeFigureData(metrics, increment, alignedEffects, d0, fullT90, runDir)
C = metrics(:, ["outcome", "model", "record_count", "earthquake_count", ...
    "outcome_scale", "fitted_r2", "r2_ci95_low", "r2_ci95_high"]);
C.delta_fitted_r2 = [nan; nan; nan; increment.delta_fitted_r2];
C.delta_ci95_low = [nan; nan; nan; increment.delta_r2_ci95_low];
C.delta_ci95_high = [nan; nan; nan; increment.delta_r2_ci95_high];
writetable(C, fullfile(runDir, "figure_data", "figure2c_fitted_r2.csv"));
D = alignedEffects; D.order = (1:height(D))';
writetable(D, fullfile(runDir, "figure_data", "figure2d_adjusted_T90_factors.csv"));
A = d0(:, ["record_id", "eid", "sys", "outage0"]);
writetable(A, fullfile(runDir, "figure_data", "figure2a_D0_records.csv"));
B = fullT90(:, ["record_id", "eid", "sys", "T90"]);
writetable(B, fullfile(runDir, "figure_data", "figure2b_T90_records.csv"));
end

function writeSupplementaryTables(metrics, increment, d0Effects, fullEffects, alignedEffects, ...
        loeoPerformancePath, loeoIncrementPath, loeoRobustnessPath, runDir)
writetable(metrics, fullfile(runDir, "supplementary", "table_fitted_performance.csv"));
writetable(increment, fullfile(runDir, "supplementary", "table_D0_increment.csv"));
writetable(d0Effects, fullfile(runDir, "supplementary", "table_D0_adjusted_effects.csv"));
writetable(fullEffects, fullfile(runDir, "supplementary", "table_T90_full_record_set_adjusted_effects.csv"));
writetable(alignedEffects, fullfile(runDir, "supplementary", "table_T90_D0_aligned_adjusted_effects.csv"));
if strlength(string(loeoPerformancePath)) > 0 && isfile(loeoPerformancePath), copyfile(loeoPerformancePath, fullfile(runDir, "supplementary", "existing_LOEO_performance.csv")); end
if strlength(string(loeoIncrementPath)) > 0 && isfile(loeoIncrementPath), copyfile(loeoIncrementPath, fullfile(runDir, "supplementary", "existing_LOEO_increment.csv")); end
if strlength(string(loeoRobustnessPath)) > 0 && isfile(loeoRobustnessPath), copyfile(loeoRobustnessPath, fullfile(runDir, "supplementary", "existing_LOEO_robustness.csv")); end
end

function writeCr2Details(details, path)
writematrix(details.covariance_cr2, path);
end

function tf = comparePreviousEstimates(oldDir, metrics, increment, fullEffects, alignedEffects)
oldDir = string(oldDir);
if ~isfolder(oldDir), error("ContextualMainFit:ReferenceRunMissing", "Reference run directory is missing: %s", oldDir); end
oldMetrics = readtable(fullfile(oldDir, "models", "primary_fitted_performance.csv"), "TextType", "string");
oldIncrement = readtable(fullfile(oldDir, "models", "primary_fitted_increment.csv"));
oldFull = readtable(fullfile(oldDir, "models", "t90_full_record_set_effects.csv"), "TextType", "string");
oldAligned = readtable(fullfile(oldDir, "models", "t90_D0_aligned_plus_D0_effects.csv"), "TextType", "string");
tf = max(abs(metrics.fitted_r2 - oldMetrics.fitted_r2)) < 1e-10 && ...
    abs(increment.delta_fitted_r2 - oldIncrement.delta_fitted_r2) < 1e-10 && ...
    max(abs(fullEffects.estimate - oldFull.estimate)) < 1e-10 && ...
    max(abs(alignedEffects.estimate - oldAligned.estimate)) < 1e-10;
end

function writeJson(path, value)
fid = fopen(path, "w");
if fid < 0, error("ContextualMainFit:JsonWrite", "Cannot write %s", path); end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, jsonencode(value, "PrettyPrint", true), "char");
end

function hash = fileHash(path)
fid = fopen(path, "r");
if fid < 0, error("ContextualMainFit:HashRead", "Cannot read %s", path); end
cleanup = onCleanup(@() fclose(fid)); bytes = fread(fid, Inf, "*uint8")';
digest = java.security.MessageDigest.getInstance("SHA-256"); digest.update(bytes);
hash = lower(reshape(dec2hex(typecast(digest.digest(), "uint8"))', 1, []));
end
