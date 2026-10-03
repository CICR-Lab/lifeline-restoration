
function make_TableS4_loss_compression_specification_sensitivity(supplementaryRoot)

if nargin < 1 || isempty(supplementaryRoot)
    scriptDir = string(fileparts(mfilename("fullpath")));
    moduleRoot = scriptDir;
    emergentRoot = string(fileparts(fileparts(moduleRoot)));
    supplementaryRoot = fullfile(emergentRoot, "outputs", "supplementary");
end
rng(1, "twister");


scriptDir = string(fileparts(mfilename("fullpath")));
moduleRoot = scriptDir;
emergentRoot = string(fileparts(fileparts(moduleRoot)));
intermediateRoot = fullfile( ...
    emergentRoot, ...
    "outputs", ...
    "intermediate");

addpath(scriptDir);

suppRoot = fullfile(supplementaryRoot, "source_data");

ensureDir(suppRoot);

inputMat = fullfile( ...
    intermediateRoot, ...
    "02_loss_compression_late_loss", ...
    "loss_tail_dataset.mat");

restorationInputMat = fullfile( ...
    intermediateRoot, ...
    "01_restoration_shape_phase_ratio", ...
    "unique3_records.mat");

mainSourcePath = fullfile( ...
    intermediateRoot, ...
    "02_loss_compression_late_loss", ...
    "SourceData_Fig3c_loss_compression.xlsx");

outXlsx = fullfile( ...
    suppRoot, ...
    "SourceData_TableS4.xlsx");

requiredFiles = [ ...
    inputMat
    restorationInputMat
    mainSourcePath];

missingFiles = requiredFiles(~isfile(requiredFiles));

if ~isempty(missingFiles)
    error( ...
        "Required Table S4 input(s) not found:\n%s", ...
        strjoin(missingFiles, newline));
end


MainFit = normalizeFitSummary(readTablePortable(mainSourcePath, "fit_summary"));
MainCoef = normalizeCoefficientSummary(readTablePortable(mainSourcePath, "ratio_summary"));

[AltFit, AltCoef] = prepareAlternativeSpecificationInputs( ...
    inputMat, ...
    restorationInputMat);
AltFit = normalizeFitSummary(AltFit);
AltCoef = normalizeCoefficientSummary(AltCoef);

PointStats = buildPointThresholdStats(inputMat);

groupOrder = [ ...
    "Pooled", ...
    "Electric power", ...
    "Water supply", ...
    "Natural gas"];

displayGroupOrder = [ ...
    "All", ...
    "Electric power", ...
    "Water supply", ...
    "Natural gas"];


DisplayTable = table();

[rows, ~] = rowsFromFitCoefficientPair( ...
    MainFit, ...
    MainCoef, ...
    "Baseline", ...
    "Lrest=b*(Dmax*tau50) (>=6 points)", ...
    "baseline_Lrest_Dmax_tau50_ge6", ...
    6, ...
    groupOrder, ...
    displayGroupOrder);

DisplayTable = [DisplayTable; rows];

for threshold = [5 7 8]

    [rows, ~] = rowsFromPointSource( ...
        PointStats, ...
        "Trajectory resolution", ...
        sprintf("Lrest=b*(Dmax*tau50) (>=%d points)", threshold), ...
        sprintf("point_threshold_ge%d", threshold), ...
        threshold, ...
        groupOrder, ...
        displayGroupOrder);

    DisplayTable = [DisplayTable; rows];
end

modelSpecs = [
    "Lrest_Dmax_tau90", "Lrest=b*(Dmax*tau90) (>=6 points)", "model_Lrest_Dmax_tau90_ge6"
    "Lrest_Dmax_tau95", "Lrest=b*(Dmax*tau95) (>=6 points)", "model_Lrest_Dmax_tau95_ge6"
    "Lrest_D0_tau50",   "Lrest=b*(D0*tau50) (>=6 points)",   "model_Lrest_D0_tau50_ge6"
    "RL_Dmax_tau50",    "RL=b*(Dmax*tau50) (>=6 points)",    "model_RL_Dmax_tau50_ge6"];

for i = 1:size(modelSpecs, 1)

    analysisId = modelSpecs(i, 1);
    specification = modelSpecs(i, 2);
    sourceKey = modelSpecs(i, 3);

    fitRows = AltFit(AltFit.analysis == analysisId, :);
    coefRows = AltCoef(AltCoef.analysis == analysisId, :);

    [rows, ~] = rowsFromFitCoefficientPair( ...
        fitRows, ...
        coefRows, ...
        "Model specification", ...
        specification, ...
        sourceKey, ...
        6, ...
        groupOrder, ...
        displayGroupOrder);

    DisplayTable = [DisplayTable; rows];
end


if isfile(outXlsx)
    delete(outXlsx);
end

writeTableCompatPortable(DisplayTable, outXlsx, "display_table");

fprintf( ...
    "Saved Supplementary Table S4 source data:\n%s\n", ...
    outXlsx);


end

function [DisplayRows, SourceRows] = rowsFromFitCoefficientPair( ...
        Fit, ...
        Coef, ...
        sensitivityDimension, ...
        specification, ...
        sourceKey, ...
        minPoints, ...
        groupOrder, ...
        displayGroupOrder)

DisplayRows = table();
SourceRows = table();

for g = 1:numel(groupOrder)

    groupName = groupOrder(g);
    displayGroup = displayGroupOrder(g);

    fitRow = selectOneGroupRow(Fit, groupName, "fit summary", sourceKey);
    coefRow = selectOneGroupRow(Coef, groupName, "coefficient summary", sourceKey);

    [displayRow, sourceRow] = makeRows( ...
        sensitivityDimension, ...
        specification, ...
        sourceKey, ...
        minPoints, ...
        displayGroup, ...
        double(fitRow.n), ...
        double(fitRow.R2), ...
        double(fitRow.slope_b), ...
        double(fitRow.bootstrap_slope_ci95_low), ...
        double(fitRow.bootstrap_slope_ci95_high), ...
        double(coefRow.median), ...
        double(coefRow.Q25), ...
        double(coefRow.Q75), ...
        double(coefRow.Q05), ...
        double(coefRow.Q95));

    DisplayRows = [DisplayRows; displayRow];
    SourceRows = [SourceRows; sourceRow];
end

end

function [DisplayRows, SourceRows] = rowsFromPointSource( ...
        PointStats, ...
        sensitivityDimension, ...
        specification, ...
        sourceKey, ...
        minPoints, ...
        groupOrder, ...
        displayGroupOrder)

DisplayRows = table();
SourceRows = table();

for g = 1:numel(groupOrder)

    groupName = groupOrder(g);
    displayGroup = displayGroupOrder(g);

    source = PointStats( ...
        PointStats.min_points == minPoints & ...
        PointStats.group_normalized == groupName, :);

    if height(source) ~= 1
        error( ...
            "%s / %s: expected exactly one point-threshold row, found %d.", ...
            sourceKey, ...
            groupName, ...
            height(source));
    end

    [displayRow, sourceRow] = makeRows( ...
        sensitivityDimension, ...
        specification, ...
        sourceKey, ...
        minPoints, ...
        displayGroup, ...
        double(source.n), ...
        double(source.R2), ...
        double(source.slope_b), ...
        double(source.bootstrap_slope_ci95_low), ...
        double(source.bootstrap_slope_ci95_high), ...
        double(source.trajectory_coefficient_median), ...
        double(source.trajectory_coefficient_Q25), ...
        double(source.trajectory_coefficient_Q75), ...
        double(source.trajectory_coefficient_Q05), ...
        double(source.trajectory_coefficient_Q95));

    DisplayRows = [DisplayRows; displayRow];
    SourceRows = [SourceRows; sourceRow];
end

end

function SourceStats = buildPointThresholdStats(inputMat)

S = load(inputMat, "LossTailRecords");

if ~isfield(S, "LossTailRecords")
    error("MAT file does not contain LossTailRecords:\n%s", inputMat);
end

T = normalizeLossTailTable(S.LossTailRecords);

requiredFields = [ ...
    "record_ID", ...
    "system", ...
    "nRestorationRatioPoints", ...
    "nUniqueRestorationValues", ...
    "Dmax", ...
    "tau50", ...
    "Dmax_tau50", ...
    "Lrest", ...
    "Lrest_over_Dmax_tau50"];

assertRequiredFields(T, requiredFields, "LossTailRecords");

thresholds = [5 7 8];
groupOrder = ["Pooled", "Power", "Water", "Gas"];
nBoot = 2000;

validBase = ...
    positiveFinite(T.Dmax) & ...
    positiveFinite(T.tau50) & ...
    positiveFinite(T.Dmax_tau50) & ...
    positiveFinite(T.Lrest) & ...
    positiveFinite(T.Lrest_over_Dmax_tau50);

pointCounts = double(T.nRestorationRatioPoints);
uniqueValueCounts = double(T.nUniqueRestorationValues);
SourceStats = table();

for ti = 1:numel(thresholds)

    threshold = thresholds(ti);

    thresholdMask = ...
        isfinite(pointCounts) & ...
        pointCounts >= threshold & ...
        isfinite(uniqueValueCounts) & ...
        uniqueValueCounts >= 3;

    for gi = 1:numel(groupOrder)

        group = groupOrder(gi);
        mask = validBase & thresholdMask;

        if group ~= "Pooled"
            mask = mask & T.system == group;
        end

        x = double(T.Dmax_tau50(mask));
        y = double(T.Lrest(mask));
        coefficient = double(T.Lrest_over_Dmax_tau50(mask));

        stats = fitThroughOriginWithBootstrap(x, y, nBoot);
        q = quantileLocal(coefficient, [0.05 0.25 0.50 0.75 0.95]);

        SourceStats = [ ...
            SourceStats
            table( ...
                threshold, ...
                group, ...
                groupLabel(group), ...
                nnz(mask), ...
                stats.slope, ...
                stats.originalFitCiLow, ...
                stats.originalFitCiHigh, ...
                stats.bootstrapCiLow, ...
                stats.bootstrapCiHigh, ...
                stats.R2, ...
                q(1), ...
                q(2), ...
                q(3), ...
                q(4), ...
                q(5), ...
                nBoot, ...
                1, ...
                'VariableNames', { ...
                    'min_points', ...
                    'group', ...
                    'Lifeline', ...
                    'n', ...
                    'slope_b', ...
                    'original_fit_slope_ci95_low', ...
                    'original_fit_slope_ci95_high', ...
                    'bootstrap_slope_ci95_low', ...
                    'bootstrap_slope_ci95_high', ...
                    'R2', ...
                    'trajectory_coefficient_Q05', ...
                    'trajectory_coefficient_Q25', ...
                    'trajectory_coefficient_median', ...
                    'trajectory_coefficient_Q75', ...
                    'trajectory_coefficient_Q95', ...
                    'bootstrap_requested_B', ...
                    'bootstrap_seed'})];
    end
end

SourceStats = normalizePointSource(SourceStats);

end

function stats = fitThroughOriginWithBootstrap(x, y, nBoot)

x = double(x(:));
y = double(y(:));

keep = isfinite(x) & isfinite(y) & x > 0 & y > 0;
x = x(keep);
y = y(keep);
n = numel(x);

stats = struct( ...
    "slope", NaN, ...
    "R2", NaN, ...
    "originalFitCiLow", NaN, ...
    "originalFitCiHigh", NaN, ...
    "bootstrapCiLow", NaN, ...
    "bootstrapCiHigh", NaN, ...
    "bootSlope", nan(nBoot, 1));

if n < 2
    return;
end

sumX2 = sum(x .^ 2);

if sumX2 <= 0
    return;
end

stats.slope = sum(x .* y) / sumX2;
prediction = stats.slope .* x;
residual = y - prediction;
SSE = sum(residual .^ 2);
SST = sum((y - mean(y)) .^ 2);

if SST > 0
    stats.R2 = 1 - SSE / SST;
end

degreesFreedom = n - 1;

if degreesFreedom > 0
    residualVariance = SSE / degreesFreedom;
    slopeSE = sqrt(residualVariance / sumX2);
    tCritical = studentTCritical95(degreesFreedom);

    if isfinite(slopeSE) && isfinite(tCritical)
        stats.originalFitCiLow = stats.slope - tCritical * slopeSE;
        stats.originalFitCiHigh = stats.slope + tCritical * slopeSE;
    end
end

for b = 1:nBoot
    sampleIndex = randi(n, n, 1);
    xb = x(sampleIndex);
    yb = y(sampleIndex);
    denominator = sum(xb .^ 2);

    if denominator > 0
        stats.bootSlope(b) = sum(xb .* yb) / denominator;
    end
end

q = quantileLocal(stats.bootSlope, [0.025 0.975]);
stats.bootstrapCiLow = q(1);
stats.bootstrapCiHigh = q(2);

end

function tCritical = studentTCritical95(degreesFreedom)

tCritical = NaN;

if ~(isfinite(degreesFreedom) && degreesFreedom > 0)
    return;
end

if exist("tinv", "file") == 2
    tCritical = tinv(0.975, degreesFreedom);
elseif exist("betaincinv", "file") == 2
    betaQuantile = betaincinv(0.05, degreesFreedom / 2, 0.5);
    if isfinite(betaQuantile) && betaQuantile > 0
        tCritical = sqrt(degreesFreedom * (1 / betaQuantile - 1));
    end
else
    tCritical = 1.959963984540054;
end

end

function q = quantileLocal(values, probabilities)

values = sort(double(values(:)));
values = values(isfinite(values));
q = nan(size(probabilities));

if isempty(values)
    return;
end

for i = 1:numel(probabilities)

    p = probabilities(i);

    if numel(values) == 1
        q(i) = values(1);
    else
        position = 1 + (numel(values) - 1) * p;
        lowerIndex = floor(position);
        upperIndex = ceil(position);

        if lowerIndex == upperIndex
            q(i) = values(lowerIndex);
        else
            q(i) = values(lowerIndex) + ...
                (values(upperIndex) - values(lowerIndex)) * ...
                (position - lowerIndex);
        end
    end
end

end

function tf = positiveFinite(x)

x = double(x);
tf = isfinite(x) & x > 0;

end

function T = normalizeLossTailTable(T)

if ~istable(T)
    error("LossTailRecords must be a table.");
end

vars = string(T.Properties.VariableNames);

if ismember("System", vars)
    T.system = normalizeSystemName(T.System);
elseif ismember("system", vars)
    T.system = normalizeSystemName(T.system);
else
    error("LossTailRecords must contain System or system.");
end

if ismember("record_ID", string(T.Properties.VariableNames))
    T.record_ID = string(T.record_ID);
end

end

function assertRequiredFields(T, requiredFields, label)

available = string(T.Properties.VariableNames);
missing = requiredFields(~ismember(requiredFields, available));

if ~isempty(missing)
    error("%s is missing required field(s): %s", ...
        label, strjoin(missing, ", "));
end

end

function s = normalizeSystemName(s)

s = string(s);
sl = lower(strtrim(s));

s(sl == "e" | sl == "power") = "Power";
s(sl == "w" | sl == "water" | sl == "water supply") = "Water";
s(sl == "g" | sl == "gas" | sl == "natural gas") = "Gas";

end

function g = groupLabel(g)

if g == "Pooled"
    g = "Pooled";
elseif g == "Power"
    g = "Electric power";
elseif g == "Water"
    g = "Water supply";
elseif g == "Gas"
    g = "Natural gas";
end

end

function [displayRow, sourceRow] = makeRows( ...
        sensitivityDimension, ...
        specification, ...
        sourceKey, ...
        minPoints, ...
        lifelineSystem, ...
        n, ...
        R2, ...
        slopeB, ...
        bootLow, ...
        bootHigh, ...
        kappaMedian, ...
        kappaQ25, ...
        kappaQ75, ...
        kappaQ05, ...
        kappaQ95)

displayRow = table( ...
    string(sensitivityDimension), ...
    string(specification), ...
    string(lifelineSystem), ...
    formatInteger(n), ...
    formatR2(R2), ...
    formatNumber(slopeB), ...
    formatToInterval(bootLow, bootHigh), ...
    formatNumber(kappaMedian), ...
    formatInterval(kappaQ25, kappaQ75), ...
    formatInterval(kappaQ05, kappaQ95), ...
    'VariableNames', { ...
        'Sensitivity_dimension', ...
        'Specification', ...
        'Lifeline_system', ...
        'Sample_size', ...
        'R2', ...
        'Fitted_slope_b', ...
        'Bootstrap_95CI_for_b', ...
        'Loss_shape_factor_median', ...
        'Loss_shape_factor_Q25_Q75', ...
        'Loss_shape_factor_Q05_Q95'});

sourceRow = table( ...
    string(sensitivityDimension), ...
    string(specification), ...
    string(sourceKey), ...
    minPoints, ...
    string(lifelineSystem), ...
    n, ...
    R2, ...
    slopeB, ...
    bootLow, ...
    bootHigh, ...
    kappaMedian, ...
    kappaQ25, ...
    kappaQ75, ...
    kappaQ05, ...
    kappaQ95, ...
    'VariableNames', { ...
        'Sensitivity_dimension', ...
        'Specification', ...
        'source_key', ...
        'min_points', ...
        'Lifeline_system', ...
        'n', ...
        'R2', ...
        'slope_b', ...
        'bootstrap_slope_ci95_low', ...
        'bootstrap_slope_ci95_high', ...
        'loss_shape_factor_median', ...
        'loss_shape_factor_Q25', ...
        'loss_shape_factor_Q75', ...
        'loss_shape_factor_Q05', ...
        'loss_shape_factor_Q95'});

end


function T = readTablePortable(path, sheetName)

if nargin < 2
    T = readtable( ...
        path, ...
        "TextType", "string", ...
        "VariableNamingRule", "preserve");
else
    T = readtable( ...
        path, ...
        "Sheet", sheetName, ...
        "TextType", "string", ...
        "VariableNamingRule", "preserve");
end

end

function T = normalizeFitSummary(T)

T = normalizeStrings(T);
T.group_normalized = normalizeGroupName(groupColumn(T));

if ismember("analysis", string(T.Properties.VariableNames))
    T.analysis = string(T.analysis);
else
    T.analysis = repmat("", height(T), 1);
end

T.n = numericColumn(T, "n");
T.R2 = numericColumn(T, "R2");
T.slope_b = numericColumn(T, "slope_b");
T.bootstrap_slope_ci95_low = numericColumn(T, "bootstrap_slope_ci95_low");
T.bootstrap_slope_ci95_high = numericColumn(T, "bootstrap_slope_ci95_high");

end

function T = normalizeCoefficientSummary(T)

T = normalizeStrings(T);
T.group_normalized = normalizeGroupName(groupColumn(T));

if ismember("analysis", string(T.Properties.VariableNames))
    T.analysis = string(T.analysis);
else
    T.analysis = repmat("", height(T), 1);
end

T.n = optionalNumericColumn(T, "n");
T.median = numericColumn(T, "median");
T.Q25 = numericColumn(T, "Q25");
T.Q75 = numericColumn(T, "Q75");
T.Q05 = numericColumn(T, "Q05");
T.Q95 = numericColumn(T, "Q95");

end

function T = normalizePointSource(T)

T = normalizeStrings(T);
T.group_normalized = normalizeGroupName(groupColumn(T));

T.min_points = numericColumn(T, "min_points");
T.n = numericColumn(T, "n");
T.R2 = numericColumn(T, "R2");
T.slope_b = numericColumn(T, "slope_b");
T.bootstrap_slope_ci95_low = numericColumn(T, "bootstrap_slope_ci95_low");
T.bootstrap_slope_ci95_high = numericColumn(T, "bootstrap_slope_ci95_high");
T.trajectory_coefficient_median = numericColumn(T, "trajectory_coefficient_median");
T.trajectory_coefficient_Q25 = numericColumn(T, "trajectory_coefficient_Q25");
T.trajectory_coefficient_Q75 = numericColumn(T, "trajectory_coefficient_Q75");
T.trajectory_coefficient_Q05 = numericColumn(T, "trajectory_coefficient_Q05");
T.trajectory_coefficient_Q95 = numericColumn(T, "trajectory_coefficient_Q95");

end

function row = selectOneGroupRow(T, groupName, tableLabel, sourceKey)

idx = T.group_normalized == groupName;

if nnz(idx) ~= 1
    error( ...
        "%s / %s / %s: expected exactly one row, found %d.", ...
        sourceKey, ...
        groupName, ...
        tableLabel, ...
        nnz(idx));
end

row = T(idx, :);

end

function groupRaw = groupColumn(T)

vars = string(T.Properties.VariableNames);

if ismember("system", vars)
    groupRaw = string(T.system);
elseif ismember("Lifeline", vars)
    groupRaw = string(T.Lifeline);
elseif ismember("lifeline", vars)
    groupRaw = string(T.lifeline);
elseif ismember("group", vars)
    groupRaw = string(T.group);
else
    error("Table has no group/system/lifeline column.");
end

end

function g = normalizeGroupName(x)

x = lower(strtrim(string(x)));
g = strings(size(x));

for i = 1:numel(x)

    switch x(i)

        case {"pooled", "all", "all systems"}
            g(i) = "Pooled";

        case {"e", "power", "electric power"}
            g(i) = "Electric power";

        case {"w", "water", "water supply"}
            g(i) = "Water supply";

        case {"g", "gas", "natural gas"}
            g(i) = "Natural gas";

        otherwise
            g(i) = strtrim(string(x(i)));
    end
end

end

function x = numericColumn(T, name)

name = string(name);
vars = string(T.Properties.VariableNames);

if ~ismember(name, vars)
    error("Required column '%s' is missing.", name);
end

x = toDouble(T.(char(name)));

end

function x = optionalNumericColumn(T, name)

name = string(name);
vars = string(T.Properties.VariableNames);

if ~ismember(name, vars)
    x = nan(height(T), 1);
else
    x = toDouble(T.(char(name)));
end

end

function x = toDouble(v)

if isnumeric(v) || islogical(v)
    x = double(v);
else
    x = str2double(string(v));
end

end

function T = normalizeStrings(T)

vars = T.Properties.VariableNames;

for i = 1:numel(vars)

    v = T.(vars{i});

    if isstring(v)
        v(ismissing(v)) = "";
        T.(vars{i}) = v;
    elseif iscellstr(v) || ischar(v)
        T.(vars{i}) = string(v);
    end
end

end


function s = formatNumber(x)

if ~isfinite(x)
    s = "";
else
    s = string(sprintf("%.3f", x));
end

end

function s = formatR2(x)

if ~isfinite(x)
    s = "";
else
    s = string(sprintf("%.3f", x));
end

end

function s = formatInteger(x)

if ~isfinite(x)
    s = "";
else
    s = string(sprintf("%d", round(x)));
end

end

function s = formatToInterval(low, high)

if ~(isfinite(low) && isfinite(high))
    s = "";
else
    s = formatNumber(low) + " to " + formatNumber(high);
end

end

function s = formatInterval(low, high)

if ~(isfinite(low) && isfinite(high))
    s = "";
else
    s = "[" + formatNumber(low) + "-" + formatNumber(high) + "]";
end

end

function ensureDir(directoryPath)

if ~isfolder(directoryPath)
    mkdir(directoryPath);
end

end

function [AllFitSummary, AllCoefficientSummary] = prepareAlternativeSpecificationInputs( ...
        inputMat, ...
        restorationInputMat)

if ~isfile(inputMat)
    error( ...
        "Prepared loss-tail MAT not found:\n%s\n" + ...
        "Run run_00_prepare_loss_tail_dataset.m first.", ...
        inputMat);
end

S = load(inputMat, "LossTailRecords");

if ~isfield(S, "LossTailRecords")
    error( ...
        "MAT file does not contain LossTailRecords:\n%s", ...
        inputMat);
end

T = normalizeLossTailTable(S.LossTailRecords);
T = addRestorationLossAlternative(T, restorationInputMat);

requiredFields = [ ...
    "record_ID", ...
    "system", ...
    "nRestorationRatioPoints", ...
    "Dmax", ...
    "D0", ...
    "tau50", ...
    "tau90", ...
    "tau95", ...
    "Lrest", ...
    "RL", ...
    "Dmax_tau50", ...
    "RL_over_Dmax_tau50", ...
    "Lrest_over_D0_tau50", ...
    "Lrest_over_Dmax_tau90", ...
    "Lrest_over_Dmax_tau95"];

assertRequiredFields(T, requiredFields, "LossTailRecords");

baselineMinPoints = 6;
nBoot = 2000;
sysOrder = ["Power", "Water", "Gas"];

pointCounts = double(T.nRestorationRatioPoints);

T = T(isfinite(pointCounts) & pointCounts >= baselineMinPoints, :);

if isempty(T)
    error( ...
        "No LossTailRecords satisfy the >=%d-point baseline.", ...
        baselineMinPoints);
end

jobs(1,1) = makeAlternativeJob( ...
    "RL_Dmax_tau50", ...
    "RL", ...
    "RL_over_Dmax_tau50", ...
    "Dmax", ...
    "tau50", ...
    "Dmax_tau50");

jobs = repmat(jobs(1), 4, 1);

jobs(2) = makeAlternativeJob( ...
    "Lrest_D0_tau50", ...
    "Lrest", ...
    "Lrest_over_D0_tau50", ...
    "D0", ...
    "tau50", ...
    "");

jobs(3) = makeAlternativeJob( ...
    "Lrest_Dmax_tau90", ...
    "Lrest", ...
    "Lrest_over_Dmax_tau90", ...
    "Dmax", ...
    "tau90", ...
    "");

jobs(4) = makeAlternativeJob( ...
    "Lrest_Dmax_tau95", ...
    "Lrest", ...
    "Lrest_over_Dmax_tau95", ...
    "Dmax", ...
    "tau95", ...
    "");

AllFitSummary = table();
AllCoefficientSummary = table();

for j = 1:numel(jobs)

    job = jobs(j);
    D = buildAlternativeCompressionData(T, sysOrder, job);

    [FitSummary, CoefficientSummary] = ...
        summarizeAlternativeCompression(D, sysOrder, nBoot);

    AllFitSummary = [ ...
        AllFitSummary
        addAlternativeAnalysisMetadata(FitSummary, job)];

    AllCoefficientSummary = [ ...
        AllCoefficientSummary
        addAlternativeCoefficientMetadata(CoefficientSummary, job)];
end

end

function T = addRestorationLossAlternative(T, restorationInputMat)

S = load(restorationInputMat, "unique3Records");
if ~isfield(S, "unique3Records")
    error("Input MAT must contain unique3Records:\n%s", restorationInputMat);
end

records = S.unique3Records;
recordIDs = string({records.record_ID});

RL = NaN(height(T), 1);

for i = 1:height(T)
    idx = find(recordIDs == string(T.record_ID(i)), 1);
    if isempty(idx)
        error("Record not found in unique3Records: %s", string(T.record_ID(i)));
    end

    rec = records(idx);
    RL(i) = computeRLFromFunctionality( ...
        rec.functionalityCurve, ...
        getNumericField(rec, "PeakOutageFraction"), ...
        getNumericField(rec, "PeakTime"), ...
        getNumericField(rec, "T100"));
end

T.RL = RL;
T.RL_over_Dmax_tau50 = safeRatioColumn(RL, double(T.Dmax_tau50));

end

function RL = computeRLFromFunctionality(C, Dmax, peakTime, T100)

recoverTol = 1e-8;
RL = NaN;

if ~(isfinite(Dmax) && Dmax >= 0 && ...
        isfinite(peakTime) && peakTime >= 0 && ...
        isfinite(T100) && T100 >= peakTime - recoverTol)
    return;
end

[t, functionality] = functionalityCurveToVectors(C);
if isempty(t)
    return;
end

loss = 1 - functionality;
requireUnitInterval(loss, "functionality-derived loss");

tAfter = t(t > peakTime + recoverTol);
lossAfter = loss(t > peakTime + recoverTol);

if peakTime > recoverTol
    t = [0; peakTime; tAfter];
    loss = [Dmax; Dmax; lossAfter];
else
    t = [0; tAfter];
    loss = [Dmax; lossAfter];
end

[tUse, lossUse, ok] = curveOnWindow(t, loss, 0, T100, recoverTol);
if ok
    RL = trapz(tUse, lossUse);
end

end

function [t, f] = functionalityCurveToVectors(C)

t = [];
f = [];

if isempty(C) || ~isstruct(C) || ...
        ~isfield(C, "day") || ~isfield(C, "value")
    return;
end

t = double([C.day]');
f = double([C.value]');

mask = isfinite(t) & isfinite(f);
t = t(mask);
f = f(mask);

if isempty(t)
    return;
end

[t, idx] = sort(t);
f = f(idx);

if any(diff(t) == 0)
    error("Functionality curve contains duplicate day values.");
end

end

function [tUse, yUse, ok] = curveOnWindow(t, y, t0, t1, tol)

ok = false;
tUse = [];
yUse = [];

if isempty(t) || isempty(y) || t0 < min(t) - tol || t1 > max(t) + tol
    return;
end

[y0, ok0] = valueAtOrInterpolate(t, y, t0, tol);
[y1, ok1] = valueAtOrInterpolate(t, y, t1, tol);

if ~(ok0 && ok1)
    return;
end

mask = t > t0 + tol & t < t1 - tol;
tUse = [t0; t(mask); t1];
yUse = [y0; y(mask); y1];
ok = numel(tUse) >= 2;

end

function [vq, ok] = valueAtOrInterpolate(t, y, tq, tol)

vq = NaN;
ok = false;

idx = find(abs(t - tq) <= tol, 1, "first");
if ~isempty(idx)
    vq = y(idx);
    ok = true;
    return;
end

left = find(t < tq, 1, "last");
right = find(t > tq, 1, "first");

if isempty(left) || isempty(right) || t(right) == t(left)
    return;
end

vq = y(left) + (tq - t(left)) * (y(right) - y(left)) / ...
    (t(right) - t(left));
ok = isfinite(vq);

end

function requireUnitInterval(x, label)

tol = 1e-8;
x = double(x(:));
x = x(isfinite(x));

if any(x < -tol | x > 1 + tol)
    error("%s contains value(s) outside [0,1].", label);
end

end

function z = safeRatioColumn(a, b)

z = NaN(size(a));
mask = isfinite(a) & isfinite(b) & b > 0;
z(mask) = a(mask) ./ b(mask);

end

function v = getNumericField(s, name)

v = NaN;

if ~isfield(s, name)
    return;
end

x = s.(name);
if isnumeric(x) || islogical(x)
    if ~isempty(x)
        v = double(x(1));
    end
else
    v = str2double(string(x));
end

end

function job = makeAlternativeJob( ...
        id, ...
        yField, ...
        coefficientField, ...
        scaleField, ...
        tauField, ...
        precomputedXField)

job = struct();
job.id = string(id);
job.yField = string(yField);
job.coefficientField = string(coefficientField);
job.scaleField = string(scaleField);
job.tauField = string(tauField);
job.precomputedXField = string(precomputedXField);
job.equation = ...
    string(yField) + ...
    " = b*(" + ...
    string(scaleField) + ...
    "*" + ...
    string(tauField) + ...
    ")";
job.coefficientDefinition = ...
    "c_i = " + ...
    string(coefficientField);

end

function D = buildAlternativeCompressionData(T, sysOrder, job)

system = normalizeSystemName(T.system);

y = double(T.(char(job.yField)));
scale = double(T.(char(job.scaleField)));
tau = double(T.(char(job.tauField)));
coefficient = double(T.(char(job.coefficientField)));

if strlength(job.precomputedXField) > 0
    x = double(T.(char(job.precomputedXField)));
else
    x = scale .* tau;
end

keep = ...
    ismember(system, sysOrder) & ...
    isfinite(y) & ...
    isfinite(scale) & ...
    isfinite(tau) & ...
    isfinite(x) & ...
    isfinite(coefficient) & ...
    y > 0 & ...
    scale > 0 & ...
    tau > 0 & ...
    x > 0 & ...
    coefficient > 0;

D = table( ...
    string(T.record_ID(keep)), ...
    system(keep), ...
    scale(keep), ...
    tau(keep), ...
    x(keep), ...
    y(keep), ...
    coefficient(keep), ...
    'VariableNames', { ...
        'record_ID', ...
        'system', ...
        'scale', ...
        'tau', ...
        'x', ...
        'y', ...
        'coefficient_ci'});

end

function [FitSummary, CoefficientSummary] = ...
        summarizeAlternativeCompression(D, sysOrder, nBoot)

groups = ["Pooled", sysOrder];

FitSummary = table();
CoefficientSummary = table();

for g = 1:numel(groups)

    if groups(g) == "Pooled"
        idx = true(height(D), 1);
    else
        idx = D.system == groups(g);
    end

    label = groupLabel(groups(g));
    stats = fitThroughOriginWithBootstrap(D.x(idx), D.y(idx), nBoot);

    FitSummary = [ ...
        FitSummary
        table( ...
            label, ...
            nnz(idx), ...
            stats.slope, ...
            stats.originalFitCiLow, ...
            stats.originalFitCiHigh, ...
            stats.bootstrapCiLow, ...
            stats.bootstrapCiHigh, ...
            stats.R2, ...
            'VariableNames', { ...
                'system', ...
                'n', ...
                'slope_b', ...
                'original_fit_slope_ci95_low', ...
                'original_fit_slope_ci95_high', ...
                'bootstrap_slope_ci95_low', ...
                'bootstrap_slope_ci95_high', ...
                'R2'})];

    q = quantileLocal( ...
        D.coefficient_ci(idx), ...
        [0.05 0.25 0.50 0.75 0.95]);

    CoefficientSummary = [ ...
        CoefficientSummary
        table( ...
            label, ...
            nnz(idx), ...
            q(3), ...
            q(2), ...
            q(4), ...
            q(1), ...
            q(5), ...
            'VariableNames', { ...
                'system', ...
                'n', ...
                'median', ...
                'Q25', ...
                'Q75', ...
                'Q05', ...
                'Q95'})];
end

end

function Out = addAlternativeAnalysisMetadata(T, job)

nRows = height(T);

Meta = table( ...
    repmat(string(job.id), nRows, 1), ...
    repmat(string(job.equation), nRows, 1), ...
    'VariableNames', { ...
        'analysis', ...
        'equation'});

Out = [Meta, T];

end

function Out = addAlternativeCoefficientMetadata(T, job)

nRows = height(T);

Meta = table( ...
    repmat(string(job.id), nRows, 1), ...
    repmat(string(job.coefficientDefinition), nRows, 1), ...
    'VariableNames', { ...
        'analysis', ...
        'coefficient_definition'});

Out = [Meta, T];

end
