function plot_FigS12(supplementaryRoot, figureDataRoot)

if nargin < 1 || isempty(supplementaryRoot)
    moduleRoot = string(fileparts(mfilename("fullpath")));
    emergentSimplicityRoot = string(fileparts(fileparts(moduleRoot)));
    supplementaryRoot = fullfile(emergentSimplicityRoot, "outputs", "supplementary");
end
if nargin >= 2 && strlength(string(figureDataRoot)) > 0
    renderFigS12FromSourceData(supplementaryRoot, figureDataRoot);
    return;
end

prepare_FigS12_contextual_data(supplementaryRoot);
fit_FigS12_contextual_stability(supplementaryRoot);
draw_FigS12_contextual_stability(supplementaryRoot);

end

function renderFigS12FromSourceData(figureDir, figureDataRoot)

sourceXlsx = fullfile(figureDataRoot, "SourceData_FigS12.xlsx");
if ~isfile(sourceXlsx)
    error("Frozen Figure S12 source workbook not found: %s", sourceXlsx);
end
ensureDirPlot(figureDir);

M = readtable(sourceXlsx, "Sheet", "model_summary", ...
    "TextType", "string", "VariableNamingRule", "preserve");
E = readtable(sourceXlsx, "Sheet", "contrast_plot_source", ...
    "TextType", "string", "VariableNamingRule", "preserve");
if ~ismember("significance", string(E.Properties.VariableNames))
    if ~ismember("p_value", string(E.Properties.VariableNames))
        error("Frozen Figure S12 contrast data require p_value or significance.");
    end
    E.significance = significanceLabels(double(E.p_value));
end
plotContextualPanels(E, M, figureDir);
fprintf("Saved Fig S12 from frozen figure data under:\n%s\n", figureDir);

end

function labels = significanceLabels(p)

p = double(p(:));
labels = strings(size(p));
labels(isfinite(p) & p < 0.001) = "***";
labels(isfinite(p) & p >= 0.001 & p < 0.01) = "**";
labels(isfinite(p) & p >= 0.01 & p < 0.05) = "*";

end

function prepare_FigS12_contextual_data(supplementaryRoot)

if nargin < 1 || isempty(supplementaryRoot)
    moduleRoot = string(fileparts(mfilename("fullpath")));
    emergentSimplicityRoot = string(fileparts(fileparts(moduleRoot)));
    supplementaryRoot = fullfile(emergentSimplicityRoot, "outputs", "supplementary");
end

moduleRoot = string(fileparts(mfilename("fullpath")));
emergentSimplicityRoot = string(fileparts(fileparts(moduleRoot)));
analysisRoot = string(fileparts(emergentSimplicityRoot));
intermediateRoot = fullfile(emergentSimplicityRoot, "outputs", "intermediate");

outputRoot = fullfile(supplementaryRoot, "source_data");
ensureDirPrepare(outputRoot);

standardizedMat = fullfile( ...
    analysisRoot, ...
    "00_Dataset", ...
    "standardized_dataset.mat");

shapeMat = fullfile( ...
    intermediateRoot, ...
    "01_restoration_shape_phase_ratio", ...
    "unique3_records.mat");

lossMat = fullfile( ...
    intermediateRoot, ...
    "02_loss_compression_late_loss", ...
    "loss_tail_dataset.mat");

pairMat = fullfile( ...
    intermediateRoot, ...
    "03_cross_lifeline_timing_coupling", ...
    "paired_timing_coupling_dataset.mat");

requiredFiles = [standardizedMat; shapeMat; lossMat; pairMat];

if any(~isfile(requiredFiles))

    error( ...
        "Required input file(s) are missing:\n%s", ...
        strjoin(requiredFiles(~isfile(requiredFiles)), newline));
end

Sstd = load(standardizedMat, "event_region_system_records");
Sshape = load(shapeMat, "unique3Records");
Sloss = load(lossMat, "LossTailRecords");
Spair = load(pairMat, "Paired");

if ~isfield(Sstd, "event_region_system_records")
    error("standardized_dataset.mat is missing event_region_system_records.");
end

if ~isfield(Sshape, "unique3Records")
    error("unique3_records.mat is missing unique3Records.");
end

if ~isfield(Sloss, "LossTailRecords")
    error("loss_tail_dataset.mat is missing LossTailRecords.");
end

if ~isfield(Spair, "Paired")
    error("paired_timing_coupling_dataset.mat is missing Paired.");
end

std = standardizeRecordId(asTable(Sstd.event_region_system_records));
shape = standardizeRecordId(asTable(Sshape.unique3Records));
loss = standardizeRecordId(asTable(Sloss.LossTailRecords));
pairsRaw = standardizePairColumns(asTable(Spair.Paired));


requireColumns( ...
    std, ...
    ["record_ID","event_ID","System","avg_PGA","IncomeGroup","population","OccTime"], ...
    "event_region_system_records");

Context = std(:, [ ...
    "record_ID", ...
    "event_ID", ...
    "System", ...
    "avg_PGA", ...
    "IncomeGroup", ...
    "population", ...
    "OccTime"]);

Context.record_ID = string(Context.record_ID);
Context.System = normalizeSystem(Context.System);
Context.avg_PGA = numericColumnPrepare(Context.avg_PGA);
Context.IncomeGroup = strtrim(string(Context.IncomeGroup));
Context.population = numericColumnPrepare(Context.population);
Context.event_year = extractEventYear(Context.OccTime);

Context.OccTime = [];

Context = Context(:, [ ...
    "record_ID", ...
    "event_ID", ...
    "System", ...
    "avg_PGA", ...
    "IncomeGroup", ...
    "event_year", ...
    "population"]);

if numel(unique(Context.record_ID)) ~= height(Context)

    error( ...
        "Standardized record_ID is not unique; contextual matching is ambiguous.");
end


requireColumns( ...
    shape, ...
    ["record_ID","nRestorationRatioPoints","nUniqueRestorationValues", ...
     "tau80_over_tau50","tau90_over_tau80","tau95_over_tau90"], ...
    "unique3Records");

Phase = shape(:, [ ...
    "record_ID", ...
    "nRestorationRatioPoints", ...
    "nUniqueRestorationValues", ...
    "tau80_over_tau50", ...
    "tau90_over_tau80", ...
    "tau95_over_tau90"]);

Phase.record_ID = string(Phase.record_ID);
Phase.nRestorationRatioPoints = numericColumnPrepare(Phase.nRestorationRatioPoints);
Phase.nUniqueRestorationValues = numericColumnPrepare(Phase.nUniqueRestorationValues);
Phase.tau80_over_tau50 = numericColumnPrepare(Phase.tau80_over_tau50);
Phase.tau90_over_tau80 = numericColumnPrepare(Phase.tau90_over_tau80);
Phase.tau95_over_tau90 = numericColumnPrepare(Phase.tau95_over_tau90);

keep = ...
    Phase.nRestorationRatioPoints >= 6 & ...
    Phase.nUniqueRestorationValues >= 3;

Phase = Phase(keep, :);

Phase.nRestorationRatioPoints = [];
Phase.nUniqueRestorationValues = [];

Phase = innerjoin(Phase, Context, "Keys", "record_ID");

Phase = Phase(:, [ ...
    "record_ID", ...
    "event_ID", ...
    "System", ...
    "tau80_over_tau50", ...
    "tau90_over_tau80", ...
    "tau95_over_tau90", ...
    "avg_PGA", ...
    "IncomeGroup", ...
    "event_year", ...
    "population"]);


sourceCompressionField = "Lrest_over_Dmax_tau50";

requireColumns( ...
    loss, ...
    ["record_ID","nRestorationRatioPoints",sourceCompressionField,"eta90"], ...
    "LossTailRecords");

Loss = loss(:, [ ...
    "record_ID", ...
    "nRestorationRatioPoints", ...
    sourceCompressionField, ...
    "eta90"]);

Loss.record_ID = string(Loss.record_ID);
Loss.nRestorationRatioPoints = numericColumnPrepare(Loss.nRestorationRatioPoints);
Loss.(sourceCompressionField) = numericColumnPrepare(Loss.(sourceCompressionField));
Loss.eta90 = numericColumnPrepare(Loss.eta90);

Loss = Loss(Loss.nRestorationRatioPoints >= 6, :);

nameIndex = find( ...
    string(Loss.Properties.VariableNames) == sourceCompressionField, ...
    1, ...
    "first");

Loss.Properties.VariableNames(nameIndex) = {'c_i'};

Loss.nRestorationRatioPoints = [];

Loss = innerjoin(Loss, Context, "Keys", "record_ID");

Loss = Loss(:, [ ...
    "record_ID", ...
    "event_ID", ...
    "System", ...
    "c_i", ...
    "eta90", ...
    "avg_PGA", ...
    "IncomeGroup", ...
    "event_year", ...
    "population"]);


requireColumns( ...
    pairsRaw, ...
    ["pair_ID","systemPair","event_ID","record_ID_E","record_ID_B", ...
     "rho90_B_over_E","rho95_B_over_E","CmaxAB90"], ...
    "Paired");

Pair = pairsRaw(:, [ ...
    "pair_ID", ...
    "systemPair", ...
    "event_ID", ...
    "record_ID_E", ...
    "record_ID_B", ...
    "rho90_B_over_E", ...
    "rho95_B_over_E", ...
    "CmaxAB90"]);

Pair.record_ID_E = string(Pair.record_ID_E);
Pair.record_ID_B = string(Pair.record_ID_B);
Pair.systemPair = normalizePair(Pair.systemPair);
Pair.rho90_B_over_E = numericColumnPrepare(Pair.rho90_B_over_E);
Pair.rho95_B_over_E = numericColumnPrepare(Pair.rho95_B_over_E);
Pair.CmaxAB90 = numericColumnPrepare(Pair.CmaxAB90);

ContextE = Context(:, [ ...
    "record_ID", ...
    "avg_PGA", ...
    "IncomeGroup", ...
    "event_year", ...
    "population"]);

ContextE.Properties.VariableNames(1) = {'record_ID_E'};

Pair = innerjoin(Pair, ContextE, "Keys", "record_ID_E");

Pair.record_ID_E = [];
Pair.record_ID_B = [];

Pair = Pair(:, [ ...
    "pair_ID", ...
    "event_ID", ...
    "systemPair", ...
    "rho90_B_over_E", ...
    "rho95_B_over_E", ...
    "CmaxAB90", ...
    "avg_PGA", ...
    "IncomeGroup", ...
    "event_year", ...
    "population"]);


outMat = fullfile( ...
    outputRoot, ...
    "contextual_stability_data.mat");

save( ...
    outMat, ...
    "Phase", ...
    "Loss", ...
    "Pair");

oldXlsx = fullfile(outputRoot, "contextual_stability_data.xlsx");
if isfile(oldXlsx)
    delete(oldXlsx);
end

fprintf("\nPrepared contextual-stability datasets:\n");
fprintf("  Phase records : %d\n", height(Phase));
fprintf("  Loss records  : %d\n", height(Loss));
fprintf("  Pair records  : %d\n", height(Pair));
fprintf("Saved:\n  %s\n", outMat);


end

function T = asTable(x)

if istable(x)

    T = x;

elseif isstruct(x)

    T = struct2table( ...
        x, ...
        "AsArray", ...
        true);

else

    error("Expected a table or struct input.");
end

end

function T = standardizeRecordId(T)

vars = string(T.Properties.VariableNames);

if ismember("record_ID", vars)
    return;
end

aliases = ["recordID","RecordID","recordId"];

idx = find( ...
    ismember(aliases, vars), ...
    1, ...
    "first");

if isempty(idx)
    return;
end

oldName = aliases(idx);

nameIndex = find( ...
    vars == oldName, ...
    1, ...
    "first");

T.Properties.VariableNames(nameIndex) = {'record_ID'};

end

function T = standardizePairColumns(T)

aliases = { ...
    "pairID", "pair_ID"; ...
    "SystemPair", "systemPair"; ...
    "recordIDE", "record_ID_E"; ...
    "recordIDB", "record_ID_B"};

for i = 1:size(aliases, 1)

    oldName = string(aliases{i,1});
    newName = string(aliases{i,2});

    vars = string(T.Properties.VariableNames);

    if ~ismember(newName, vars) && ismember(oldName, vars)

        nameIndex = find( ...
            vars == oldName, ...
            1, ...
            "first");

        T.Properties.VariableNames(nameIndex) = {char(newName)};
    end
end

end

function requireColumns(T, required, label)

required = string(required);
vars = string(T.Properties.VariableNames);

missing = required(~ismember(required, vars));

if ~isempty(missing)

    error( ...
        "%s is missing required column(s): %s", ...
        label, ...
        strjoin(missing, ", "));
end

end

function x = numericColumnPrepare(v)

if isnumeric(v) || islogical(v)

    x = double(v);

else

    x = str2double(string(v));
end

end

function s = normalizeSystem(s)

s = lower(strtrim(string(s)));

s(ismember(s, ["e","electric power","power"])) = "Power";
s(ismember(s, ["w","water","water supply"])) = "Water";
s(ismember(s, ["g","gas","natural gas"])) = "Gas";

end

function s = normalizePair(s)

s = upper(strtrim(string(s)));
s = replace( ...
    s, ...
    [" ","/","-","_","\rightarrow"], ...
    "");

s(s == "EW" | s == "E>W") = "EW";
s(s == "EG" | s == "E>G") = "EG";

end

function y = extractEventYear(v)

n = numel(v);
y = nan(n,1);

if isdatetime(v)

    y = year(v);
    return;
end

for i = 1:n

    if iscell(v)
        value = v{i};
    else
        value = v(i);
    end

    if isdatetime(value)

        y(i) = year(value);
        continue;
    end

    if isnumeric(value) && isscalar(value) && isfinite(value)

        try
            y(i) = year(datetime(value, "ConvertFrom", "excel"));
        catch
        end

        continue;
    end

    txt = strtrim(string(value));

    if txt == "" || ismissing(txt)
        continue;
    end

    try

        y(i) = year(datetime(txt));

    catch

        tok = regexp( ...
            char(txt), ...
            '(18|19|20)\d{2}', ...
            'match', ...
            'once');

        if ~isempty(tok)
            y(i) = str2double(tok);
        end
    end
end

end

function ensureDirPrepare(d)

if exist(d, "dir") ~= 7
    mkdir(d);
end

end



function fit_FigS12_contextual_stability(supplementaryRoot)

if nargin < 1 || isempty(supplementaryRoot)
    moduleRoot = string(fileparts(mfilename("fullpath")));
    emergentSimplicityRoot = string(fileparts(fileparts(moduleRoot)));
    supplementaryRoot = fullfile(emergentSimplicityRoot, "outputs", "supplementary");
end

cfg.centerYear = 2000;
cfg.systemOrder = ["Power","Water","Gas"];
cfg.pairOrder = ["EW","EG"];

cfg.incomeOrder = ["Non-high income","High income"];
cfg.incomeReference = "Non-high income";

cfg.dfMethod = "satterthwaite";
cfg.bootstrapB = 2000;
cfg.bootstrapSeed = 20260903;
cfg.bootstrapMinSuccessFraction = 0.90;
cfg.bootstrapPercentiles = [2.5 97.5];

moduleRoot = string(fileparts(mfilename("fullpath")));

outputRoot = fullfile(supplementaryRoot, "source_data");

ensureDirFit(outputRoot);

inputMat = fullfile(outputRoot, "contextual_stability_data.mat");

if ~isfile(inputMat)
    error("Missing prepared contextual data:\n%s\nRun run_00 first.", inputMat);
end

S = load(inputMat, "Phase", "Loss", "Pair");
Phase = S.Phase;
Loss = S.Loss;
Pair = S.Pair;


indicators = buildIndicators(Phase, Loss, Pair, cfg);

ModelSummary = table();
ContextualContrasts = table();

for f = 1:numel(indicators)

    fam = indicators(f);
    T = fam.data;

    fprintf("\n============================================================\n");
    fprintf("Contextual indicator: %s\n", fam.name);
    fprintf("============================================================\n");

    if isempty(T)
        warning("No complete observations for indicator %s.", fam.name);
        continue;
    end

    nObs = height(T);
    nEvents = numel(unique(string(T.event_ID)));
    if fam.unitVar ~= ""
        nUnits = numel(unique(string(T.(fam.unitVar))));
    else
        nUnits = nObs;
    end

    fprintf("Observations: %d | units: %d | earthquakes: %d\n", ...
        nObs, nUnits, nEvents);

    if nEvents < 5
        warning("%s contains fewer than five earthquakes; inference is not reliable.", fam.name);
    end

    M0 = fitLmeChecked(T, fam.formula0);
    M1 = fitLmeChecked(T, fam.formula1);

    R2m0 = marginalR2(M0, T);
    R2m1 = marginalR2(M1, T);
    deltaR2m = R2m1 - R2m0;

    thisDiag = computeResidualDiagnostics(M1, T, fam);

    dist = computeEmpiricalDistribution(T.raw_response, T.response, fam);

    thisCoefficients = extractContextualCoefficients( ...
        M1, fam, nObs, nUnits, nEvents, cfg);

    thisContrasts = makeContextualContrasts( ...
        thisCoefficients, fam, dist, T);

    boot = bootstrapContextualContrasts(T, fam, thisContrasts, cfg);
    thisContrasts.CI95_boot_low = boot.CI95_boot_low;
    thisContrasts.CI95_boot_high = boot.CI95_boot_high;
    thisContrasts.analysis_shift_CI95_boot_low = boot.analysis_shift_CI95_boot_low;
    thisContrasts.analysis_shift_CI95_boot_high = boot.analysis_shift_CI95_boot_high;
    thisContrasts.bootstrap_success_n = boot.bootstrap_success_n;
    thisContrasts.bootstrap_success_fraction = boot.bootstrap_success_fraction;

    thisContrasts.CI95_low = thisContrasts.CI95_boot_low;
    thisContrasts.CI95_high = thisContrasts.CI95_boot_high;

    ContextualContrasts = [ContextualContrasts; thisContrasts];

    summaryRow = table( ...
        string(fam.name), ...
        string(fam.display), ...
        string(fam.rawVariable), ...
        string(fam.responseTransform), ...
        string(fam.formula0), ...
        string(fam.formula1), ...
        nObs, ...
        nUnits, ...
        nEvents, ...
        string(fam.incomeReference), ...
        fam.nZero, ...
        fam.nOne, ...
        fam.boundaryFraction, ...
        string(fam.boundaryHandling), ...
        R2m0, ...
        R2m1, ...
        deltaR2m, ...
        dist.rawQ25, ...
        dist.rawQ75, ...
        dist.rawIQR, ...
        dist.analysisQ25, ...
        dist.analysisQ75, ...
        dist.analysisIQR, ...
        dist.multiplicativeIQRRatio, ...
        thisDiag.residual_skewness, ...
        thisDiag.residual_excess_kurtosis, ...
        thisDiag.prop_abs_stdres_gt2, ...
        thisDiag.prop_abs_stdres_gt3, ...
        thisDiag.corr_absres_fitted, ...
        thisDiag.n_fitted_outside_bounds, ...
        height(thisContrasts), ...
        'VariableNames', { ...
            'indicator', ...
            'display', ...
            'raw_variable', ...
            'response_transform', ...
            'M0_formula', ...
            'M1_formula', ...
            'n', ...
            'n_units', ...
            'n_earthquakes', ...
            'income_reference', ...
            'n_exact_zero', ...
            'n_exact_one', ...
            'boundary_fraction', ...
            'boundary_handling', ...
            'R2m_M0', ...
            'R2m_M1', ...
            'delta_R2m', ...
            'raw_Q25', ...
            'raw_Q75', ...
            'raw_IQR', ...
            'analysis_Q25', ...
            'analysis_Q75', ...
            'analysis_IQR', ...
            'analysis_IQR_multiplicative_ratio', ...
            'residual_skewness', ...
            'residual_excess_kurtosis', ...
            'prop_abs_stdres_gt2', ...
            'prop_abs_stdres_gt3', ...
            'corr_absres_fitted', ...
            'n_fitted_outside_bounds', ...
            'n_contextual_contrasts'});

    ModelSummary = [ModelSummary; summaryRow];
end


resultsMat = fullfile(outputRoot, "contextual_stability_results.mat");

save( ...
    resultsMat, ...
    "ModelSummary", ...
    "ContextualContrasts");

fprintf("\n============================================================\n");
fprintf("Contextual-stability analysis complete.\n");
fprintf("Saved:\n  %s\n", resultsMat);
fprintf("============================================================\n");

end


function indicators = buildIndicators(Phase, Loss, Pair, cfg)

indicators = struct( ...
    "name", {}, ...
    "display", {}, ...
    "data", {}, ...
    "unitVar", {}, ...
    "rawVariable", {}, ...
    "formula0", {}, ...
    "formula1", {}, ...
    "responseTransform", {}, ...
    "incomeReference", {}, ...
    "nZero", {}, ...
    "nOne", {}, ...
    "boundaryFraction", {}, ...
    "boundaryHandling", {});

recordSpecs = { ...
    Phase, "phase_tau80_tau50", "tau80/tau50", "tau80_over_tau50"; ...
    Phase, "phase_tau90_tau80", "tau90/tau80", "tau90_over_tau80"; ...
    Phase, "phase_tau95_tau90", "tau95/tau90", "tau95_over_tau90"; ...
    Loss,  "compression_ci",    "Compression coefficient c_i", "c_i"};

for i = 1:size(recordSpecs,1)
    T = recordSpecs{i,1};
    T = addCommonPredictors(T, cfg);
    T.System = categoricalWithReference(string(T.System), cfg.systemOrder, "Power");
    T.event_ID = categorical(string(T.event_ID));
    T.record_ID = categorical(string(T.record_ID));

    rawVar = string(recordSpecs{i,4});
    raw = numericColumnFit(T.(rawVar));
    keep = isfinite(raw) & raw > 0;
    T = T(keep,:);
    raw = raw(keep);
    T.raw_response = raw;
    T.response = log(raw);

    T = completeContextRows(T, ["response","System","event_ID","record_ID"]);
    T = standardizeFixedCategories(T, cfg);

    formula0 = "response ~ System + (1|event_ID)";
    formula1 = "response ~ System + avg_PGA + event_decade + " + ...
        "IncomeGroup + log2_population + (1|event_ID)";

    indicators(end+1) = makeIndicator( ...
        string(recordSpecs{i,2}), string(recordSpecs{i,3}), T, ...
        "record_ID", rawVar, formula0, formula1, "log", cfg, ...
        nan, nan, nan, "not_bounded");
end

T = Loss;
T = addCommonPredictors(T, cfg);
T.System = categoricalWithReference(string(T.System), cfg.systemOrder, "Power");
T.event_ID = categorical(string(T.event_ID));
T.record_ID = categorical(string(T.record_ID));

raw = numericColumnFit(T.eta90);
keep = isfinite(raw) & raw > 0;
T = T(keep,:);
T.raw_response = raw(keep);
T.response = log(T.raw_response);
T = completeContextRows(T, ["response","System","event_ID","record_ID"]);
T = standardizeFixedCategories(T, cfg);
diag.nZero = sum(raw == 0, "omitnan");
diag.nOne = sum(raw == 1, "omitnan");
diag.boundaryFraction = mean(raw == 0 | raw == 1, "omitnan");
diag.transformName = "log";
diag.boundaryHandling = "positive_log_transform";

formula0 = "response ~ System + (1|event_ID)";
formula1 = "response ~ System + avg_PGA + event_decade + " + ...
    "IncomeGroup + log2_population + (1|event_ID)";

indicators(end+1) = makeIndicator( ...
    "tail_eta90", "Late-loss fraction eta90", T, "record_ID", "eta90", ...
    formula0, formula1, diag.transformName, cfg, ...
    diag.nZero, diag.nOne, diag.boundaryFraction, diag.boundaryHandling);

pairSpecs = { ...
    Pair, "pair_rho90", "Paired timing ratio rho90", "rho90_B_over_E"; ...
    Pair, "pair_rho95", "Paired timing ratio rho95", "rho95_B_over_E"};

for i = 1:size(pairSpecs,1)
    T = pairSpecs{i,1};
    T = addCommonPredictors(T, cfg);
    T.systemPair = categoricalWithReference(string(T.systemPair), cfg.pairOrder, "EW");
    T.event_ID = categorical(string(T.event_ID));
    T.pair_ID = categorical(string(T.pair_ID));

    rawVar = string(pairSpecs{i,4});
    raw = numericColumnFit(T.(rawVar));
    keep = isfinite(raw) & raw > 0;
    T = T(keep,:);
    raw = raw(keep);
    T.raw_response = raw;
    T.response = log(raw);

    T = completeContextRows(T, ["response","systemPair","event_ID","pair_ID"]);
    T = standardizeFixedCategories(T, cfg);

    formula0 = "response ~ systemPair + (1|event_ID)";
    formula1 = "response ~ systemPair + avg_PGA + event_decade + " + ...
        "IncomeGroup + log2_population + (1|event_ID)";

    indicators(end+1) = makeIndicator( ...
        string(pairSpecs{i,2}), string(pairSpecs{i,3}), T, ...
        "pair_ID", rawVar, formula0, formula1, "log", cfg, ...
        nan, nan, nan, "not_bounded");
end


end

function fam = makeIndicator(name, display, T, unitVar, rawVariable, ...
    formula0, formula1, responseTransform, cfg, nZero, nOne, ...
    boundaryFraction, boundaryHandling)

fam.name = char(name);
fam.display = char(display);
fam.data = T;
fam.unitVar = string(unitVar);
fam.rawVariable = char(rawVariable);
fam.formula0 = char(formula0);
fam.formula1 = char(formula1);
fam.responseTransform = char(responseTransform);
fam.nZero = nZero;
fam.nOne = nOne;
fam.boundaryFraction = boundaryFraction;
fam.boundaryHandling = char(boundaryHandling);

if ismember("IncomeGroup", string(T.Properties.VariableNames)) && ~isempty(T)
    cats = string(categories(T.IncomeGroup));
    if any(cats == cfg.incomeReference)
        fam.incomeReference = cfg.incomeReference;
    elseif ~isempty(cats)
        fam.incomeReference = cats(1);
    else
        fam.incomeReference = "";
    end
else
    fam.incomeReference = "";
end

end

function C = extractContextualCoefficients(M, fam, nObs, nUnits, nEvents, cfg)

[beta, betaNames, stats] = fixedEffects(M, "DFMethod", char(cfg.dfMethod));
names = string(betaNames.Name);
beta = double(beta(:));
ci = double(coefCI(M, "DFMethod", char(cfg.dfMethod)));
pValue = double(stats.pValue(:));

contextNames = ["avg_PGA","event_decade","IncomeGroup_High income","log2_population"];
isContext = ismember(names, contextNames);
idx = find(isContext);

C = table();
if isempty(idx)
    return;
end

C.indicator = repmat(string(fam.name), numel(idx), 1);
C.indicator_display = repmat(string(fam.display), numel(idx), 1);
C.response_transform = repmat(string(fam.responseTransform), numel(idx), 1);
C.n = repmat(nObs, numel(idx), 1);
C.n_units = repmat(nUnits, numel(idx), 1);
C.n_earthquakes = repmat(nEvents, numel(idx), 1);
C.term = names(idx);
C.term_display = arrayfun(@contextualCoefficientLabel, names(idx));
C.estimate = beta(idx);
C.CI95_model_low = ci(idx,1);
C.CI95_model_high = ci(idx,2);
C.model_df_method = repmat(string(cfg.dfMethod), numel(idx), 1);
C.p_value = pValue(idx);
C.significance = arrayfun(@significanceLabel, C.p_value);

end

function C = makeContextualContrasts(coefTable, fam, dist, T)

C = coefTable;
if isempty(C)
    return;
end

C.context = strings(height(C),1);
C.predictor_SD = nan(height(C),1);
C.predictor_contrast = nan(height(C),1);
C.contrast_definition = strings(height(C),1);

C.coefficient_estimate = C.estimate;
C.coefficient_CI95_model_low = C.CI95_model_low;
C.coefficient_CI95_model_high = C.CI95_model_high;

for i = 1:height(C)
    term = string(C.term(i));

    if term == "avg_PGA"
        x = double(T.avg_PGA);
        sx = finiteStd(x);
        C.context(i) = "PGA (2 SD)";
        C.predictor_SD(i) = sx;
        C.predictor_contrast(i) = 2*sx;
        C.contrast_definition(i) = "2 SD of avg_PGA";

    elseif term == "event_decade"
        x = double(T.event_decade);
        sx = finiteStd(x);
        C.context(i) = "Event year (2 SD)";
        C.predictor_SD(i) = sx;
        C.predictor_contrast(i) = 2*sx;
        C.contrast_definition(i) = "2 SD of event year";

    elseif term == "log2_population"
        x = double(T.log2_population);
        sx = finiteStd(x);
        C.context(i) = "Population (2 SD)";
        C.predictor_SD(i) = sx;
        C.predictor_contrast(i) = 2*sx;
        C.contrast_definition(i) = "2 SD of log2(population)";

    elseif term == "IncomeGroup_High income"
        C.context(i) = "High vs non-high income";
        C.predictor_SD(i) = NaN;
        C.predictor_contrast(i) = 1;
        C.contrast_definition(i) = "binary 0-to-1 contrast";

    else
        C.context(i) = term;
        C.predictor_contrast(i) = 1;
        C.contrast_definition(i) = "unit contrast";
    end
end

C.analysis_shift_estimate = C.coefficient_estimate .* C.predictor_contrast;
C.analysis_shift_CI95_low = C.coefficient_CI95_model_low .* C.predictor_contrast;
C.analysis_shift_CI95_high = C.coefficient_CI95_model_high .* C.predictor_contrast;

C.effect_scale = repmat("original_scale_multiplicative_ratio", height(C),1);
C.backtransform_type = repmat("multiplicative_ratio", height(C), 1);
C.backtransformed_estimate = exp(C.analysis_shift_estimate);
C.backtransformed_CI95_model_low = exp(C.analysis_shift_CI95_low);
C.backtransformed_CI95_model_high = exp(C.analysis_shift_CI95_high);

C.estimate = C.backtransformed_estimate;
C.CI95_model_low = C.backtransformed_CI95_model_low;
C.CI95_model_high = C.backtransformed_CI95_model_high;

C.original_scale_percent_change = 100 .* (C.backtransformed_estimate - 1);

C.analysis_IQR = repmat(dist.analysisIQR, height(C), 1);
C.raw_Q25 = repmat(dist.rawQ25, height(C), 1);
C.raw_Q75 = repmat(dist.rawQ75, height(C), 1);
C.raw_IQR = repmat(dist.rawIQR, height(C), 1);
C.analysis_Q25 = repmat(dist.analysisQ25, height(C), 1);
C.analysis_Q75 = repmat(dist.analysisQ75, height(C), 1);

end

function B = bootstrapContextualContrasts(T, fam, contrasts, cfg)

nContrast = height(contrasts);
B = table();
B.CI95_boot_low = nan(nContrast,1);
B.CI95_boot_high = nan(nContrast,1);
B.analysis_shift_CI95_boot_low = nan(nContrast,1);
B.analysis_shift_CI95_boot_high = nan(nContrast,1);
B.bootstrap_success_n = zeros(nContrast,1);
B.bootstrap_success_fraction = zeros(nContrast,1);

if isempty(T) || nContrast == 0
    return;
end

TbootBase = T;
TbootBase.event_ID = string(TbootBase.event_ID);
eventLabels = unique(TbootBase.event_ID, "stable");
nEvents = numel(eventLabels);
if nEvents < 2
    warning("%s has fewer than two earthquakes; cluster bootstrap is unavailable.", fam.name);
    return;
end

nBoot = cfg.bootstrapB;
bootShift = nan(nBoot,nContrast);
bootRatio = nan(nBoot,nContrast);

oldRng = rng;
cleanupObj = onCleanup(@() rng(oldRng));
nameOffset = sum(double(char(string(fam.name))));
rng(double(cfg.bootstrapSeed) + nameOffset, "twister");

for b = 1:nBoot
    draw = randi(nEvents,nEvents,1);
    Tb = TbootBase([],:);

    for k = 1:nEvents
        sourceEvent = eventLabels(draw(k));
        block = TbootBase(TbootBase.event_ID == sourceEvent,:);
        block.event_ID = repmat("boot_event_" + string(k), height(block), 1);
        Tb = [Tb; block];
    end

    Tb.event_ID = categorical(Tb.event_ID);
    Tb = standardizeFixedCategories(Tb, cfg);

    try
        Mb = fitLmeChecked(Tb, fam.formula1);
        [bb, bn] = fixedEffects(Mb);
        bNames = string(bn.Name);
        bb = double(bb(:));

        for j = 1:nContrast
            jj = find(bNames == string(contrasts.term(j)), 1, "first");
            if isempty(jj) || ~isfinite(bb(jj))
                continue;
            end

            shift = bb(jj) .* contrasts.predictor_contrast(j);
            if ~isfinite(shift)
                continue;
            end

            bootShift(b,j) = shift;
            if string(fam.responseTransform) == "log"
                bootRatio(b,j) = exp(shift);
            else
                bootRatio(b,j) = shift;
            end
        end
    catch
    end
end

for j = 1:nContrast
    okShift = isfinite(bootShift(:,j));
    okRatio = isfinite(bootRatio(:,j));
    nSuccess = sum(okShift & okRatio);
    frac = nSuccess / nBoot;
    B.bootstrap_success_n(j) = nSuccess;
    B.bootstrap_success_fraction(j) = frac;

    if frac < cfg.bootstrapMinSuccessFraction
        warning('%s | %s: earthquake-cluster bootstrap success %.1f%% (%d/%d), below the prespecified %.1f%% threshold; CI set to NaN.', ...
            fam.name, char(string(contrasts.context(j))), 100*frac, nSuccess, nBoot, ...
            100*cfg.bootstrapMinSuccessFraction);
        continue;
    end

    B.analysis_shift_CI95_boot_low(j) = prctile(bootShift(okShift,j), cfg.bootstrapPercentiles(1));
    B.analysis_shift_CI95_boot_high(j) = prctile(bootShift(okShift,j), cfg.bootstrapPercentiles(2));
    B.CI95_boot_low(j) = prctile(bootRatio(okRatio,j), cfg.bootstrapPercentiles(1));
    B.CI95_boot_high(j) = prctile(bootRatio(okRatio,j), cfg.bootstrapPercentiles(2));
end

end

function s = finiteStd(x)
x = double(x(:));
x = x(isfinite(x));
if numel(x) < 2
    s = NaN;
else
    s = std(x,0);
end
end

function dist = computeEmpiricalDistribution(raw, analysis, fam)

raw = double(raw(:));
analysis = double(analysis(:));
raw = raw(isfinite(raw));
analysis = analysis(isfinite(analysis));

rq = prctile(raw,[25 75]);
aq = prctile(analysis,[25 75]);

dist.rawQ25 = rq(1);
dist.rawQ75 = rq(2);
dist.rawIQR = rq(2)-rq(1);
dist.analysisQ25 = aq(1);
dist.analysisQ75 = aq(2);
dist.analysisIQR = aq(2)-aq(1);

if string(fam.responseTransform) == "log"
    dist.multiplicativeIQRRatio = exp(dist.analysisIQR);
else
    dist.multiplicativeIQRRatio = NaN;
end

end

function label = contextualCoefficientLabel(name)

label = string(name);
label = replace(label, "IncomeGroup_High income", "High income vs non-high income");
label = replace(label, "avg_PGA", "Average PGA");
label = replace(label, "log2_population", "log2(Population)");
label = replace(label, "event_decade", "Event year (+10 years)");

end

function label = significanceLabel(p)

if ~isfinite(p)
    label = "";
elseif p < 0.001
    label = "***";
elseif p < 0.01
    label = "**";
elseif p < 0.05
    label = "*";
else
    label = "";
end

end

function T = addCommonPredictors(T, cfg)

T.avg_PGA = numericColumnFit(T.avg_PGA);
T.population = numericColumnFit(T.population);
T.event_year = numericColumnFit(T.event_year);
T.avg_PGA(~isfinite(T.avg_PGA) | T.avg_PGA < 0) = NaN;

T.log2_population = nan(height(T),1);
ok = isfinite(T.population) & T.population > 0;
T.log2_population(ok) = log2(T.population(ok));

T.event_decade = (T.event_year - cfg.centerYear) ./ 10;
T.IncomeGroup = categoricalWithReference( ...
    binaryIncomeGroup(string(T.IncomeGroup)), ...
    cfg.incomeOrder, cfg.incomeReference);

end

function T = completeContextRows(T, structuralVars)

requiredNumeric = ["response","avg_PGA","event_decade","log2_population"];
keep = true(height(T),1);

for i = 1:numel(requiredNumeric)
    v = numericColumnFit(T.(requiredNumeric(i)));
    keep = keep & isfinite(v);
end

keep = keep & ~isundefined(T.IncomeGroup);

for i = 1:numel(structuralVars)
    name = structuralVars(i);
    if ~ismember(name, string(T.Properties.VariableNames))
        continue;
    end
    v = T.(name);
    if isnumeric(v)
        keep = keep & isfinite(v);
    elseif iscategorical(v)
        keep = keep & ~isundefined(v);
    else
        keep = keep & ~ismissing(v);
    end
end

T = T(keep,:);

end

function T = standardizeFixedCategories(T, cfg)

if ismember("System", string(T.Properties.VariableNames))
    T.System = categoricalWithReference(string(T.System), cfg.systemOrder, "Power");
end
if ismember("systemPair", string(T.Properties.VariableNames))
    T.systemPair = categoricalWithReference(string(T.systemPair), cfg.pairOrder, "EW");
end
T.IncomeGroup = categoricalWithReference( ...
    binaryIncomeGroup(string(T.IncomeGroup)), cfg.incomeOrder, cfg.incomeReference);
T.event_ID = categorical(string(T.event_ID));
if ismember("record_ID", string(T.Properties.VariableNames))
    T.record_ID = categorical(string(T.record_ID));
end
if ismember("pair_ID", string(T.Properties.VariableNames))
    T.pair_ID = categorical(string(T.pair_ID));
end

end

function M = fitLmeChecked(T, formula)

lastwarn('');
M = fitlme(T, char(formula), ...
    "FitMethod", "REML", ...
    "DummyVarCoding", "reference");

[warnMsg,~] = lastwarn;
warnText = lower(string(warnMsg));
if contains(warnText,"rank deficient") || contains(warnText,"rank-deficient")
    error("Rank-deficient mixed model:\n%s", formula);
end

end

function D = computeResidualDiagnostics(M, T, fam)

fitted = double(predict(M, T, "Conditional", true));
y = double(T.response);
resid = y - fitted;

ok = isfinite(y) & isfinite(fitted) & isfinite(resid);
y = y(ok);
fitted = fitted(ok);
resid = resid(ok);

n = numel(resid);
if n == 0
    mu = NaN; sd = NaN; skew = NaN; exkurt = NaN;
    p2 = NaN; p3 = NaN; corAbs = NaN;
else
    mu = mean(resid);
    sd = std(resid,0);
    if isfinite(sd) && sd > 0
        z = (resid - mu) ./ sd;
    else
        z = nan(size(resid));
    end

    centered = resid - mu;
    m2 = mean(centered.^2);
    if isfinite(m2) && m2 > 0
        m3 = mean(centered.^3);
        m4 = mean(centered.^4);
        skew = m3 / (m2^(3/2));
        exkurt = m4 / (m2^2) - 3;
    else
        skew = NaN;
        exkurt = NaN;
    end

    p2 = mean(abs(z) > 2, "omitnan");
    p3 = mean(abs(z) > 3, "omitnan");

    if numel(fitted) >= 3 && std(abs(resid),0) > 0 && std(fitted,0) > 0
        C = corrcoef(abs(resid), fitted, "Rows", "complete");
        if all(size(C) == [2 2])
            corAbs = C(1,2);
        else
            corAbs = NaN;
        end
    else
        corAbs = NaN;
    end
end

if string(fam.responseTransform) == "identity"
    outside = fitted < 0 | fitted > 1;
    nOutside = sum(outside);
    propOutside = mean(outside);
elseif string(fam.name) == "tail_eta90" && string(fam.responseTransform) == "log"
    fittedRaw = exp(fitted);
    outside = fittedRaw <= 0 | fittedRaw > 1;
    nOutside = sum(outside);
    propOutside = mean(outside);
else
    nOutside = 0;
    propOutside = 0;
end

D = table( ...
    string(fam.name), string(fam.display), string(fam.responseTransform), ...
    n, mean(fitted,"omitnan"), min(fitted,[],"omitnan"), max(fitted,[],"omitnan"), ...
    mu, sd, skew, exkurt, p2, p3, corAbs, nOutside, propOutside, ...
    'VariableNames', { ...
        'indicator','display','response_transform','n', ...
        'fitted_mean','fitted_min','fitted_max', ...
        'residual_mean','residual_sd','residual_skewness', ...
        'residual_excess_kurtosis','prop_abs_stdres_gt2', ...
        'prop_abs_stdres_gt3','corr_absres_fitted', ...
        'n_fitted_outside_bounds','prop_fitted_outside_bounds'});

end

function R2m = marginalR2(M, T)

fixedPrediction = predict(M,T,"Conditional",false);
fixedVariance = var(double(fixedPrediction),1,"omitnan");
[psi,mse] = covarianceParameters(M);
randomVariance = 0;

for r = 1:numel(psi)
    P = double(psi{r});
    if isempty(P)
        continue;
    end
    if ~isscalar(P)
        error("marginalR2 expects random-intercept-only models.");
    end
    randomVariance = randomVariance + P;
end

mse = double(mse);
denominator = fixedVariance + randomVariance + mse;
if ~isfinite(denominator) || denominator <= 0
    R2m = NaN;
else
    R2m = fixedVariance/denominator;
end

end

function x = numericColumnFit(v)
if isnumeric(v) || islogical(v)
    x = double(v);
else
    x = str2double(string(v));
end
end

function labels = binaryIncomeGroup(labels)
labels = strtrim(string(labels));
labels(ismissing(labels)) = "";
isHigh = strcmpi(labels,"High income") | strcmpi(labels,"High");
isKnownNonHigh = strcmpi(labels,"Low income") | strcmpi(labels,"Low") | ...
    strcmpi(labels,"Lower middle income") | strcmpi(labels,"Lower_middle") | ...
    strcmpi(labels,"Upper middle income") | strcmpi(labels,"Upper_middle") | ...
    strcmpi(labels,"Non-high income");
out = strings(size(labels));
out(isHigh) = "High income";
out(isKnownNonHigh) = "Non-high income";
labels = out;
end

function c = categoricalWithReference(labels, preferredOrder, reference)

labels = strtrim(string(labels));
labels(ismissing(labels)) = "";
present = unique(labels(labels ~= ""), "stable");
preferredOrder = string(preferredOrder);
reference = string(reference);

if any(present == reference)
    ordered = [reference, preferredOrder(preferredOrder ~= reference & ismember(preferredOrder,present))];
else
    ordered = preferredOrder(ismember(preferredOrder,present));
end
ordered = [ordered, present(~ismember(present,ordered))'];

if isempty(ordered)
    c = categorical(labels);
else
    c = categorical(labels,cellstr(ordered),"Ordinal",false);
end

end


function ensureDirFit(d)
if exist(d,"dir") ~= 7
    mkdir(d);
end
end


function draw_FigS12_contextual_stability(supplementaryRoot)

if nargin < 1 || isempty(supplementaryRoot)
    moduleRoot = string(fileparts(mfilename("fullpath")));
    emergentSimplicityRoot = string(fileparts(fileparts(moduleRoot)));
    supplementaryRoot = fullfile(emergentSimplicityRoot, "outputs", "supplementary");
end
scriptDir = string(fileparts(mfilename("fullpath")));
moduleRoot = scriptDir;
sourceDataDir = fullfile(supplementaryRoot, "source_data");
figureDir = fullfile(supplementaryRoot, "figures");
ensureDirPlot(sourceDataDir);
ensureDirPlot(figureDir);

resultsMat = fullfile(sourceDataDir,"contextual_stability_results.mat");
if ~isfile(resultsMat)
    error("Missing results file: %s\nRun run_01 first.",resultsMat);
end

S = load(resultsMat, ...
    "ModelSummary","ContextualContrasts");

ModelSummary = S.ModelSummary;
ContextualContrasts = S.ContextualContrasts;

removeObsoleteOutputs(figureDir);

sourceXlsx = fullfile(sourceDataDir, "SourceData_FigS12.xlsx");
if isfile(sourceXlsx)
    delete(sourceXlsx);
end

writetable(modelSummarySource(ModelSummary), sourceXlsx, "Sheet", "model_summary");
writetable(contextualResultsSource(ContextualContrasts), sourceXlsx, "Sheet", "contextual_results");
writetable(contextualContrastSource(ContextualContrasts), sourceXlsx, "Sheet", "contrast_plot_source");

oldR2Source = fullfile(sourceDataDir,"Contextual_R2_plot_source.csv");
if isfile(oldR2Source)
    delete(oldR2Source);
end

oldContrastCsv = fullfile(sourceDataDir,"Contextual_contrasts_plot_source.csv");
if isfile(oldContrastCsv)
    delete(oldContrastCsv);
end

oldResultsXlsx = fullfile(sourceDataDir,"Contextual_stability_results.xlsx");
if isfile(oldResultsXlsx)
    delete(oldResultsXlsx);
end

plotContextualPanels( ...
    ContextualContrasts, ...
    ModelSummary, ...
    figureDir);

preparedMat = fullfile(sourceDataDir, "contextual_stability_data.mat");
if isfile(preparedMat)
    delete(preparedMat);
end
if isfile(resultsMat)
    delete(resultsMat);
end

fprintf("\nSaved contextual-stability figures to:\n  %s\n",figureDir);


end

function plotContextualPanels(E,M,outDir)

indicatorOrder = [ ...
    "phase_tau80_tau50", ...
    "phase_tau90_tau80", ...
    "phase_tau95_tau90", ...
    "compression_ci", ...
    "tail_eta90", ...
    "pair_rho90", ...
    "pair_rho95"];

panelLetters = ["a","b","c","d","e","f","g","h"];

fig = figure( ...
    "Color","w", ...
    "Units","pixels", ...
    "Position",[60 35 1900 1420], ...
    "Renderer","painters");

tl = tiledlayout( ...
    fig, ...
    3, ...
    3, ...
    "TileSpacing","loose", ...
    "Padding","compact");

panelTiles = [1 2 3 4 5 6 7];

for i = 1:numel(indicatorOrder)

    ax = nexttile(tl,panelTiles(i));

    D = E( ...
        string(E.indicator) == indicatorOrder(i), ...
        :);

    drawContrastPanel(ax,D);
    addPanelLetter(ax,panelLetters(i));

end

ax = nexttile(tl,8,[1 2]);

drawR2Panel( ...
    ax, ...
    M, ...
    indicatorOrder);

addPanelLetter(ax,panelLetters(8));

outBase = fullfile( ...
    outDir, ...
    "FigS12_contextual_stability");

exportFigure(fig,outBase);

end

function drawContrastPanel(ax,D)

hold(ax,"on");
box(ax,"off");

if isempty(D)

    text( ...
        ax, ...
        0.5, ...
        0.5, ...
        "No estimates", ...
        "Units","normalized", ...
        "HorizontalAlignment","center", ...
        "FontName","Arial", ...
        "FontSize",16);

    return;

end

D = orderContextRows(D);

y = 1 + (0:height(D)-1)' * 1.00;

allX = [ ...
    double(D.CI95_low); ...
    double(D.CI95_high); ...
    1];

allX = allX(isfinite(allX));
if isempty(allX)
    allX = [0.5;1.5];
end

xMin = min(allX);
xMax = max(allX);
span = max(xMax-xMin,0.2);
xMin = max(0, xMin - 0.10*span);
xMax = xMax + 0.28*span;

xline( ...
    ax, ...
    1, ...
    "-", ...
    "Color",[0.40 0.40 0.40], ...
    "LineWidth",0.9, ...
    "HandleVisibility","off");

for i = 1:height(D)

    plot( ...
        ax, ...
        [D.CI95_low(i) D.CI95_high(i)], ...
        [y(i) y(i)], ...
        "-", ...
        "Color",[0.20 0.43 0.68], ...
        "LineWidth",1.9);

    plot( ...
        ax, ...
        D.estimate(i), ...
        y(i), ...
        "o", ...
        "MarkerSize",8.0, ...
        "MarkerFaceColor",[0.20 0.43 0.68], ...
        "MarkerEdgeColor","w", ...
        "LineWidth",0.7);

    coefLabel = sprintf( ...
        "%.3f%s", ...
        double(D.estimate(i)), ...
        char(string(D.significance(i))));

    labelX = ...
        D.CI95_high(i) + ...
        0.022*span;

    text( ...
        ax, ...
        labelX, ...
        y(i), ...
        coefLabel, ...
        "FontName","Arial", ...
        "FontSize",15.0, ...
        "VerticalAlignment","middle", ...
        "HorizontalAlignment","left");

end

labels = string(D.context);
labels(labels == "High vs non-high income") = "Income group";
labels(labels == "Event year (2 SD)") = "Earthquake year (2 SD)";
labels(labels == "Population (2 SD)") = "Population (2 SD)";
labels(labels == "PGA (2 SD)") = "PGA (2 SD)";

set( ...
    ax, ...
    "YTick",y, ...
    "YTickLabel",cellstr(labels), ...
    "YDir","reverse", ...
    "TickLength",[0 0], ...
    "FontName","Arial", ...
    "FontSize",15.0, ...
    "Layer","top", ...
    "TickLabelInterpreter","none");

ylim(ax,[0.60 y(end)+0.40]);
xlim(ax,[xMin xMax]);

grid(ax,"off");

indicatorLabel = contrastIndicatorLabel(string(D.indicator(1)));

xlabel( ...
    ax, ...
    "Multiplicative change in " + indicatorLabel, ...
    "FontName","Arial", ...
    "FontSize",15.5, ...
    "Interpreter","latex");


end

function label = contrastIndicatorLabel(indicator)

switch string(indicator)
    case "phase_tau80_tau50"
        label = "$\tau_{80}/\tau_{50}$";
    case "phase_tau90_tau80"
        label = "$\tau_{90}/\tau_{80}$";
    case "phase_tau95_tau90"
        label = "$\tau_{95}/\tau_{90}$";
    case "compression_ci"
        label = "$\kappa$";
    case "tail_eta90"
        label = "$\eta_{90}$";
    case "pair_rho90"
        label = "$\rho_{90}$";
    case "pair_rho95"
        label = "$\rho_{95}$";
    otherwise
        label = string(indicator);
end

end

function D = orderContextRows(D)

wanted = [ ...
    "Event year (2 SD)", ...
    "High vs non-high income", ...
    "Population (2 SD)", ...
    "PGA (2 SD)"];

rank = 999*ones(height(D),1);

for i = 1:height(D)

    hit = find( ...
        wanted == string(D.context(i)), ...
        1);

    if ~isempty(hit)
        rank(i) = hit;
    end

end

[~,idx] = sort(rank);

D = D(idx,:);

end


function drawR2Panel(ax,M,indicatorOrder)

rows = nan(numel(indicatorOrder),1);

for i = 1:numel(indicatorOrder)

    hit = find( ...
        string(M.indicator) == indicatorOrder(i), ...
        1);

    if ~isempty(hit)
        rows(i) = hit;
    end

end

rows = rows(isfinite(rows));

M = M(rows,:);

hold(ax,"on");
box(ax,"off");

if isempty(M)

    text( ...
        ax, ...
        0.5, ...
        0.5, ...
        "No R2 results", ...
        "Units","normalized", ...
        "HorizontalAlignment","center", ...
        "FontName","Arial", ...
        "FontSize",16);

    return;

end

y = 1:height(M);

offset = 0.17;

barh( ...
    ax, ...
    y-offset, ...
    double(M.R2m_M0), ...
    0.30, ...
    "FaceColor",[0.72 0.72 0.72], ...
    "EdgeColor","none", ...
    "DisplayName","Structural model");

barh( ...
    ax, ...
    y+offset, ...
    double(M.R2m_M1), ...
    0.30, ...
    "FaceColor",[0.20 0.43 0.68], ...
    "EdgeColor","none", ...
    "DisplayName","Contextual model");

for i = 1:height(M)

    x = max( ...
        M.R2m_M0(i), ...
        M.R2m_M1(i));

    text( ...
        ax, ...
        x+0.004, ...
        y(i), ...
        sprintf( ...
            "\\Delta R_m^2 = %.3f", ...
            M.delta_R2m(i)), ...
        "FontName","Arial", ...
        "FontSize",12.5, ...
        "VerticalAlignment","middle", ...
        "Interpreter","tex");

end

labels = shortIndicatorLabels( ...
    string(M.indicator));

set( ...
    ax, ...
    "YTick",y, ...
    "YTickLabel",cellstr(labels), ...
    "YDir","reverse", ...
    "FontName","Arial", ...
    "FontSize",14.5, ...
    "TickLength",[0 0], ...
    "TickLabelInterpreter","latex");

xMax = max( ...
    [M.R2m_M0;M.R2m_M1], ...
    [], ...
    "omitnan");

xlim( ...
    ax, ...
    [0 max(0.20,xMax+0.12)]);

xlabel( ...
    ax, ...
    "R^2", ...
    "FontName","Arial", ...
    "FontSize",15.5, ...
    "Interpreter","tex");

legend( ...
    ax, ...
    "Location","southeast", ...
    "Box","off", ...
    "FontName","Arial", ...
    "FontSize",12.5);

grid(ax,"off");

end

function labels = shortIndicatorLabels(indicators)

labels = indicators;

labels(labels == "phase_tau80_tau50") = "$\tau_{80}/\tau_{50}$";
labels(labels == "phase_tau90_tau80") = "$\tau_{90}/\tau_{80}$";
labels(labels == "phase_tau95_tau90") = "$\tau_{95}/\tau_{90}$";

labels(labels == "compression_ci") = "$\kappa$";
labels(labels == "tail_eta90") = "$\eta_{90}$";
labels(labels == "pair_rho90") = "$\rho_{90}$";
labels(labels == "pair_rho95") = "$\rho_{95}$";

end

function addPanelLetter(ax,letter)

text( ...
    ax, ...
    -0.08, ...
    1.05, ...
    string(letter), ...
    "Units","normalized", ...
    "FontName","Arial", ...
    "FontSize",18, ...
    "FontWeight","bold", ...
    "HorizontalAlignment","left", ...
    "VerticalAlignment","bottom", ...
    "Clipping","off");

end


function removeObsoleteOutputs(outDir)

names = [ ...
    "ExtendedData_contextual_coefficients", ...
    "ExtendedData_contextual_stability_bound", ...
    "ExtendedData_contextual_R2_change", ...
    "ExtendedData_contextual_R2_comparison"];

exts = [ ...
    ".fig", ...
    ".png", ...
    ".pdf"];

for i = 1:numel(names)

    for j = 1:numel(exts)

        p = fullfile( ...
            outDir, ...
            names(i)+exts(j));

        if isfile(p)
            delete(p);
        end

    end

end

oldTables = [ ...
    "Contextual_coefficients_plot_source.csv", ...
    "Contextual_effects_plot_source.csv", ...
    "Contextual_R2_change_plot_source.csv"];

for i = 1:numel(oldTables)

    p = fullfile( ...
        outDir, ...
        oldTables(i));

    if isfile(p)
        delete(p);
    end

end

end

function exportFigure(fig,outBase)

cleanupObj = onCleanup(@() closeFigureIfValid(fig));

savefig( ...
    fig, ...
    outBase+".fig");

exportgraphics( ...
    fig, ...
    outBase+".png", ...
    "Resolution",600, ...
    "BackgroundColor","white");

try

    exportgraphics( ...
        fig, ...
        outBase+".pdf", ...
        "ContentType","vector", ...
        "BackgroundColor","white");

catch

    print( ...
        fig, ...
        char(outBase+".pdf"), ...
        "-dpdf", ...
        "-painters");

end

end

function closeFigureIfValid(fig)
if isgraphics(fig)
    close(fig);
end
end

function T = contextualContrastSource(C)

keepVars = [ ...
    "indicator","indicator_display","response_transform","context", ...
    "n","n_units","n_earthquakes","predictor_SD","predictor_contrast", ...
    "effect_scale","estimate","CI95_low","CI95_high","model_df_method", ...
    "CI95_model_low","CI95_model_high","CI95_boot_low","CI95_boot_high", ...
    "p_value","analysis_shift_estimate","analysis_shift_CI95_low","analysis_shift_CI95_high", ...
    "analysis_shift_CI95_boot_low","analysis_shift_CI95_boot_high", ...
    "bootstrap_success_n","bootstrap_success_fraction", ...
    "backtransform_type","original_scale_percent_change"];

keepVars = keepVars(ismember(keepVars, string(C.Properties.VariableNames)));
T = C(:, keepVars);

end

function T = modelSummarySource(M)

T = M(:, [ ...
    "indicator","display","response_transform","n","n_earthquakes", ...
    "R2m_M0","R2m_M1","delta_R2m", ...
    "residual_skewness","residual_excess_kurtosis", ...
    "prop_abs_stdres_gt2","prop_abs_stdres_gt3", ...
    "corr_absres_fitted"]);

end

function T = contextualResultsSource(C)

T = C(:, [ ...
    "indicator","indicator_display","response_transform","context", ...
    "n","n_earthquakes","predictor_SD","predictor_contrast", ...
    "effect_scale","estimate","CI95_low","CI95_high","model_df_method", ...
    "CI95_model_low","CI95_model_high","CI95_boot_low","CI95_boot_high", ...
    "p_value","analysis_shift_estimate","analysis_shift_CI95_low","analysis_shift_CI95_high", ...
    "analysis_shift_CI95_boot_low","analysis_shift_CI95_boot_high", ...
    "bootstrap_success_n","bootstrap_success_fraction", ...
    "backtransform_type","original_scale_percent_change"]);

end

function ensureDirPlot(d)

if exist(d,"dir") ~= 7
    mkdir(d);
end

end

