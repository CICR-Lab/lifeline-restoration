function run_01_fit_restoration_shape_models()




scriptDir = string(fileparts(mfilename("fullpath")));
moduleRoot = scriptDir;
emergentSimplicityRoot = string(fileparts(fileparts(moduleRoot)));
outputDir = fullfile( ...
    emergentSimplicityRoot, ...
    "outputs", ...
    "intermediate", ...
    "01_restoration_shape_phase_ratio");

ensureDir(outputDir);


opts = struct();

opts.maeThreshold = 0.05;
opts.minFitPoints = 5;
opts.mainMinPoints = 6;
opts.nRestarts = 25;


opts.initPositiveMultipliers = [0.1, 0.3, 1, 3, 10];
opts.initUnitIntervalFractions = [0.1, 0.3, 0.5, 0.7, 0.9];


opts.positiveLowerBound = 1e-12;

opts.lambdaBounds = [opts.positiveLowerBound, Inf];
opts.kBounds = [opts.positiveLowerBound, Inf];
opts.bBounds = [opts.positiveLowerBound, Inf];
opts.cBounds = [opts.positiveLowerBound, Inf];
opts.wBounds = [opts.positiveLowerBound, Inf];
opts.gammaBounds = [-1 + opts.positiveLowerBound, Inf];
opts.aBounds = [opts.positiveLowerBound, 1];

opts.modelOrder = [ ...
    "Weibull", ...
    "Gompertz", ...
    "Log-logistic", ...
    "Exponential", ...
    "Mechanical analogy", ...
    "Opabola", ...
    "Zorn", ...
    "Trigonometric"];

opts.selectedSystems = normalizeSelectedSystems( ...
    ["Power", "Water", "Gas"]);

curves = loadPreparedRecords(moduleRoot);


curves = curves( ...
    [curves.nRestorationRatioPoints] >= opts.minFitPoints & ...
    [curves.nUniqueRestorationValues] >= 3);

curves = filterCurvesBySystem( ...
    curves, ...
    opts.selectedSystems);

if isempty(curves)
    error( ...
        "No records satisfy the fitting requirements after system filtering.");
end


fprintf( ...
    "\n=== Fitting restoration-shape models: >=%d points and >=3 unique values ===\n", ...
    opts.minFitPoints);

result = fitCurveDataset( ...
    curves, ...
    "restoration_shape_ge5", ...
    opts);

mainSummary = buildRestorationShapeSummaryGe6( ...
    result, ...
    opts);

result.RestorationShapeSummaryGe6 = mainSummary;

badfitCandidates = buildBadFitWeibullGompertzLoglogisticCandidateTable( ...
    result, ...
    curves, ...
    opts);


writeFitOutputs( ...
    result, ...
    outputDir);

restorationSummaryPath = fullfile( ...
    outputDir, ...
    "restoration_shape_ge5_fit_result.mat");

badfitPath = fullfile( ...
    outputDir, ...
    "badfit_weibull_gompertz_loglogistic_ge6_candidates.xlsx");

sourceDataPath = fullfile( ...
    outputDir, ...
    "SourceData_Fig3a_restoration_shape.xlsx");

writeRestorationShapeSourceData( ...
    mainSummary, ...
    sourceDataPath);

if isfile(badfitPath)
    delete(badfitPath);
end

writetable( ...
    badfitCandidates, ...
    badfitPath, ...
    "Sheet", "badfit_candidates");

if ~isfile(restorationSummaryPath)
    error( ...
        "Restoration-shape MAT file was not created:\n%s", ...
        restorationSummaryPath);
end

nGe5 = nnz( ...
    [curves.nRestorationRatioPoints] >= opts.minFitPoints);

nGe6 = nnz( ...
    [curves.nRestorationRatioPoints] >= opts.mainMinPoints);

fprintf( ...
    "Fitted sensitivity cohort (>=5 points): %d\n", ...
    nGe5);

fprintf( ...
    "Main cohort (>=6 points): %d\n", ...
    nGe6);

fprintf( ...
    "Weibull/Gompertz/Log-logistic bad-fit candidates " + ...
    "(all three finite MAE > %.3f): %d\n", ...
    opts.maeThreshold, ...
    height(badfitCandidates));

fprintf( ...
    "Saved outputs under:\n%s\n", ...
    outputDir);

fprintf( ...
    "Saved restoration-shape MAT output:\n%s\n", ...
    restorationSummaryPath);

fprintf( ...
    "Saved restoration-shape source data:\n%s\n", ...
    sourceDataPath);

fprintf("\nDone.\n");

end



function result = fitCurveDataset(curves, datasetLabel, opts)

summaryRows = {};
fitRows = {};
bestRows = {};

for k = 1:numel(curves)

    rec = curves(k);

    [tAbs, rObs] = preparedCurveToArrays(rec.restorationratioCurve);

    if numel(tAbs) < 2

        summaryRows{end + 1, 1} = makeSkippedSummary( ...
            rec, ...
            datasetLabel, ...
            "too_few_restoration_points");

        continue;
    end

    peakTime = getNumericField(rec, "PeakTime");
    tau100 = getNumericField(rec, "tau100");

    if ~isfinite(peakTime)

        summaryRows{end + 1, 1} = makeSkippedSummary( ...
            rec, ...
            datasetLabel, ...
            "missing_PeakTime");

        continue;
    end

    tau = tAbs - peakTime;

    hasTau100 = isfinite(tau100) && tau100 > 0;

    fits = struct([]);

    modelOrder = buildModelOrderForCurve( ...
        opts.modelOrder, ...
        hasTau100);

    for m = 1:numel(modelOrder)

        modelName = modelOrder(m);

        fit = fitOneModelR( ...
            tau, ...
            rObs, ...
            modelName, ...
            opts, ...
            tau100);

        met = computeMetrics( ...
            rObs, ...
            fit.rhat, ...
            fit.SSE, ...
            fit.m, ...
            fit.maskFit);

        fitRow = makeFitRow( ...
            rec, ...
            datasetLabel, ...
            modelName, ...
            fit, ...
            met, ...
            peakTime);

        fitRows{end + 1, 1} = fitRow;

        fits = appendStruct(fits, fitRow);
    end

    bestRow = chooseBestModelRow(fits);

    bestRows{end + 1, 1} = bestRow;

    summaryRow = makeSummaryRow( ...
        rec, ...
        datasetLabel, ...
        fits, ...
        bestRow, ...
        peakTime, ...
        opts.maeThreshold);

    summaryRows{end + 1, 1} = summaryRow;


end

result = struct();

result.CurveSummary = vertcatOrEmpty( ...
    summaryRows, ...
    emptySummaryTable());

result.CurveFitParams = vertcatOrEmpty( ...
    fitRows, ...
    emptyFitTable());

result.BestModelPerCurve = vertcatOrEmpty( ...
    bestRows, ...
    emptyBestTable());

result.datasetLabel = datasetLabel;
result.opts = opts;

end


function writeFitOutputs(result, workingDir)

ensureDir(workingDir);

out = struct("fitResult", result);
save( ...
    fullfile(workingDir, "restoration_shape_ge5_fit_result.mat"), ...
    "-struct", ...
    "out");

end

function T = buildRestorationShapeSummaryGe6(result, opts)

models = opts.modelOrder;

groups = [ ...
    "All systems", ...
    "Power", ...
    "Water", ...
    "Gas"];

groupCodes = [ ...
    "ALL", ...
    "E", ...
    "W", ...
    "G"];

maeThresholds = [ ...
    0.025, ...
    0.050, ...
    0.075, ...
    0.100];

CS = result.CurveSummary;
FP = result.CurveFitParams;
BM = result.BestModelPerCurve;

if isempty(CS)

    T = emptyRestorationShapeSummaryTable();
    return;
end

mainMask = ...
    CS.n_points >= opts.mainMinPoints & ...
    ~startsWith(string(CS.bestModel), "skipped:");

CS = CS(mainMask, :);

rows = cell(0, 1);

for g = 1:numel(groups)

    if groups(g) == "All systems"

        ids = unique( ...
            string(CS.record_ID), ...
            "stable");

    else

        systemMask = ...
            canonicalizeSystemVector(CS.System) == groupCodes(g);

        ids = unique( ...
            string(CS.record_ID(systemMask)), ...
            "stable");
    end

    nEligible = numel(ids);

    B = BM( ...
        ismember(string(BM.record_ID), ids), ...
        :);

    validWinner = ismember( ...
        string(B.bestModel), ...
        models);

    for m = 1:numel(models)

        mdl = models(m);

        F = FP( ...
            ismember(string(FP.record_ID), ids) & ...
            string(FP.model) == mdl, ...
            :);

        finiteMae = isfinite(F.MAE);

        winnerCount = nnz( ...
            validWinner & ...
            string(B.bestModel) == mdl);

        rows{end + 1, 1} = struct( ...
            "group", groups(g), ...
            "system_code", groupCodes(g), ...
            "model", mdl, ...
            "nEligible", nEligible, ...
            "winnerRate", ...
                safeDivide(winnerCount, nEligible), ...
            "shareMAE_le_0p025", ...
                safeDivide( ...
                    nnz(finiteMae & ...
                    F.MAE <= maeThresholds(1)), ...
                    nEligible), ...
            "shareMAE_le_0p05", ...
                safeDivide( ...
                    nnz(finiteMae & ...
                    F.MAE <= maeThresholds(2)), ...
                    nEligible), ...
            "shareMAE_le_0p075", ...
                safeDivide( ...
                    nnz(finiteMae & ...
                    F.MAE <= maeThresholds(3)), ...
                    nEligible), ...
            "shareMAE_le_0p10", ...
                safeDivide( ...
                    nnz(finiteMae & ...
                    F.MAE <= maeThresholds(4)), ...
                    nEligible));
    end
end

if isempty(rows)

    T = emptyRestorationShapeSummaryTable();

else

    T = struct2table(vertcat(rows{:}));

    T = T(:, { ...
        'group', ...
        'system_code', ...
        'model', ...
        'nEligible', ...
        'winnerRate', ...
        'shareMAE_le_0p025', ...
        'shareMAE_le_0p05', ...
        'shareMAE_le_0p075', ...
        'shareMAE_le_0p10'});
end

end

function T = emptyRestorationShapeSummaryTable()

T = table( ...
    strings(0, 1), ...
    strings(0, 1), ...
    strings(0, 1), ...
    zeros(0, 1), ...
    nan(0, 1), ...
    nan(0, 1), ...
    nan(0, 1), ...
    nan(0, 1), ...
    nan(0, 1), ...
    'VariableNames', { ...
        'group', ...
        'system_code', ...
        'model', ...
        'nEligible', ...
        'winnerRate', ...
        'shareMAE_le_0p025', ...
        'shareMAE_le_0p05', ...
        'shareMAE_le_0p075', ...
        'shareMAE_le_0p10'});

end

function writeRestorationShapeSourceData(summary, filename)

SourceData = buildRestorationShapeSourceData(summary);
CoreWinnerRate = buildCoreWinnerRateSummary(SourceData);

if isfile(filename)
    delete(filename);
end

writetable(SourceData, filename, "Sheet", "Fig3a_source_data");
writetable(CoreWinnerRate, filename, "Sheet", "Core_model_winner_rate");

if ~isfile(filename)
    error("Restoration-shape source-data workbook was not created:\n%s", filename);
end

end

function T = buildRestorationShapeSourceData(summary)

summary = summary(string(summary.system_code) ~= "ALL", :);

systems = ["Power", "Water", "Gas"];
systemLabels = ["Electric power", "Water supply", "Natural gas"];
models = [ ...
    "Weibull", ...
    "Gompertz", ...
    "Log-logistic", ...
    "Exponential", ...
    "Mechanical analogy", ...
    "Opabola", ...
    "Zorn", ...
    "Trigonometric"];

rows = cell(0, 1);

for s = 1:numel(systems)
    S = summary(string(summary.group) == systems(s), :);

    for m = 1:numel(models)
        R = S(string(S.model) == models(m), :);

        if isempty(R)
            continue;
        end

        rows{end + 1, 1} = struct( ...
            "system", systemLabels(s), ...
            "model", models(m), ...
            "nEligible", double(R.nEligible(1)), ...
            "winnerRate", double(R.winnerRate(1)), ...
            "shareMAE_le_0p05", double(R.shareMAE_le_0p05(1)));
    end
end

if isempty(rows)
    T = table( ...
        strings(0, 1), ...
        strings(0, 1), ...
        zeros(0, 1), ...
        nan(0, 1), ...
        nan(0, 1), ...
        "VariableNames", { ...
            'system', ...
            'model', ...
            'nEligible', ...
            'winnerRate', ...
            'shareMAE_le_0p05'});
else
    T = struct2table(vertcat(rows{:}));
    T = T(:, { ...
        'system', ...
        'model', ...
        'nEligible', ...
        'winnerRate', ...
        'shareMAE_le_0p05'});
end

end

function T = buildCoreWinnerRateSummary(SourceData)

systems = ["Electric power", "Water supply", "Natural gas"];
coreModels = ["Weibull", "Gompertz", "Log-logistic"];
rows = cell(0, 1);

for s = 1:numel(systems)
    S = SourceData(string(SourceData.system) == systems(s), :);

    if isempty(S)
        continue;
    end

    rates = nan(1, numel(coreModels));

    for m = 1:numel(coreModels)
        R = S(string(S.model) == coreModels(m), :);
        if ~isempty(R)
            rates(m) = double(R.winnerRate(1));
        end
    end

    rows{end + 1, 1} = struct( ...
        "system", systems(s), ...
        "nEligible", double(S.nEligible(1)), ...
        "Weibull_winnerRate", rates(1), ...
        "Gompertz_winnerRate", rates(2), ...
        "LogLogistic_winnerRate", rates(3), ...
        "coreModelWinnerRateSum", sum(rates, "omitnan"));
end

if isempty(rows)
    T = table( ...
        strings(0, 1), ...
        zeros(0, 1), ...
        nan(0, 1), ...
        nan(0, 1), ...
        nan(0, 1), ...
        nan(0, 1), ...
        "VariableNames", { ...
            'system', ...
            'nEligible', ...
            'Weibull_winnerRate', ...
            'Gompertz_winnerRate', ...
            'LogLogistic_winnerRate', ...
            'coreModelWinnerRateSum'});
else
    T = struct2table(vertcat(rows{:}));
    T = T(:, { ...
        'system', ...
        'nEligible', ...
        'Weibull_winnerRate', ...
        'Gompertz_winnerRate', ...
        'LogLogistic_winnerRate', ...
        'coreModelWinnerRateSum'});
end

end

function T = buildBadFitWeibullGompertzLoglogisticCandidateTable( ...
        result, ...
        curves, ...
        opts)

CS = result.CurveSummary;
FP = result.CurveFitParams;

if isempty(CS)

    T = emptyBadFitCandidateTable();
    return;
end

mask = ...
    CS.n_points >= opts.mainMinPoints & ...
    logical(CS.badFitWeibullGompertzLoglogistic);

B = CS(mask, :);

rows = cell(0, 1);

models = badFitModelNames();

modelKeys = [ ...
    "weibull", ...
    "gompertz", ...
    "loglogistic"];

for i = 1:height(B)

    recordID = string(B.record_ID(i));

    recIdx = find( ...
        string({curves.record_ID}) == recordID, ...
        1);

    if isempty(recIdx)
        continue;
    end

    rec = curves(recIdx);

    fitByModel = cell(1, numel(models));
    maeValues = nan(1, numel(models));
    convergedValues = false(1, numel(models));

    for m = 1:numel(models)

        F = FP( ...
            string(FP.record_ID) == recordID & ...
            string(FP.model) == models(m), ...
            :);

        if isempty(F)

            fitByModel{m} = [];

        else

            fitByModel{m} = F(1, :);
            maeValues(m) = F.MAE(1);
            convergedValues(m) = logical(F.converged(1));
        end
    end

    if ~( ...
            all(convergedValues) && ...
            all(isfinite(maeValues)) && ...
            all(maeValues > opts.maeThreshold))

        continue;
    end

    [bestMAE, bestIdx] = min(maeValues);

    row = struct();

    row.record_ID = recordID;

    row.exceptionClass = "";
    row.exceptionClassLabel = "";
    row.classNote = "";

    row.event_ID = getNumericField(rec, "event_ID");
    row.regionIndex = getNumericField(rec, "regionIndex");
    row.regionName = getStringField(rec, "regionName");
    row.System = getStringField(rec, "System");
    row.Country = getStringField(rec, "Country");
    row.magnitude = getNumericField(rec, "magnitude");

    row.badfitCriterion = ...
        "Weibull, Gompertz and Log-logistic all converged and " + ...
        "all have finite MAE > " + ...
        string(opts.maeThreshold);

    for m = 1:numel(models)

        key = modelKeys(m);
        F = fitByModel{m};

        row.("MAE_" + key) = ...
            tableScalar(F, "MAE");

        row.("RMSE_" + key) = ...
            tableScalar(F, "RMSE");

        row.("R2_" + key) = ...
            tableScalar(F, "R2");

    end

    row.bestMAE_Weibull_Gompertz_Loglogistic = bestMAE;
    row.bestModel_Weibull_Gompertz_Loglogistic = models(bestIdx);

    rows{end + 1, 1} = row;
end

if isempty(rows)

    T = emptyBadFitCandidateTable();

else

    T = struct2table(vertcat(rows{:}));
end

end

function T = emptyBadFitCandidateTable()

row = struct();

row.record_ID = "";

row.exceptionClass = "";
row.exceptionClassLabel = "";
row.classNote = "";

row.event_ID = NaN;
row.regionIndex = NaN;
row.regionName = "";
row.System = "";
row.Country = "";
row.magnitude = NaN;

row.badfitCriterion = "";

modelKeys = [ ...
    "weibull", ...
    "gompertz", ...
    "loglogistic"];

for m = 1:numel(modelKeys)

    key = modelKeys(m);

    row.("MAE_" + key) = NaN;
    row.("RMSE_" + key) = NaN;
    row.("R2_" + key) = NaN;
end

row.bestMAE_Weibull_Gompertz_Loglogistic = NaN;
row.bestModel_Weibull_Gompertz_Loglogistic = "";

T = struct2table(row);
T(1, :) = [];

end

function x = tableScalar(T, variableName)

x = NaN;

if isempty(T) || ...
        ~ismember( ...
        string(variableName), ...
        string(T.Properties.VariableNames))

    return;
end

v = T.(char(variableName));

if ~isempty(v)
    x = v(1);
end

end

function y = safeDivide(a, b)

if b > 0
    y = a / b;
else
    y = NaN;
end

end

function codes = canonicalizeSystemVector(systemNames)

systemNames = string(systemNames);
codes = strings(size(systemNames));

for i = 1:numel(systemNames)
    codes(i) = canonicalSystemCode(systemNames(i));
end

end


function fit = fitOneModelR(t, r, modelName, opts, tau100)

fit = struct();

fit.model = string(modelName);
fit.maskFit = true(size(t));
fit.converged = false;
fit.atBoundary = false;
fit.boundaryParameter = "";
fit.boundarySide = "";
fit.boundaryLimit = "";
fit.boundaryScaledDistance = NaN;
fit.boundaryTolerance = 1e-6;
fit.m = 2;

switch string(modelName)

    case "Exponential"

        fit.m = 1;

        if numel(t) < opts.minFitPoints
            fit = makeInvalidFit(fit, t);
            return;
        end

        lb = opts.lambdaBounds(1);
        ub = opts.lambdaBounds(2);
        paramNames = "lambda";

        lambda0 = initLambdaFromR(t, r);

        p0 = lambda0 .* opts.initPositiveMultipliers(:);

        best = runMultistart( ...
            @(p, tt) modelExponential(tt, p(1)), ...
            t, ...
            r, ...
            p0, ...
            lb, ...
            ub, ...
            opts.nRestarts);

        fit.params = defaultParams();
        fit.params.lambda = best.p(1);

        fit.rhat = modelExponential(t, best.p(1));
        fit.SSE = best.SSE;
        fit.converged = best.exitflag > 0;
        [ ...
            fit.atBoundary, ...
            fit.boundaryParameter, ...
            fit.boundarySide, ...
            fit.boundaryLimit, ...
            fit.boundaryScaledDistance, ...
            fit.boundaryTolerance] = getBoundaryDiagnostics( ...
                best.p, lb, ub, paramNames);

    case "Weibull"

        if numel(t) < opts.minFitPoints
            fit = makeInvalidFit(fit, t);
            return;
        end

        lb = [ ...
            opts.lambdaBounds(1), ...
            opts.kBounds(1)];

        ub = [ ...
            opts.lambdaBounds(2), ...
            opts.kBounds(2)];
        paramNames = ["lambda", "k"];

        lambda0 = initLambdaFromR(t, r);

        p0 = makeCartesianStarts( ...
            lambda0 .* opts.initPositiveMultipliers, ...
            opts.initPositiveMultipliers);

        best = runMultistart( ...
            @(p, tt) modelWeibull(tt, p), ...
            t, ...
            r, ...
            p0, ...
            lb, ...
            ub, ...
            opts.nRestarts);

        fit.params = defaultParams();

        fit.params.lambda = best.p(1);
        fit.params.k = best.p(2);

        fit.rhat = modelWeibull(t, best.p);
        fit.SSE = best.SSE;
        fit.converged = best.exitflag > 0;
        [ ...
            fit.atBoundary, ...
            fit.boundaryParameter, ...
            fit.boundarySide, ...
            fit.boundaryLimit, ...
            fit.boundaryScaledDistance, ...
            fit.boundaryTolerance] = getBoundaryDiagnostics( ...
                best.p, lb, ub, paramNames);

    case "Gompertz"

        if numel(t) < opts.minFitPoints
            fit = makeInvalidFit(fit, t);
            return;
        end

        lb = [ ...
            opts.bBounds(1), ...
            opts.cBounds(1)];

        ub = [ ...
            opts.bBounds(2), ...
            opts.cBounds(2)];
        paramNames = ["b", "c"];

        tau50Init = initT50FromR(t, r);
        b0 = 1 / tau50Init;
        c0 = log(2) / expm1(1);

        p0 = makeCartesianStarts( ...
            b0 .* opts.initPositiveMultipliers, ...
            c0 .* opts.initPositiveMultipliers);

        best = runMultistart( ...
            @(p, tt) modelGompertz(tt, p), ...
            t, ...
            r, ...
            p0, ...
            lb, ...
            ub, ...
            opts.nRestarts);

        fit.params = defaultParams();

        fit.params.b = best.p(1);
        fit.params.c = best.p(2);

        fit.rhat = modelGompertz(t, best.p);
        fit.SSE = best.SSE;
        fit.converged = best.exitflag > 0;
        [ ...
            fit.atBoundary, ...
            fit.boundaryParameter, ...
            fit.boundarySide, ...
            fit.boundaryLimit, ...
            fit.boundaryScaledDistance, ...
            fit.boundaryTolerance] = getBoundaryDiagnostics( ...
                best.p, lb, ub, paramNames);

    case "Log-logistic"

        if numel(t) < opts.minFitPoints
            fit = makeInvalidFit(fit, t);
            return;
        end

        lb = [ ...
            opts.lambdaBounds(1), ...
            opts.kBounds(1)];

        ub = [ ...
            opts.lambdaBounds(2), ...
            opts.kBounds(2)];
        paramNames = ["lambda", "k"];


        lambda0 = initT50FromR(t, r);

        p0 = makeCartesianStarts( ...
            lambda0 .* opts.initPositiveMultipliers, ...
            opts.initPositiveMultipliers);

        best = runMultistart( ...
            @(p, tt) modelLogLogistic(tt, p), ...
            t, ...
            r, ...
            p0, ...
            lb, ...
            ub, ...
            opts.nRestarts);

        fit.params = defaultParams();

        fit.params.lambda = best.p(1);
        fit.params.k = best.p(2);

        fit.rhat = modelLogLogistic(t, best.p);
        fit.SSE = best.SSE;
        fit.converged = best.exitflag > 0;
        [ ...
            fit.atBoundary, ...
            fit.boundaryParameter, ...
            fit.boundarySide, ...
            fit.boundaryLimit, ...
            fit.boundaryScaledDistance, ...
            fit.boundaryTolerance] = getBoundaryDiagnostics( ...
                best.p, lb, ub, paramNames);

    case "Mechanical analogy"

        fit.m = 1;

        if numel(t) < opts.minFitPoints
            fit = makeInvalidFit(fit, t);
            return;
        end

        lb = opts.wBounds(1);
        ub = opts.wBounds(2);
        paramNames = "w";


        tau50Init = initT50FromR(t, r);
        x50Mechanical = 1.67834699001666;
        w0 = x50Mechanical / tau50Init;

        p0 = w0 .* opts.initPositiveMultipliers(:);

        best = runMultistart( ...
            @(p, tt) modelMechanicalAnalogy(tt, p(1)), ...
            t, ...
            r, ...
            p0, ...
            lb, ...
            ub, ...
            opts.nRestarts);

        fit.params = defaultParams();
        fit.params.w = best.p(1);

        fit.rhat = modelMechanicalAnalogy(t, best.p(1));
        fit.SSE = best.SSE;
        fit.converged = best.exitflag > 0;
        [ ...
            fit.atBoundary, ...
            fit.boundaryParameter, ...
            fit.boundarySide, ...
            fit.boundaryLimit, ...
            fit.boundaryScaledDistance, ...
            fit.boundaryTolerance] = getBoundaryDiagnostics( ...
                best.p, lb, ub, paramNames);

    case "Opabola"

        fit.m = 1;

        if numel(t) < opts.minFitPoints
            fit = makeInvalidFit(fit, t);
            return;
        end

        requireTau100(tau100, modelName);

        lb = opts.gammaBounds(1);
        ub = opts.gammaBounds(2);
        paramNames = "gamma";


        p0 = opts.initPositiveMultipliers(:) - 1;

        best = runMultistart( ...
            @(p, tt) modelOpabola(tt, p(1), tau100), ...
            t, ...
            r, ...
            p0, ...
            lb, ...
            ub, ...
            opts.nRestarts);

        fit.params = defaultParams();

        fit.params.gamma = best.p(1);
        fit.params.tau100 = tau100;

        fit.rhat = modelOpabola(t, best.p(1), tau100);
        fit.SSE = best.SSE;
        fit.converged = best.exitflag > 0;
        [ ...
            fit.atBoundary, ...
            fit.boundaryParameter, ...
            fit.boundarySide, ...
            fit.boundaryLimit, ...
            fit.boundaryScaledDistance, ...
            fit.boundaryTolerance] = getBoundaryDiagnostics( ...
                best.p, lb, ub, paramNames);

    case "Zorn"

        if numel(t) < opts.minFitPoints
            fit = makeInvalidFit(fit, t);
            return;
        end

        requireTau100(tau100, modelName);

        lb = [-Inf, -Inf];
        ub = [ Inf,  Inf];


        e0 = opts.initPositiveMultipliers;
        rho0 = opts.initUnitIntervalFractions;
        [eGrid, rhoGrid] = ndgrid(e0, rho0);

        q0 = [ ...
            log(eGrid(:)), ...
            log(rhoGrid(:) ./ (1 - rhoGrid(:)))];

        best = runMultistart( ...
            @(q, tt) modelZornConstrainedQ(tt, q, tau100), ...
            t, ...
            r, ...
            q0, ...
            lb, ...
            ub, ...
            opts.nRestarts);

        [eBest, fBest] = zornConstrainedParams(best.p);

        fit.params = defaultParams();

        fit.params.e = eBest;
        fit.params.f = fBest;
        fit.params.tau100 = tau100;

        fit.rhat = modelZorn( ...
            t, ...
            [eBest, fBest], ...
            tau100);

        fit.SSE = best.SSE;
        fit.converged = best.exitflag > 0;
        [ ...
            fit.atBoundary, ...
            fit.boundaryParameter, ...
            fit.boundarySide, ...
            fit.boundaryLimit, ...
            fit.boundaryScaledDistance, ...
            fit.boundaryTolerance] = getZornBoundaryDiagnostics( ...
                eBest, ...
                fBest);

    case "Trigonometric"

        fit.m = 1;

        if numel(t) < opts.minFitPoints
            fit = makeInvalidFit(fit, t);
            return;
        end

        requireTau100(tau100, modelName);

        lb = opts.aBounds(1);
        ub = opts.aBounds(2);
        paramNames = "a";


        p0 = opts.initUnitIntervalFractions(:);

        best = runMultistart( ...
            @(p, tt) modelTrigonometric(tt, p(1), tau100), ...
            t, ...
            r, ...
            p0, ...
            lb, ...
            ub, ...
            opts.nRestarts);

        fit.params = defaultParams();

        fit.params.a = best.p(1);
        fit.params.tau100 = tau100;

        fit.rhat = modelTrigonometric(t, best.p(1), tau100);
        fit.SSE = best.SSE;
        fit.converged = best.exitflag > 0;
        [ ...
            fit.atBoundary, ...
            fit.boundaryParameter, ...
            fit.boundarySide, ...
            fit.boundaryLimit, ...
            fit.boundaryScaledDistance, ...
            fit.boundaryTolerance] = getBoundaryDiagnostics( ...
                best.p, lb, ub, paramNames);

    otherwise

        error( ...
            "Unknown model: %s", ...
            modelName);
end

end

function fit = makeInvalidFit(fit, t)

fit.params = defaultParams();
fit.rhat = nan(size(t));
fit.SSE = NaN;
fit.converged = false;
fit.atBoundary = false;
fit.boundaryParameter = "";
fit.boundarySide = "";
fit.boundaryLimit = "";
fit.boundaryScaledDistance = NaN;
fit.boundaryTolerance = 1e-6;

end

function p = defaultParams()

p = struct( ...
    "lambda", NaN, ...
    "k", NaN, ...
    "b", NaN, ...
    "c", NaN, ...
    "w", NaN, ...
    "gamma", NaN, ...
    "e", NaN, ...
    "f", NaN, ...
    "a", NaN, ...
    "tau100", NaN);

end

function y = modelExponential(t, lambda)

y = 1 - exp(-t ./ lambda);
y = min(max(y, 0), 1);

end

function y = modelWeibull(t, p)

lambda = p(1);
k = p(2);

y = 1 - exp(-(t ./ lambda) .^ k);
y = min(max(y, 0), 1);

end

function y = modelGompertz(t, p)

b = p(1);
c = p(2);

y = 1 - exp(-c * (exp(b * t) - 1));
y = min(max(y, 0), 1);

end

function y = modelLogLogistic(t, p)

lambda = p(1);
k = p(2);

y = zeros(size(t));

mask = t > 0;

z = k .* ( ...
    log(t(mask)) - ...
    log(lambda));

y(mask) = 1 ./ (1 + exp(-z));
y = min(max(y, 0), 1);

end

function y = modelMechanicalAnalogy(t, w)

y = 1 - exp(-w .* t) .* (1 + w .* t);
y = min(max(y, 0), 1);

end

function y = modelOpabola(t, gamma, tau100)

x = t ./ tau100;

y = ...
    (x .* (1 + gamma)) ./ ...
    (1 + x .* gamma);

y = min(max(y, 0), 1);

end

function y = modelZorn(t, p, tau100)

e = p(1);
f = p(2);

x = t ./ tau100;

y = ...
    (x .* (1 + e)) ./ ...
    (x .^ f + e);

y = min(max(y, 0), 1);

end

function y = modelZornConstrainedQ(t, q, tau100)

[e, f] = zornConstrainedParams(q);

y = modelZorn( ...
    t, ...
    [e, f], ...
    tau100);

end

function [e, f] = zornConstrainedParams(q)

e = exp(q(1));

s = 1 / (1 + exp(-q(2)));

f = (e + 1) * s;

end

function y = modelTrigonometric(t, a, tau100)

x = t ./ tau100;

y = zeros(size(t));

maskIn = ...
    x >= 0 & ...
    x <= 1;

if abs(a) < sqrt(eps)

    y(maskIn) = x(maskIn) .^ 2;

else

    den = sin(pi * a / 2);

    y(maskIn) = ...
        (sin(pi * a .* x(maskIn) / 2) ./ den) .^ 2;
end

y(x > 1) = 1;
y(x < 0) = 0;

y = min(max(y, 0), 1);

end

function modelOrder = buildModelOrderForCurve(modelOrder, hasTau100)

tau100Models = [ ...
    "Opabola", ...
    "Zorn", ...
    "Trigonometric"];

if ~hasTau100
    modelOrder = modelOrder( ...
        ~ismember(modelOrder, tau100Models));
end

end

function requireTau100(tau100, modelName)

if ~(isfinite(tau100) && tau100 > 0)

    error( ...
        "Model %s requires a finite positive tau100.", ...
        string(modelName));
end

end

function lambda0 = initLambdaFromR(t, r)

target = 1 - exp(-1);
lambda0 = initTAtR(t, r, target);

if isnan(lambda0)
    lambda0 = fallbackPositiveTimeScale(t);
end

end

function tau50Init = initT50FromR(t, r)

tau50Init = initTAtR(t, r, 0.5);

if isnan(tau50Init)
    tau50Init = fallbackPositiveTimeScale(t);
end

end

function t0 = fallbackPositiveTimeScale(t)

tPositive = t(isfinite(t) & t > 0);

if isempty(tPositive)
    t0 = 1;
else
    t0 = max(max(tPositive) / 2, 1e-3);
end

end

function p0 = makeCartesianStarts(v1, v2)


v1 = double(v1(:));
v2 = double(v2(:));

[g1, g2] = ndgrid(v1, v2);
p0 = [g1(:), g2(:)];

end

function tEstimate = initTAtR(t, r, target)

tEstimate = NaN;

mask = isfinite(t) & isfinite(r);

t = t(mask);
r = r(mask);

if numel(t) < 2
    return;
end

[t, order] = sort(t);
r = r(order);

idx = find(r >= target, 1, "first");

if isempty(idx)
    return;
end

if idx == 1

    tEstimate = t(1);

else

    tEstimate = interp1( ...
        r(idx - 1:idx), ...
        t(idx - 1:idx), ...
        target, ...
        "linear", ...
        "extrap");
end

if ~isfinite(tEstimate) || tEstimate <= 0
    tEstimate = NaN;
end

end


function best = runMultistart( ...
        modelFun, ...
        t, ...
        r, ...
        p0, ...
        lb, ...
        ub, ...
        nRestarts)

bestConverged = struct( ...
    "p", [], ...
    "SSE", Inf, ...
    "exitflag", -1);

bestAny = struct( ...
    "p", [], ...
    "SSE", Inf, ...
    "exitflag", -1);

if exist("lsqcurvefit", "file") ~= 2
    error( ...
        "RestorationFit:MissingOptimizer", ...
        [ ...
        "This analysis requires lsqcurvefit from MATLAB Optimization Toolbox. " + ...
        "No fallback optimizer is used." ...
        ]);
end


opt = optimoptions( ...
    "lsqcurvefit", ...
    "Display", "off", ...
    "MaxIterations", 2000, ...
    "MaxFunctionEvaluations", 10000, ...
    "FunctionTolerance", 1e-10, ...
    "StepTolerance", 1e-10, ...
    "OptimalityTolerance", 1e-8);

p0 = double(p0);

if isscalar(lb)

    p0 = p0(:);

elseif isvector(p0)

    p0 = p0(:)';
end


nBaseStarts = size(p0, 1);
nRestarts = max(nBaseStarts, max(1, round(nRestarts)));

nRuntimeFailures = 0;
firstRuntimeError = [];

for restartIndex = 1:nRestarts

    baseIndex = mod(restartIndex - 1, nBaseStarts) + 1;
    baseInit = p0(baseIndex, :);

    if restartIndex <= nBaseStarts
        pInit = baseInit;
    else
        pInit = jitterInit( ...
            baseInit, ...
            lb, ...
            ub, ...
            restartIndex);
    end

    try

        [p, ~, residual, exitflag] = lsqcurvefit( ...
            modelFun, ...
            pInit, ...
            t, ...
            r, ...
            lb, ...
            ub, ...
            opt);

        if any(~isfinite(residual))
            SSE = Inf;
        else
            SSE = sum(residual .^ 2);
        end

        if isfinite(SSE) && SSE < bestAny.SSE

            bestAny = struct( ...
                "p", p, ...
                "SSE", SSE, ...
                "exitflag", exitflag);
        end

        if exitflag > 0 && ...
                isfinite(SSE) && ...
                SSE < bestConverged.SSE

            bestConverged = struct( ...
                "p", p, ...
                "SSE", SSE, ...
                "exitflag", exitflag);
        end

    catch ME

        nRuntimeFailures = nRuntimeFailures + 1;

        if isempty(firstRuntimeError)
            firstRuntimeError = ME;
        end
    end
end

if nRuntimeFailures == nRestarts

    error( ...
        "RestorationFit:AllRestartsFailed", ...
        [ ...
        "All %d multistart attempts raised runtime errors. " + ...
        "First error [%s]: %s" ...
        ], ...
        nRestarts, ...
        string(firstRuntimeError.identifier), ...
        string(firstRuntimeError.message));

elseif nRuntimeFailures > 0

    warning( ...
        "RestorationFit:PartialRestartFailures", ...
        [ ...
        "%d of %d multistart attempts raised runtime errors and were skipped. " + ...
        "First error [%s]: %s" ...
        ], ...
        nRuntimeFailures, ...
        nRestarts, ...
        string(firstRuntimeError.identifier), ...
        string(firstRuntimeError.message));
end

if ~isempty(bestConverged.p)

    best = bestConverged;

elseif ~isempty(bestAny.p)

    best = bestAny;

else

    best.p = p0(1, :);

    rh = modelFun(best.p, t);

    if any(~isfinite(rh))
        best.SSE = Inf;
    else
        best.SSE = sum((r - rh) .^ 2);
    end

    best.exitflag = -1;
end

end

function pInit = jitterInit(p0, lb, ub, seed)


boundEps = 1e-12;

rng(1000 + seed, "twister");

pInit = p0;

for i = 1:numel(p0)

    if isfinite(lb(i)) && lb(i) >= 0

        scale = 10 ^ (0.45 * randn);

        pInit(i) = ...
            max(p0(i), max(lb(i), boundEps)) * scale;

    else

        pInit(i) = ...
            p0(i) + 0.5 * randn;
    end

    if isfinite(lb(i))
        pInit(i) = max( ...
            pInit(i), ...
            lb(i) + boundEps);
    end

    if isfinite(ub(i))
        pInit(i) = min( ...
            pInit(i), ...
            ub(i) - boundEps);
    end
end

end

function [tf, parameterText, sideText, limitText, minScaledDistance, boundaryTol] = ...
        getBoundaryDiagnostics(p, lb, ub, parameterNames)


boundaryTol = 1e-6;

p = double(p(:));
lb = double(lb(:));
ub = double(ub(:));
parameterNames = string(parameterNames(:));

if numel(parameterNames) ~= numel(p)
    error( ...
        "RestorationFit:BoundaryNameMismatch", ...
        "Number of parameter names must match number of fitted parameters.");
end

hitParameters = strings(0, 1);
hitSides = strings(0, 1);
hitLimits = strings(0, 1);
minScaledDistance = Inf;

for i = 1:numel(p)

    scale = max(1, abs(p(i)));

    if isfinite(lb(i))
        dLower = abs(p(i) - lb(i)) / scale;
        minScaledDistance = min(minScaledDistance, dLower);

        if dLower <= boundaryTol
            hitParameters(end + 1, 1) = parameterNames(i);
            hitSides(end + 1, 1) = "lower";
            hitLimits(end + 1, 1) = string(sprintf("%.15g", lb(i)));
        end
    end

    if isfinite(ub(i))
        dUpper = abs(ub(i) - p(i)) / scale;
        minScaledDistance = min(minScaledDistance, dUpper);

        if dUpper <= boundaryTol
            hitParameters(end + 1, 1) = parameterNames(i);
            hitSides(end + 1, 1) = "upper";
            hitLimits(end + 1, 1) = string(sprintf("%.15g", ub(i)));
        end
    end
end

if isinf(minScaledDistance)
    minScaledDistance = NaN;
end

tf = ~isempty(hitParameters);
parameterText = strjoin(hitParameters, ";");
sideText = strjoin(hitSides, ";");
limitText = strjoin(hitLimits, ";");

end

function [tf, parameterText, sideText, limitText, minScaledDistance, boundaryTol] = ...
        getZornBoundaryDiagnostics(e, f)

boundaryTol = 1e-6;

e = double(e);
f = double(f);

distances = [e, f, e + 1 - f];
[minScaledDistance, idx] = min(distances ./ max(1, abs([e, f, e + 1])));

if ~isfinite(minScaledDistance)
    minScaledDistance = NaN;
end

tf = isfinite(minScaledDistance) && minScaledDistance <= boundaryTol;
parameterText = "";
sideText = "";
limitText = "";

if tf
    switch idx
        case 1
            parameterText = "e";
            sideText = "lower";
            limitText = "0";
        case 2
            parameterText = "f";
            sideText = "lower";
            limitText = "0";
        case 3
            parameterText = "f";
            sideText = "upper";
            limitText = "e+1";
    end
end

end


function met = computeMetrics(y, yhat, SSE, m, maskFit)

if nargin < 5 || isempty(maskFit)
    maskFit = true(size(y));
end

res = ...
    y(maskFit) - ...
    yhat(maskFit);

res = res(isfinite(res));

n = numel(res);

met = struct();

if n == 0

    met.R2 = NaN;
    met.RMSE = NaN;
    met.MAE = NaN;
    met.MaxAbsError = NaN;
    met.SSE = NaN;
    met.AICc = NaN;
    met.BIC = NaN;

    return;
end

SSEuse = sumNoNan(res .^ 2);

if nargin >= 3 && isfinite(SSE)
    SSEuse = SSE;
end

yy = y(maskFit);

SST = sumNoNan( ...
    (yy - meanNoNan(yy)) .^ 2);

met.SSE = SSEuse;

if isfinite(SSEuse) && SST > 0
    met.R2 = 1 - SSEuse / SST;
else
    met.R2 = NaN;
end

met.RMSE = sqrt(meanNoNan(res .^ 2));
met.MAE = meanNoNan(abs(res));
met.MaxAbsError = max(abs(res));

K = m + 1;

if n > K + 1 && ...
        isfinite(SSEuse) && ...
        SSEuse >= 0

    SSEforIC = max(SSEuse, realmin("double"));

    AIC = ...
        n * log(SSEforIC / n) + ...
        2 * K;

    met.BIC = ...
        n * log(SSEforIC / n) + ...
        K * log(n);

else

    AIC = NaN;
    met.BIC = NaN;
end

if n - K - 1 > 0

    met.AICc = ...
        AIC + ...
        (2 * K * (K + 1)) / ...
        (n - K - 1);

else

    met.AICc = NaN;
end

end


function fitRow = makeFitRow( ...
        rec, ...
        datasetLabel, ...
        modelName, ...
        fit, ...
        met, ...
        peakTime)

p = fit.params;
recordTau100 = getNumericField(rec, "tau100");

fitRow = struct( ...
    "dataset", string(datasetLabel), ...
    "record_ID", string(rec.record_ID), ...
    "event_ID", getNumericField(rec, "event_ID"), ...
    "regionIndex", getNumericField(rec, "regionIndex"), ...
    "regionName", getStringField(rec, "regionName"), ...
    "System", getStringField(rec, "System"), ...
    "Country", getStringField(rec, "Country"), ...
    "magnitude", getNumericField(rec, "magnitude"), ...
    "PeakTime", peakTime, ...
    "n_points", getNumericField( ...
        rec, ...
        "nRestorationRatioPoints"), ...
    "sys", canonicalSystemCode( ...
        getStringField(rec, "System")), ...
    "model", string(modelName), ...
    "n_params", fit.m, ...
    "lambda", p.lambda, ...
    "k", p.k, ...
    "b", p.b, ...
    "c", p.c, ...
    "w", p.w, ...
    "gamma", p.gamma, ...
    "e", p.e, ...
    "f", p.f, ...
    "a", p.a, ...
    "tau100", recordTau100, ...
    "hasTau100", isfinite(recordTau100), ...
    "R2", met.R2, ...
    "RMSE", met.RMSE, ...
    "MAE", met.MAE, ...
    "MaxAbsError", met.MaxAbsError, ...
    "SSE", met.SSE, ...
    "AICc", met.AICc, ...
    "BIC", met.BIC, ...
    "converged", logical(fit.converged), ...
    "atBoundary", logical(fit.atBoundary), ...
    "boundaryParameter", string(fit.boundaryParameter), ...
    "boundarySide", string(fit.boundarySide), ...
    "boundaryLimit", string(fit.boundaryLimit), ...
    "boundaryScaledDistance", fit.boundaryScaledDistance, ...
    "boundaryTolerance", fit.boundaryTolerance);

end

function bestRow = chooseBestModelRow(fits)

converged = fits([fits.converged]);

if isempty(converged)

    bestRow = makeBestRowFromCandidate( ...
        makeInvalidBestCandidate(fits(1)));

    return;
end

valid = isfinite([converged.AICc]);

if any(valid)

    candidates = converged(valid);

    [~, order] = sortrows( ...
        [ ...
        [candidates.AICc]', ...
        [candidates.n_params]']);

    best = candidates(order(1));

else

    best = makeInvalidBestCandidate(converged(1));
end

bestRow = makeBestRowFromCandidate(best);

end

function bestRow = makeBestRowFromCandidate(best)

bestRow = struct( ...
    "dataset", best.dataset, ...
    "record_ID", best.record_ID, ...
    "event_ID", best.event_ID, ...
    "regionIndex", best.regionIndex, ...
    "regionName", best.regionName, ...
    "System", best.System, ...
    "Country", best.Country, ...
    "magnitude", best.magnitude, ...
    "PeakTime", best.PeakTime, ...
    "n_points", best.n_points, ...
    "sys", best.sys, ...
    "bestModel", best.model, ...
    "bestAICc", best.AICc, ...
    "bestBIC", best.BIC, ...
    "bestRMSE", best.RMSE, ...
    "bestMAE", best.MAE, ...
    "bestMaxAbsError", best.MaxAbsError, ...
    "bestR2", best.R2);

end

function best = makeInvalidBestCandidate(template)

best = template;

best.model = "invalid:no_finite_AICc";
best.n_params = NaN;
best.AICc = NaN;
best.BIC = NaN;
best.RMSE = NaN;
best.MAE = NaN;
best.MaxAbsError = NaN;
best.R2 = NaN;
best.converged = false;
best.atBoundary = false;

end

function summaryRow = makeSummaryRow( ...
        rec, ...
        datasetLabel, ...
        fits, ...
        bestRow, ...
        peakTime, ...
        maeThreshold)

badFitModels = badFitModelNames();

weib = pickFit(fits, badFitModels(1));
gomp = pickFit(fits, badFitModels(2));
llog = pickFit(fits, badFitModels(3));

isBadFit = isBadFitByMae( ...
    [ ...
    weib.MAE, ...
    gomp.MAE, ...
    llog.MAE], ...
    [ ...
    weib.converged, ...
    gomp.converged, ...
    llog.converged], ...
    maeThreshold);

summaryRow = struct( ...
    "dataset", string(datasetLabel), ...
    "record_ID", string(rec.record_ID), ...
    "event_ID", getNumericField(rec, "event_ID"), ...
    "regionIndex", getNumericField(rec, "regionIndex"), ...
    "regionName", getStringField(rec, "regionName"), ...
    "System", getStringField(rec, "System"), ...
    "Country", getStringField(rec, "Country"), ...
    "magnitude", getNumericField(rec, "magnitude"), ...
    "PeakTime", peakTime, ...
    "tau100", getNumericField(rec, "tau100"), ...
    "n_points", getNumericField( ...
        rec, ...
        "nRestorationRatioPoints"), ...
    "sys", canonicalSystemCode( ...
        getStringField(rec, "System")), ...
    "n_unique_values", getNumericField( ...
        rec, ...
        "nUniqueRestorationValues"), ...
    "badFitWeibullGompertzLoglogistic", logical(isBadFit), ...
    "bestModel", bestRow.bestModel, ...
    "bestMAE", bestRow.bestMAE, ...
    "weibull_MAE", weib.MAE, ...
    "gompertz_MAE", gomp.MAE, ...
    "loglogistic_MAE", llog.MAE, ...
    "weibull_R2", weib.R2, ...
    "gompertz_R2", gomp.R2, ...
    "loglogistic_R2", llog.R2);

end

function isBad = isBadFitByMae(maeValues, convergedValues, threshold)

isBad = ...
    all(logical(convergedValues)) && ...
    all(isfinite(maeValues)) && ...
    all(maeValues > threshold);

end

function row = makeSkippedSummary(rec, datasetLabel, reason)

row = struct( ...
    "dataset", string(datasetLabel), ...
    "record_ID", getStringField(rec, "record_ID"), ...
    "event_ID", getNumericField(rec, "event_ID"), ...
    "regionIndex", getNumericField(rec, "regionIndex"), ...
    "regionName", getStringField(rec, "regionName"), ...
    "System", getStringField(rec, "System"), ...
    "Country", getStringField(rec, "Country"), ...
    "magnitude", getNumericField(rec, "magnitude"), ...
    "PeakTime", getNumericField(rec, "PeakTime"), ...
    "tau100", getNumericField(rec, "tau100"), ...
    "n_points", getNumericField( ...
        rec, ...
        "nRestorationRatioPoints"), ...
    "sys", canonicalSystemCode( ...
        getStringField(rec, "System")), ...
    "n_unique_values", getNumericField( ...
        rec, ...
        "nUniqueRestorationValues"), ...
    "badFitWeibullGompertzLoglogistic", false, ...
    "bestModel", "skipped:" + string(reason), ...
    "bestMAE", NaN, ...
    "weibull_MAE", NaN, ...
    "gompertz_MAE", NaN, ...
    "loglogistic_MAE", NaN, ...
    "weibull_R2", NaN, ...
    "gompertz_R2", NaN, ...
    "loglogistic_R2", NaN);

end

function f = pickFit(fits, modelName)

idx = find( ...
    string({fits.model}) == string(modelName), ...
    1);

if isempty(idx)

    f = struct( ...
        "MAE", NaN, ...
        "R2", NaN, ...
        "converged", false);

else

    f = fits(idx);
end

end


function T = emptySummaryTable()

T = struct2table( ...
    makeSkippedSummary( ...
    emptyRec(), ...
    "", ...
    ""));

T(1, :) = [];

end

function T = emptyFitTable()

dummyFit = struct( ...
    "params", defaultParams(), ...
    "converged", false, ...
    "atBoundary", false, ...
    "boundaryParameter", "", ...
    "boundarySide", "", ...
    "boundaryLimit", "", ...
    "boundaryScaledDistance", NaN, ...
    "boundaryTolerance", 1e-6, ...
    "m", NaN);

met = struct( ...
    "R2", NaN, ...
    "RMSE", NaN, ...
    "MAE", NaN, ...
    "MaxAbsError", NaN, ...
    "SSE", NaN, ...
    "AICc", NaN, ...
    "BIC", NaN);

T = struct2table( ...
    makeFitRow( ...
    emptyRec(), ...
    "", ...
    "Weibull", ...
    dummyFit, ...
    met, ...
    NaN));

T(1, :) = [];

end

function T = emptyBestTable()

dummyFit = struct( ...
    "params", defaultParams(), ...
    "converged", false, ...
    "atBoundary", false, ...
    "boundaryParameter", "", ...
    "boundarySide", "", ...
    "boundaryLimit", "", ...
    "boundaryScaledDistance", NaN, ...
    "boundaryTolerance", 1e-6, ...
    "m", NaN);

met = struct( ...
    "R2", NaN, ...
    "RMSE", NaN, ...
    "MAE", NaN, ...
    "MaxAbsError", NaN, ...
    "SSE", NaN, ...
    "AICc", NaN, ...
    "BIC", NaN);

dummyFitRow = makeFitRow( ...
    emptyRec(), ...
    "", ...
    "Weibull", ...
    dummyFit, ...
    met, ...
    NaN);

T = struct2table( ...
    chooseBestModelRow(dummyFitRow));

T(1, :) = [];

end

function rec = emptyRec()

rec = struct( ...
    "record_ID", "", ...
    "event_ID", NaN, ...
    "regionIndex", NaN, ...
    "regionName", "", ...
    "System", "", ...
    "Country", "", ...
    "magnitude", NaN, ...
    "nRestorationRatioPoints", NaN, ...
    "nUniqueRestorationValues", NaN, ...
    "PeakTime", NaN, ...
    "tau100", NaN, ...
    "restorationratioCurve", []);

end

function models = badFitModelNames()

models = [ ...
    "Weibull", ...
    "Gompertz", ...
    "Log-logistic"];

end

function curves = loadPreparedRecords(moduleRoot)

emergentSimplicityRoot = string(fileparts(fileparts(moduleRoot)));
outputDir = fullfile( ...
    emergentSimplicityRoot, ...
    "outputs", ...
    "intermediate", ...
    "01_restoration_shape_phase_ratio");
inputMat = fullfile( ...
    outputDir, ...
    "unique3_records.mat");

if ~isfile(inputMat)
    error( ...
        "Prepared dataset not found:\n%s\n" + ...
        "Run run_00_prepare_restoration_shape_phase_ratio_data.m first.", ...
        inputMat);
end

S = load(inputMat, "unique3Records");

if ~isfield(S, "unique3Records")
    error( ...
        "MAT file does not contain unique3Records:\n%s", ...
        inputMat);
end

curves = S.unique3Records;

if isempty(curves)
    error( ...
        "unique3Records is empty:\n%s", ...
        inputMat);
end

requiredFields = [ ...
    "record_ID", ...
    "System", ...
    "restorationratioCurve", ...
    "nRestorationRatioPoints", ...
    "nUniqueRestorationValues", ...
    "PeakTime", ...
    "tau100"];

availableFields = string(fieldnames(curves));

missingFields = requiredFields( ...
    ~ismember(requiredFields, availableFields));

if ~isempty(missingFields)
    error( ...
        "unique3Records is missing required field(s): %s", ...
        strjoin(missingFields, ", "));
end

end

function curves = filterCurvesBySystem(curves, selectedSystems)

if isempty(curves)
    return;
end

systems = strings(numel(curves), 1);

for i = 1:numel(curves)
    systems(i) = normalizeSystemName( ...
        getStringField(curves(i), "System"));
end

keep = ismember(systems, selectedSystems);

curves = curves(keep);

end

function systems = normalizeSelectedSystems(systemsIn)

systemsIn = upper( ...
    strtrim(string(systemsIn(:)')));

systems = strings(size(systemsIn));

for i = 1:numel(systemsIn)

    switch systemsIn(i)

        case {"E", "POWER"}

            systems(i) = "Power";

        case {"W", "WATER"}

            systems(i) = "Water";

        case {"G", "GAS"}

            systems(i) = "Gas";

        otherwise

            error( ...
                "Unknown selected system: %s", ...
                systemsIn(i));
    end
end

systems = unique( ...
    systems, ...
    "stable");

end

function s = normalizeSystemName(s)

s = string(s);

sl = lower(strtrim(s));

if any(sl == ["e", "power"])

    s = "Power";

elseif any(sl == ["w", "water"])

    s = "Water";

elseif any(sl == ["g", "gas"])

    s = "Gas";
end

end

function sys = canonicalSystemCode(systemName)

s = upper(strtrim(string(systemName)));

sys = s;

if any(s == ["POWER", "E"])

    sys = "E";

elseif any(s == ["WATER", "W"])

    sys = "W";

elseif any(s == ["GAS", "G"])

    sys = "G";
end

end

function [t, r] = preparedCurveToArrays(C)

if isempty(C)
    t = [];
    r = [];
    return;
end

if isstruct(C)

    t = nan(numel(C), 1);
    r = nan(numel(C), 1);
    keep = true(numel(C), 1);

    for i = 1:numel(C)
        if isfield(C, "day_censor_code") && ...
                getNumericField(C(i), "day_censor_code") ~= 0

            keep(i) = false;
            continue;
        end

        t(i) = getNumericField(C(i), "day");
        r(i) = getNumericField(C(i), "value");
    end

    keep = keep & isfinite(t) & isfinite(r);
    t = t(keep);
    r = r(keep);

else

    C = double(C);

    if size(C, 2) < 2
        t = [];
        r = [];
        return;
    end

    t = C(:, 1);
    r = C(:, 2);

    keep = isfinite(t) & isfinite(r);
    t = t(keep);
    r = r(keep);
end

end

function x = getNumericField(S, fieldName)

x = NaN;

if ~isstruct(S) || ...
        ~isfield(S, fieldName) || ...
        isempty(S.(fieldName))

    return;
end

v = S.(fieldName);

if isnumeric(v) || islogical(v)

    if isscalar(v) && isfinite(double(v))
        x = double(v);
    end

elseif isstring(v) || ischar(v)

    y = str2double(string(v));

    if isscalar(y) && isfinite(y)
        x = y;
    end
end

end

function s = getStringField(S, fieldName)

s = "";

if ~isstruct(S) || ...
        ~isfield(S, fieldName) || ...
        isempty(S.(fieldName))

    return;
end

v = S.(fieldName);

if isstring(v)

    v(ismissing(v)) = "";

    s = string(v);

    if numel(s) > 1
        s = strjoin(s(:).', " ");
    end

elseif ischar(v)

    s = string(v);

elseif isnumeric(v) || islogical(v)

    if isscalar(v) && isfinite(double(v))
        s = string(v);
    end
end

end

function out = appendStruct(out, item)

if isempty(out)
    out = item;
else
    out(end + 1, 1) = item;
end

end

function T = vertcatOrEmpty(rows, emptyT)

if isempty(rows)

    T = emptyT;

else

    T = struct2table(vertcat(rows{:}));
end

end

function ensureDir(p)

if ~isfolder(p)
    mkdir(p);
end

end

function x = sumNoNan(v)

v = v(isfinite(v));

if isempty(v)
    x = 0;
else
    x = sum(v);
end

end

function x = meanNoNan(v)

v = v(isfinite(v));

if isempty(v)
    x = NaN;
else
    x = mean(v);
end

end


function [Tnew, nPreserved] = preserveBadFitManualAnnotations( ...
        Tnew, ...
        filename, ...
        sheetName)

nPreserved = 0;

manualFields = [ ...
    "exceptionClass", ...
    "exceptionClassLabel", ...
    "classNote"];

for i = 1:numel(manualFields)
    fieldName = manualFields(i);

    if ~ismember(fieldName, string(Tnew.Properties.VariableNames))
        Tnew.(char(fieldName)) = strings(height(Tnew), 1);
    else
        Tnew.(char(fieldName)) = annotationColumnToString( ...
            Tnew.(char(fieldName)));
    end
end

if ~isfile(filename)
    return;
end

try

    Told = readtable( ...
        filename, ...
        "Sheet", char(sheetName), ...
        "TextType", "string", ...
        "VariableNamingRule", "preserve");

catch readError

    error( ...
        "Existing bad-fit workbook could not be read safely:\n%s\n" + ...
        "The workbook was NOT overwritten, so manual classifications " + ...
        "are protected.\nOriginal read error: %s", ...
        filename, ...
        readError.message);
end

oldVars = string(Told.Properties.VariableNames);

if ~ismember("record_ID", oldVars)

    error( ...
        "Existing bad-fit workbook has no record_ID column:\n%s\n" + ...
        "The workbook was NOT overwritten because manual annotations " + ...
        "cannot be matched safely.", ...
        filename);
end

Told.record_ID = strtrim(string(Told.record_ID));
Tnew.record_ID = strtrim(string(Tnew.record_ID));

oldIDs = Told.record_ID;
validOldID = ...
    ~ismissing(oldIDs) & ...
    strlength(oldIDs) > 0;

oldIDsValid = oldIDs(validOldID);

if numel(unique(oldIDsValid)) ~= numel(oldIDsValid)

    duplicateIDs = unique( ...
        oldIDsValid( ...
            countStringOccurrences(oldIDsValid) > 1));

    error( ...
        "Existing bad-fit workbook contains duplicate record_ID values: %s\n" + ...
        "The workbook was NOT overwritten. Resolve duplicates first.", ...
        strjoin(duplicateIDs, ", "));
end

availableManualFields = manualFields( ...
    ismember(manualFields, oldVars));

if isempty(availableManualFields)
    return;
end

for i = 1:numel(availableManualFields)

    fieldName = availableManualFields(i);

    Told.(char(fieldName)) = annotationColumnToString( ...
        Told.(char(fieldName)));
end

oldIndex = containers.Map( ...
    'KeyType', 'char', ...
    'ValueType', 'double');

for i = 1:height(Told)

    id = Told.record_ID(i);

    if ismissing(id) || strlength(id) == 0
        continue;
    end

    oldIndex(char(id)) = i;
end

for i = 1:height(Tnew)

    id = Tnew.record_ID(i);

    if ismissing(id) || ...
            strlength(id) == 0 || ...
            ~isKey(oldIndex, char(id))

        continue;
    end

    oldRow = oldIndex(char(id));
    preservedThisRecord = false;

    for j = 1:numel(availableManualFields)

        fieldName = availableManualFields(j);

        oldValue = Told.(char(fieldName))(oldRow);

        if ~ismissing(oldValue) && strlength(strtrim(oldValue)) > 0

            Tnew.(char(fieldName))(i) = oldValue;
            preservedThisRecord = true;
        end
    end

    if preservedThisRecord
        nPreserved = nPreserved + 1;
    end
end

end

function s = annotationColumnToString(x)

if isstring(x)

    s = x;

elseif ischar(x)

    s = string(cellstr(x));

elseif iscell(x)

    s = strings(numel(x), 1);

    for i = 1:numel(x)
        s(i) = annotationScalarToString(x{i});
    end

elseif isnumeric(x) || islogical(x)

    s = strings(numel(x), 1);

    for i = 1:numel(x)

        if isfinite(double(x(i)))
            s(i) = string(x(i));
        end
    end

else

    s = strings(numel(x), 1);

    for i = 1:numel(x)
        s(i) = annotationScalarToString(x(i));
    end
end

s = s(:);
s(ismissing(s)) = "";

badLiteral = ...
    lower(strtrim(s)) == "nan" | ...
    lower(strtrim(s)) == "<missing>";

s(badLiteral) = "";

end

function s = annotationScalarToString(x)

s = "";

if isempty(x)
    return;
end

if isnumeric(x) || islogical(x)

    if isscalar(x) && isfinite(double(x))
        s = string(x);
    end

    return;
end

try

    sx = string(x);

    if isempty(sx)
        return;
    end

    sx = sx(1);

    if ismissing(sx)
        return;
    end

    sx = strtrim(sx);

    if lower(sx) == "nan" || lower(sx) == "<missing>"
        return;
    end

    s = sx;

catch
end

end

function counts = countStringOccurrences(x)

x = string(x(:));
counts = zeros(size(x));

for i = 1:numel(x)
    counts(i) = nnz(x == x(i));
end

end


