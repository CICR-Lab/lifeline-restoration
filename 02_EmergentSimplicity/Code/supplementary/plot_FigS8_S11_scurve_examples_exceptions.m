
function plot_FigS8_S11_scurve_examples_exceptions(supplementaryRoot, figureDataRoot)

if nargin < 1 || isempty(supplementaryRoot)
    scriptDir = string(fileparts(mfilename("fullpath")));
    moduleRoot = scriptDir;
    emergentRoot = string(fileparts(fileparts(moduleRoot)));
    supplementaryRoot = fullfile(emergentRoot, "outputs", "supplementary");
end
figureOnly = nargin >= 2 && strlength(string(figureDataRoot)) > 0;

scriptDir = string(fileparts(mfilename("fullpath")));
moduleRoot = scriptDir;
emergentRoot = string(fileparts(fileparts(moduleRoot)));
if figureOnly
    intermediateRoot = string(figureDataRoot);
else
    intermediateRoot = fullfile( ...
        emergentRoot, ...
        "outputs", ...
        "intermediate", ...
        "01_restoration_shape_phase_ratio");
end

addpath(scriptDir);

curveMat = fullfile( ...
    intermediateRoot, ...
    "unique3_records.mat");

fitMat = fullfile( ...
    intermediateRoot, ...
    "restoration_shape_ge5_fit_result.mat");

badfitPath = fullfile( ...
    intermediateRoot, ...
    "badfit_weibull_gompertz_loglogistic_ge6_candidates.xlsx");

suppOutputRoot = supplementaryRoot;
if figureOnly
    figureDir = suppOutputRoot;
    sourceDataDir = string(figureDataRoot);
else
    figureDir = fullfile(suppOutputRoot, "figures");
    sourceDataDir = fullfile(suppOutputRoot, "source_data");
end
manualInputPath = fullfile( ...
    emergentRoot, ...
    "Code", ...
    "exception_classes.xlsx");
sourceWorkbook = fullfile( ...
    sourceDataDir, ...
    "SourceData_FigS8_S11.xlsx");
if figureOnly
    % Use classifications captured with the selected figure-data source.
    manualInputPath = sourceWorkbook;
end

ensureDir(figureDir);
ensureDir(sourceDataDir);

if ~isfile(curveMat)
    error( ...
        "Curve MAT not found:\n%s\n" + ...
        "Run run_00_prepare_well_resolved_curves.m first.", ...
        curveMat);
end

if ~isfile(fitMat)
    error( ...
        "Curve-fit MAT not found:\n%s\n" + ...
        "Run run_01_fit_restoration_shape_models.m first.", ...
        fitMat);
end

if ~isfile(badfitPath)
    error( ...
        "Bad-fit candidate workbook not found:\n%s\n" + ...
        "Run run_01_fit_restoration_shape_models.m first.", ...
        badfitPath);
end

if ~isfile(manualInputPath)
    error( ...
        "Manual exception-class input not found:\n%s\n" + ...
        "Create it from the confirmed FigS9_S11_exceptions classifications " + ...
        "and follow the README instructions.", ...
        manualInputPath);
end


figS4RecordIDs = [ ...
    "E1_R2_Power"
    "E15_R6_Power"
    "E193_R12_Power"
    "E2_R1_Water"
    "E2_R8_Water"
    "E3_R2_Water"
    "E10_R1_Gas"
    "E24_R2_Gas"
    "E6_R7_Gas"];
if figureOnly
    selectedRecords = readtable(sourceWorkbook, "Sheet", "FigS8_records", ...
        "TextType", "string", "VariableNamingRule", "preserve");
    assertRequiredColumns(selectedRecords, "record_ID", "Figure S8 source data");
    figS4RecordIDs = string(selectedRecords.record_ID);
end


S = load(curveMat, "unique3Records");

if ~isfield(S, "unique3Records")
    error( ...
        "MAT file does not contain unique3Records:\n%s", ...
        curveMat);
end

curves = S.unique3Records;

if isempty(curves)
    error("unique3Records is empty.");
end

Fit = load(fitMat, "fitResult");

if ~isfield(Fit, "fitResult") || ~isfield(Fit.fitResult, "CurveFitParams")
    error("MAT file does not contain fitResult.CurveFitParams:\n%s", fitMat);
end

F = Fit.fitResult.CurveFitParams;

F = normalizeFitTable(F);

requiredFitVars = [ ...
    "record_ID", ...
    "System", ...
    "model", ...
    "MAE", ...
    "RMSE", ...
    "R2", ...
    "lambda", ...
    "k", ...
    "b", ...
    "c"];

assertRequiredColumns( ...
    F, ...
    requiredFitVars, ...
    "restoration_shape_ge5_curve_fit_params.csv");

Bad = readtable( ...
    badfitPath, ...
    "Sheet", "badfit_candidates", ...
    "TextType", "string", ...
    "VariableNamingRule", "preserve");

Bad = normalizeBadfitTable(Bad);
Bad = applyExceptionClassInput(Bad, manualInputPath);

requiredBadVars = [ ...
    "record_ID", ...
    "System", ...
    "exceptionClass", ...
    "MAE_weibull", ...
    "MAE_gompertz", ...
    "MAE_loglogistic"];

assertRequiredColumns( ...
    Bad, ...
    requiredBadVars, ...
    "bad-fit candidate workbook");

validateExceptionClasses(Bad);


models = [ ...
    "Weibull", ...
    "Gompertz", ...
    "Log-logistic"];

modelLabels = models;
lineStyles = ["-", "--", ":"];

systemOrder = ["Power", "Water", "Gas"];
systemNames = ["Electric power", "Water supply", "Natural gas"];

systemColors = [ ...
    0.70 0.14 0.16
    0.17 0.46 0.72
    0.78 0.55 0.18];

dataColor = [0.70 0.14 0.16];
fitColor = [0.17 0.46 0.72];

maeThreshold = 0.05;


Sel = buildSelectedRecordTable( ...
    figS4RecordIDs, ...
    curves, ...
    F, ...
    models);

Sel.isAdequateAll3 = ...
    Sel.MAE_weibull <= maeThreshold & ...
    Sel.MAE_gompertz <= maeThreshold & ...
    Sel.MAE_loglogistic <= maeThreshold;

if any(~Sel.isAdequateAll3)

    warning( ...
        "Some fixed Fig. S8 records no longer have MAE <= %.3f for all " + ...
        "three core models: %s", ...
        maeThreshold, ...
        strjoin( ...
            Sel.record_ID(~Sel.isAdequateAll3), ...
            ", "));
end


Bad = sortExceptionBySystemAndWorkbookOrder( ...
    Bad, ...
    systemOrder);


plotCurveGridFixed( ...
    Sel, ...
    curves, ...
    F, ...
    models, ...
    modelLabels, ...
    lineStyles, ...
    systemOrder, ...
    systemNames, ...
    systemColors, ...
    dataColor, ...
    fitColor, ...
    fullfile(figureDir, "FigS8_scurve_examples"), ...
    3, ...
    3, ...
    "FigS8");

classNames = [ ...
    "Long tail", ...
    "Near threshold", ...
    "Unknown"];

figureNames = [ ...
    "FigS9_exception_long_tail", ...
    "FigS10_exception_near_threshold", ...
    "FigS11_exception_unknown"];

fixedRows = [2, 2, 3];
fixedCols = [2, 3, 5];

for c = 1:numel(classNames)

    D = Bad(Bad.exceptionClass == classNames(c), :);

    capacity = fixedRows(c) * fixedCols(c);

    if height(D) > capacity

        error( ...
            "%s contains %d records but the fixed grid has only %d tiles.", ...
            classNames(c), ...
            height(D), ...
            capacity);
    end

    plotCurveGridFixed( ...
        D, ...
        curves, ...
        F, ...
        models, ...
        modelLabels, ...
        lineStyles, ...
        systemOrder, ...
        systemNames, ...
        systemColors, ...
        dataColor, ...
        fitColor, ...
        fullfile(figureDir, figureNames(c)), ...
        fixedRows(c), ...
        fixedCols(c), ...
        "FigS" + string(c + 8));
end

writeSourceWorkbook = ~figureOnly;
if writeSourceWorkbook && isfile(sourceWorkbook)
    try
        delete(sourceWorkbook);
    catch
        writeSourceWorkbook = false;
    end
    if isfile(sourceWorkbook)
        writeSourceWorkbook = false;
    end
end

if writeSourceWorkbook
    writeTableCompatPortable( ...
        Sel, ...
        sourceWorkbook, ...
        "FigS8_records");

    writeTableCompatPortable( ...
        exceptionSourceTable(Bad), ...
        sourceWorkbook, ...
        "FigS9_S11_exceptions");

elseif ~figureOnly
    warning( ...
        "Source workbook is open or locked and was left unchanged:\n%s", ...
        sourceWorkbook);
end

fprintf( ...
    "Saved Fig S8-S11 outputs under:\n%s\n", ...
    suppOutputRoot);


end

function T = normalizeFitTable(T)

T = normalizeStrings(T);

T = renameVarIfPresent(T, "curve_id", "record_ID");
T = renameVarIfPresent(T, "system", "System");

T.record_ID = strtrim(string(T.record_ID));
T.System = canonicalSystemName(T.System);
T.model = canonicalModelName(T.model);

end

function T = normalizeBadfitTable(T)

T = normalizeStrings(T);

T = renameVarIfPresent(T, "curveID", "record_ID");
T = renameVarIfPresent(T, "curve_id", "record_ID");
T = renameVarIfPresent(T, "system", "System");
T = renameVarIfPresent(T, "exceptionclass", "exceptionClass");
T = renameVarIfPresent(T, "exception_class", "exceptionClass");
T = renameVarIfPresent(T, "exceptionclasslabel", "exceptionClassLabel");
T = renameVarIfPresent(T, "exception_class_label", "exceptionClassLabel");

T.record_ID = strtrim(string(T.record_ID));
T.System = canonicalSystemName(T.System);

if ~ismember( ...
        "exceptionClass", ...
        string(T.Properties.VariableNames))

    T.exceptionClass = strings(height(T), 1);
end

if ismember( ...
        "exceptionClassLabel", ...
        string(T.Properties.VariableNames))

    labelClass = standardClassName(T.exceptionClassLabel);
    useLabel = ...
        ~ismissing(labelClass) & ...
        strlength(strtrim(labelClass)) > 0;

    T.exceptionClass(useLabel) = labelClass(useLabel);
end

T.exceptionClass = standardClassName(T.exceptionClass);
T.workbookOrder = transpose((1:height(T)));

end

function T = applyExceptionClassInput(T, inputPath)

A = readExceptionClassInput(inputPath);
if isempty(A)
    error( ...
        "Manual exception-class input is empty or has no FigS9_S11_exceptions sheet:\n%s", ...
        inputPath);
end

A = normalizeStrings(A);
A = renameVarIfPresent(A, "curveID", "record_ID");
A = renameVarIfPresent(A, "curve_id", "record_ID");
A = renameVarIfPresent(A, "exceptionclass", "exceptionClass");
A = renameVarIfPresent(A, "exception_class", "exceptionClass");
A = renameVarIfPresent(A, "exceptionclasslabel", "exceptionClassLabel");
A = renameVarIfPresent(A, "exception_class_label", "exceptionClassLabel");

assertRequiredColumns(A, "record_ID", "exception-class input");

A.record_ID = strtrim(string(A.record_ID));

if numel(unique(A.record_ID)) ~= height(A)
    duplicateIDs = unique(A.record_ID( ...
        arrayfun(@(x) nnz(A.record_ID == x) > 1, A.record_ID)));
    error( ...
        "Manual exception-class input contains duplicate record_ID values: %s", ...
        strjoin(duplicateIDs, ", "));
end

extraIDs = setdiff(A.record_ID, T.record_ID);
if ~isempty(extraIDs)
    error( ...
        "Manual exception-class input contains record(s) absent from the current " + ...
        "candidate set: %s", ...
        strjoin(extraIDs, ", "));
end

if ~ismember("exceptionClass", string(A.Properties.VariableNames))
    A.exceptionClass = strings(height(A), 1);
end
if ~ismember("exceptionClassLabel", string(A.Properties.VariableNames))
    A.exceptionClassLabel = strings(height(A), 1);
end
if ~ismember("plotOrder", string(A.Properties.VariableNames))
    A.plotOrder = nan(height(A), 1);
end
A.rowOrder = transpose((1:height(A)));

classFromText = standardClassName(A.exceptionClass);
classFromLabel = standardClassName(A.exceptionClassLabel);
useLabel = ~ismissing(classFromLabel) & strlength(strtrim(classFromLabel)) > 0;
classFromText(useLabel) = classFromLabel(useLabel);

[matched, loc] = ismember(T.record_ID, A.record_ID);
if any(~matched)
    error( ...
        "Exception-class input is missing %d bad-fit record(s): %s", ...
        nnz(~matched), ...
        strjoin(T.record_ID(~matched), ", "));
end

T.exceptionClass = classFromText(loc);
plotOrder = double(A.plotOrder(loc));
rowOrder = double(A.rowOrder(loc));
plotOrder(~isfinite(plotOrder)) = rowOrder(~isfinite(plotOrder));
T.plotOrder = plotOrder;

missingMask = ismissing(T.exceptionClass) | strlength(strtrim(T.exceptionClass)) == 0;
if any(missingMask)
    error( ...
        "The FigS9_S11_exceptions sheet still contains %d unclassified record(s): %s\n" + ...
        "Complete exceptionClass as Long tail, Near threshold or Unknown.", ...
        nnz(missingMask), ...
        strjoin(T.record_ID(missingMask), ", "));
end

end

function A = readExceptionClassInput(inputPath)

A = table();
if ~isfile(inputPath)
    return;
end

try
    A = readtable( ...
        inputPath, ...
        "Sheet", "FigS9_S11_exceptions", ...
        "TextType", "string", ...
        "VariableNamingRule", "preserve");
catch
end

end

function T = exceptionSourceTable(T)

removeVars = [ ...
    "exceptionClassLabel", ...
    "classNote", ...
    "badfitCriterion", ...
    "bestMAE_Weibull_Gompertz_Loglogistic", ...
    "bestModel_Weibull_Gompertz_Loglogistic", ...
    "plotOrder", ...
    "plotOrderForPlot", ...
    "workbookOrder"];

vars = string(T.Properties.VariableNames);
T(:, intersect(removeVars, vars, "stable")) = [];

end

function validateExceptionClasses(T)

raw = strtrim(string(T.exceptionClass));

missingMask = ...
    ismissing(raw) | ...
    strlength(raw) == 0;

if any(missingMask)

    error( ...
        "The bad-fit workbook still contains %d unclassified record(s): %s\n" + ...
        "Complete exceptionClass using Long tail, Near threshold or Unknown.", ...
        nnz(missingMask), ...
        strjoin(T.record_ID(missingMask), ", "));
end

allowed = [ ...
    "Long tail", ...
    "Near threshold", ...
    "Unknown"];

invalidMask = ~ismember(raw, allowed);

if any(invalidMask)

    error( ...
        "Unrecognized exceptionClass value(s): %s", ...
        strjoin(unique(raw(invalidMask)), ", "));
end

end


function T = buildSelectedRecordTable( ...
        recordIDs, ...
        curves, ...
        F, ...
        models)

n = numel(recordIDs);

record_ID = string(recordIDs(:));
inputOrder = transpose((1:n));

System = strings(n, 1);

MAE_weibull = nan(n, 1);
MAE_gompertz = nan(n, 1);
MAE_loglogistic = nan(n, 1);

for i = 1:n

    rec = getCurve(curves, record_ID(i));

    System(i) = canonicalSystemName( ...
        getStringField(rec, "System"));

    mae = nan(1, numel(models));

    for m = 1:numel(models)

        row = fitRow( ...
            F, ...
            record_ID(i), ...
            models(m));

        if isempty(row)

            error( ...
                "No %s fit was found for fixed Fig. S4 record %s.", ...
                models(m), ...
                record_ID(i));
        end

        mae(m) = firstScalar(row.MAE);
    end

    MAE_weibull(i) = mae(1);
    MAE_gompertz(i) = mae(2);
    MAE_loglogistic(i) = mae(3);
end

T = table( ...
    record_ID, ...
    inputOrder, ...
    System, ...
    MAE_weibull, ...
    MAE_gompertz, ...
    MAE_loglogistic);

end


function D = sortExceptionBySystemAndWorkbookOrder(D, systemOrder)

if ismember("plotOrder", string(D.Properties.VariableNames)) && ...
        any(isfinite(double(D.plotOrder)))

    D.plotOrderForPlot = double(D.plotOrder);
    D.plotOrderForPlot(~isfinite(D.plotOrderForPlot)) = inf;
    D = sortrows( ...
        D, ...
        ["plotOrderForPlot", "workbookOrder"]);
    D.plotOrderForPlot = [];
    return;
end

systemRank = nan(height(D), 1);

for i = 1:numel(systemOrder)
    systemRank(D.System == systemOrder(i)) = i;
end

systemRank(~isfinite(systemRank)) = ...
    numel(systemOrder) + 1;

D.systemRankForPlot = systemRank;

D = sortrows( ...
    D, ...
    ["systemRankForPlot", "workbookOrder"]);

D.systemRankForPlot = [];

end


function plotCurveGridFixed( ...
        D, ...
        curves, ...
        F, ...
        models, ...
        modelLabels, ...
        lineStyles, ...
        systemOrder, ...
        systemNames, ...
        systemColors, ...
        dataColor, ...
        fitColor, ...
        outBase, ...
        nRows, ...
        nCols, ...
        figTag)

D = reorderRecordsForFigure( ...
    D, ...
    figTag);

[figWidth, figHeight, tileSpacingMode, paddingMode] = ...
    figureLayoutSettings( ...
        figTag, ...
        nRows, ...
        nCols);

fig = figure( ...
    "Color", "w", ...
    "Units", "pixels", ...
    "Position", [60 40 figWidth figHeight], ...
    "Visible", "off", ...
    "Resize", "off");

set( ...
    fig, ...
    "DefaultAxesFontName", "Arial", ...
    "DefaultTextFontName", "Arial");

tl = tiledlayout( ...
    fig, ...
    nRows, ...
    nCols, ...
    "TileSpacing", tileSpacingMode, ...
    "Padding", paddingMode);

legendHandles = gobjects(0);

legendLabels = [ ...
    "Data", ...
    "Weibull", ...
    "Gompertz", ...
    "Log-logistic"];

for i = 1:(nRows * nCols)

    ax = nexttile(tl);
    hold(ax, "on");

    if i > height(D)

        axis(ax, "off");
        continue;
    end

    rec = getCurve( ...
        curves, ...
        D.record_ID(i));

    recordSystem = getRecordSystem(D, i, rec);

    systemIndex = find( ...
        systemOrder == recordSystem, ...
        1);

    if isempty(systemIndex)
        systemIndex = 1;
    end

    reservedBox = legendReservedBox( ...
        figTag, ...
        i);

    legendHandles = drawOneCurvePanel( ...
        ax, ...
        rec, ...
        D(i, :), ...
        F, ...
        models, ...
        modelLabels, ...
        lineStyles, ...
        systemNames(systemIndex), ...
        systemColors(systemIndex, :), ...
        dataColor, ...
        fitColor, ...
        legendHandles, ...
        reservedBox, ...
        figTag, ...
        i, ...
        nCols);

    if figTag == "FigS8" && i == 5

        placeLegend( ...
            ax, ...
            legendHandles, ...
            legendLabels, ...
            "northwest");

    elseif figTag == "FigS11" && i == 8

        placeLegend( ...
            ax, ...
            legendHandles, ...
            legendLabels, ...
            "northwest");

    elseif figTag == "FigS9" && i == 2

        placeLegend( ...
            ax, ...
            legendHandles, ...
            legendLabels, ...
            "northeast");

    elseif figTag == "FigS10" && i == 2

        placeLegend( ...
            ax, ...
            legendHandles, ...
            legendLabels, ...
            "southeast");
    end
end

set( ...
    findall(fig, "Type", "axes"), ...
    "FontName", "Arial", ...
    "FontSize", 15, ...
    "LineWidth", 0.8, ...
    "Layer", "top");

set( ...
    findall(fig, "Type", "text"), ...
    "FontName", "Arial", ...
    "FontSize", 15);

set(findall(fig, "-property", "FontName"), "FontName", "Arial");
set(findall(fig, "-property", "FontSize"), "FontSize", 15);

set(fig, "Visible", "on");
drawnow;

saveas(fig, char(outBase + ".fig"), "fig");

exportgraphics( ...
    fig, ...
    char(outBase + ".png"), ...
    "Resolution", 300);

try

    exportgraphics( ...
        fig, ...
        char(outBase + ".pdf"), ...
        "ContentType", "vector");

catch

    print( ...
        fig, ...
        char(outBase + ".pdf"), ...
        "-dpdf", ...
        "-painters");
end

close(fig);

end

function systemName = getRecordSystem(D, rowIndex, rec)

systemName = "";

if ismember("System", string(D.Properties.VariableNames))
    systemName = canonicalSystemName(D.System(rowIndex));
end

if strlength(systemName) == 0
    systemName = canonicalSystemName(getStringField(rec, "System"));
end

end

function legendHandles = drawOneCurvePanel( ...
        ax, ...
        rec, ...
        rowInfo, ...
        F, ...
        models, ...
        modelLabels, ...
        lineStyles, ...
        systemName, ...
        systemColor, ...
        dataColor, ...
        fitColor, ...
        legendHandles, ...
        reservedBox, ...
        figTag, ...
        tileIndex, ...
        nCols)

[tauObs, rObs] = prepareCurveForPlot(rec);

if isempty(tauObs)

    axis(ax, "off");

    text( ...
        ax, ...
        0.5, ...
        0.5, ...
        "No usable restoration-ratio data", ...
        "HorizontalAlignment", "center", ...
        "VerticalAlignment", "middle", ...
        "FontName", "Arial", ...
        "FontSize", 15);

    return;
end

hData = scatter( ...
    ax, ...
    tauObs, ...
    rObs, ...
    18, ...
    dataColor, ...
    "filled", ...
    "MarkerFaceAlpha", 0.82, ...
    "MarkerEdgeColor", "w", ...
    "LineWidth", 0.3, ...
    "DisplayName", "Data");

modelHandles = gobjects(numel(models), 1);

tauGrid = linspace( ...
    min(tauObs), ...
    max(tauObs), ...
    240)';

for m = 1:numel(models)

    row = fitRow( ...
        F, ...
        rowInfo.record_ID, ...
        models(m));

    if isempty(row)
        continue;
    end

    y = modelPrediction( ...
        models(m), ...
        tauGrid, ...
        row);

    modelHandles(m) = plot( ...
        ax, ...
        tauGrid, ...
        y, ...
        lineStyles(m), ...
        "Color", fitColor, ...
        "LineWidth", 1.35, ...
        "DisplayName", modelLabels(m));
end

if isempty(legendHandles)

    keep = isgraphics(modelHandles);

    legendHandles = [ ...
        hData; ...
        modelHandles(keep)];
end

xlabel(ax, "$\tau$ (days)", "Interpreter", "latex", "FontName", "Arial", "FontSize", 15);
ylabel(ax, "$R(\tau)$", "Interpreter", "latex", "FontName", "Arial", "FontSize", 15);
title(ax, systemName, "Interpreter", "none", "FontName", "Arial", "FontSize", 15, "FontWeight", "normal");

addPanelLetter(ax, tileIndex, figTag);

ylim(ax, [0 1.03]);

if min(tauObs) == max(tauObs)

    xlim( ...
        ax, ...
        [min(tauObs) - 0.5, max(tauObs) + 0.5]);

else

    xlim( ...
        ax, ...
        [min(tauObs), max(tauObs) + eps]);
end

box(ax, "off");
grid(ax, "off");

mae = extractCoreMae( ...
    rowInfo, ...
    F, ...
    models);

maeText = sprintf( ...
    "Wei %.3f\nGom %.3f\nLL %.3f", ...
    mae(1), ...
    mae(2), ...
    mae(3));

[maeX, maeY, maeHorizontalAlignment, maeVerticalAlignment] = ...
    maeTextPlacement(figTag, tileIndex, nCols);

text( ...
    ax, ...
    maeX, ...
    maeY, ...
    maeText, ...
    "Units", "normalized", ...
    "HorizontalAlignment", maeHorizontalAlignment, ...
    "VerticalAlignment", maeVerticalAlignment, ...
    "Interpreter", "none", ...
    "FontName", "Arial", ...
    "FontSize", 15, ...
    "FontWeight", "normal", ...
    "Color", [0.10 0.10 0.10], ...
    "BackgroundColor", "w", ...
    "Margin", 0.8, ...
    "Clipping", "on");

reservedBox = reservedBox;
systemColor = systemColor;

end


function addPanelLetter(ax, tileIndex, figTag)

letters = "abcdefghijklmnopqrstuvwxyz";

if tileIndex <= strlength(letters)
    panelLabel = extractBetween(letters, tileIndex, tileIndex);
else
    panelLabel = string(tileIndex);
end

labelY = 1.105;

if string(figTag) == "FigS9" || string(figTag) == "FigS10"
    labelY = 1.045;
end

text( ...
    ax, ...
    -0.135, ...
    labelY, ...
    panelLabel, ...
    "Units", "normalized", ...
    "HorizontalAlignment", "left", ...
    "VerticalAlignment", "top", ...
    "Interpreter", "none", ...
    "FontName", "Arial", ...
    "FontSize", 17, ...
    "FontWeight", "bold", ...
    "Color", [0.05 0.05 0.05], ...
    "Clipping", "off");

end

function [x, y, hAlign, vAlign] = maeTextPlacement(figTag, tileIndex, nCols)

figTag = string(figTag);
rowIndex = ceil(double(tileIndex) / double(nCols));

x = 0.055;
y = 0.965;
hAlign = "left";
vAlign = "top";

if figTag == "FigS8" || figTag == "FigS9" || ...
        (figTag == "FigS11" && rowIndex == 2)

    x = 0.945;
    y = 0.055;
    hAlign = "right";
    vAlign = "bottom";
end

end

function mae = extractCoreMae(rowInfo, F, models)

mae = nan(1, numel(models));

fieldNames = [ ...
    "MAE_weibull", ...
    "MAE_gompertz", ...
    "MAE_loglogistic"];

for m = 1:numel(models)

    fieldName = fieldNames(m);

    if ismember( ...
            fieldName, ...
            string(rowInfo.Properties.VariableNames))

        value = rowInfo.(char(fieldName));

        if ~isempty(value)
            mae(m) = firstScalar(value);
        end
    end

    if ~isfinite(mae(m))

        row = fitRow( ...
            F, ...
            rowInfo.record_ID, ...
            models(m));

        if ~isempty(row)
            mae(m) = firstScalar(row.MAE);
        end
    end
end

end

function [tau, r] = prepareCurveForPlot(rec)

C = cleanCurve( ...
    getRawCurve(rec));

if isempty(C)

    tau = [];
    r = [];
    return;
end

tAbs = C(:, 1);
r = C(:, 2);

peakTime = getNumericField(rec, "PeakTime");
tau100 = getNumericField(rec, "tau100");

if ~isfinite(peakTime)

    tau = [];
    r = [];
    return;
end

[tAbs, r] = truncateAtPreparedTau100( ...
    tAbs, ...
    r, ...
    peakTime, ...
    tau100);

tau = max(tAbs - peakTime, 0);

end

function raw = getRawCurve(rec)

raw = [];

if isstruct(rec) && ...
        isfield(rec, "restorationratioCurve")

    raw = rec.restorationratioCurve;
end

end

function C = cleanCurve(C)

if isempty(C)
    C = zeros(0, 2);
    return;
end

if isstruct(C)

    buf = nan(numel(C), 2);

    for i = 1:numel(C)

        if isfield(C, "day")
            buf(i, 1) = getNumericField(C(i), "day");
        end

        if isfield(C, "value")
            buf(i, 2) = getNumericField(C(i), "value");
        end

        if isfield(C, "day_censor_code") && ...
                getNumericField(C(i), "day_censor_code") ~= 0

            buf(i, :) = NaN;
        end
    end

    C = buf;
end

C = double(C);

if size(C, 2) < 2

    C = zeros(0, 2);
    return;
end

C = C(:, 1:2);

mask = ...
    isfinite(C(:, 1)) & ...
    isfinite(C(:, 2));

C = C(mask, :);

if isempty(C)
    return;
end

C(:, 2) = min(max(C(:, 2), 0), 1);

[~, order] = sort(C(:, 1));
C = C(order, :);

[~, ia] = unique( ...
    C(:, 1), ...
    "last");

C = C(sort(ia), :);

end

function [t, r] = truncateAtPreparedTau100( ...
        t, ...
        r, ...
        peakTime, ...
        tau100)

if ~( ...
        isfinite(peakTime) && ...
        isfinite(tau100) && ...
        tau100 >= 0)

    return;
end

T100 = peakTime + tau100;
timeComparisonTol = 1e-10;
timeTol = timeComparisonTol * max(1, abs(T100));

keep = t <= T100 + timeTol;

t = t(keep);
r = r(keep);

end


function boxNorm = legendReservedBox(figTag, tileIndex)
boxNorm = [];
figTag = figTag;
tileIndex = tileIndex;
end


function D = reorderRecordsForFigure(D, figTag)

if figTag ~= "FigS11"
    return;
end

D = swapRowsByRecordID( ...
    D, ...
    "E2_R2_Water", ...
    "E9_R1_Water");

D = swapRowsByRecordID( ...
    D, ...
    "E13_R10_Water", ...
    "E35_R6_Water");

D = swapRowsByRecordID( ...
    D, ...
    "E6_R12_Water", ...
    "E13_R18_Water");

D = swapRowsByRecordID( ...
    D, ...
    "E3_R6_Water", ...
    "E13_R31_Water");

D = swapRowsByIndex(D, 6, 11);

end

function [figWidth, figHeight, tileSpacingMode, paddingMode] = ...
        figureLayoutSettings(figTag, nRows, nCols)

figWidth = max(1480, 270 * nCols + 180);
figHeight = 360 * nRows + 180;

tileSpacingMode = "loose";
paddingMode = "compact";

switch string(figTag)

    case "FigS8"

        figWidth = max(1500, 280 * nCols + 190);
        figHeight = 545 * nRows + 235;

    case "FigS9"

        figWidth = max(1480, 275 * nCols + 185);
        figHeight = 525 * nRows + 230;

    case "FigS10"

        figWidth = max(1480, 275 * nCols + 185);
        figHeight = 495 * nRows + 215;
        tileSpacingMode = "compact";

    case "FigS11"

        figWidth = max(1560, 295 * nCols + 190);
        figHeight = 485 * nRows + 230;
end

end




function D = swapRowsByRecordID(D, idA, idB)

ids = string(D.record_ID);

ia = find(ids == string(idA), 1, "first");
ib = find(ids == string(idB), 1, "first");

if ~isempty(ia) && ~isempty(ib)

    order = 1:height(D);
    order([ia, ib]) = order([ib, ia]);

    D = D(order, :);
end

end

function D = swapRowsByIndex(D, ia, ib)

if height(D) >= max(ia, ib)
    order = 1:height(D);
    order([ia, ib]) = order([ib, ia]);
    D = D(order, :);
end

end

function lgd = placeLegend( ...
        ax, ...
        legendHandles, ...
        legendLabels, ...
        location)

legendHandles = legendHandles( ...
    isgraphics(legendHandles));

if isempty(legendHandles)

    lgd = matlab.graphics.illustration.Legend.empty;
    return;
end

labels = legendLabels( ...
    1:min(numel(legendLabels), numel(legendHandles)));

lgd = legend( ...
    ax, ...
    legendHandles, ...
    labels, ...
    "Location", location, ...
    "Box", "off", ...
    "FontName", "Arial", ...
    "FontSize", 15, ...
    "Interpreter", "none");

lgd.ItemTokenSize = [18 10];

end


function row = fitRow(F, recordID, modelName)

recordID = strtrim(string(recordID));
modelName = canonicalModelName(modelName);

mask = ...
    strtrim(string(F.record_ID)) == recordID & ...
    canonicalModelName(F.model) == modelName;

row = F(mask, :);

if height(row) > 1

    error( ...
        "Multiple fit rows found for record %s and model %s.", ...
        recordID, ...
        modelName);
end

end

function y = modelPrediction(modelName, tau, row)

tau = max(tau, 0);

switch canonicalModelName(modelName)

    case "Weibull"

        lambda = firstScalar(row.lambda);
        k = firstScalar(row.k);

        y = 1 - exp(-((tau ./ lambda) .^ k));

    case "Gompertz"

        b = firstScalar(row.b);
        c = firstScalar(row.c);

        y = 1 - exp( ...
            -c .* ...
            (exp(b .* tau) - 1));

    case "Log-logistic"

        lambda = firstScalar(row.lambda);
        k = firstScalar(row.k);

        y = zeros(size(tau));

        mask = ...
            tau > 0 & ...
            lambda > 0 & ...
            k > 0;

        z = ...
            k .* ...
            (log(tau(mask)) - log(lambda));

        y(mask) = ...
            1 ./ ...
            (1 + exp(-z));

    otherwise

        y = nan(size(tau));
end

y(~isfinite(y)) = NaN;
y = min(max(y, 0), 1);

end

function v = firstScalar(x)

if isempty(x)
    v = NaN;
else
    v = double(x(1));
end

end


function rec = getCurve(curves, recordID)

ids = strings(numel(curves), 1);

for i = 1:numel(curves)
    ids(i) = getStringField(curves(i), "record_ID");
end

idx = ids == string(recordID);

if ~any(idx)

    error( ...
        "Curve not found in unique3Records for record_ID: %s", ...
        recordID);
end

rec = curves(find(idx, 1, "first"));

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


function name = canonicalModelName(x)

x = lower(strtrim(string(x)));
name = strings(size(x));

for i = 1:numel(x)

    switch x(i)

        case "weibull"

            name(i) = "Weibull";

        case "gompertz"

            name(i) = "Gompertz";

        case {"loglogistic", "log-logistic", "log logistic"}

            name(i) = "Log-logistic";

        otherwise

            name(i) = string(x(i));
    end
end

end

function sys = canonicalSystemName(x)

x = lower(strtrim(string(x)));
sys = strings(size(x));

for i = 1:numel(x)

    switch x(i)

        case {"e", "power"}

            sys(i) = "Power";

        case {"w", "water"}

            sys(i) = "Water";

        case {"g", "gas"}

            sys(i) = "Gas";

        otherwise

            sys(i) = string(x(i));
    end
end

end

function cls = standardClassName(x)

x = lower(strtrim(string(x)));
cls = strings(size(x));

for i = 1:numel(x)

    if ismissing(x(i)) || strlength(x(i)) == 0

        cls(i) = "";

    elseif x(i) == "1" || ...
            contains(x(i), "long tail") || ...
            x(i) == "tail"

        cls(i) = "Long tail";

    elseif x(i) == "2" || ...
            contains(x(i), "near threshold") || ...
            x(i) == "threshold"

        cls(i) = "Near threshold";

    elseif x(i) == "3" || ...
            contains(x(i), "unknown")

        cls(i) = "Unknown";

    else

        cls(i) = string(x(i));
    end
end

end

function T = normalizeStrings(T)

for i = 1:numel(T.Properties.VariableNames)

    name = T.Properties.VariableNames{i};
    v = T.(name);

    if iscellstr(v) || ischar(v) || isstring(v)
        T.(name) = string(v);
    end
end

end

function T = renameVarIfPresent(T, oldName, newName)

vars = string(T.Properties.VariableNames);

if any(vars == oldName) && ...
        ~any(vars == newName)

    idx = find( ...
        vars == oldName, ...
        1, ...
        "first");

    T.Properties.VariableNames{idx} = ...
        char(newName);
end

end

function assertRequiredColumns(T, required, label)

available = string(T.Properties.VariableNames);

missing = required( ...
    ~ismember(required, available));

if ~isempty(missing)

    error( ...
        "%s is missing required column(s): %s", ...
        label, ...
        strjoin(missing, ", "));
end

end

function ensureDir(p)

if exist(p, "dir") ~= 7
    mkdir(p);
end

end
