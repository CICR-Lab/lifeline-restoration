function run_00_prepare_restoration_shape_phase_ratio_data()




moduleRoot = string(fileparts(mfilename("fullpath")));

emergentSimplicityRoot = string(fileparts(fileparts(moduleRoot)));
analysisRoot = string(fileparts(emergentSimplicityRoot));

inputMat = fullfile( ...
    analysisRoot, ...
    "00_Dataset", ...
    "standardized_dataset.mat");

outputRoot = fullfile( ...
    emergentSimplicityRoot, ...
    "outputs", ...
    "intermediate", ...
    "01_restoration_shape_phase_ratio");

outMat = fullfile( ...
    outputRoot, ...
    "unique3_records.mat");

sourceWorkbook = fullfile( ...
    outputRoot, ...
    "SourceData_Fig3b_phase_ratios.xlsx");

broadContextMat = fullfile( ...
    outputRoot, ...
    "D0_Tp_broad_context.mat");

ensureDir(outputRoot);
calculate_D0_Tp();

fprintf("Module root:\n%s\n\n", moduleRoot);
fprintf("Prepared MAT output:\n%s\n\n", outMat);
fprintf("Phase-ratio source data:\n%s\n\n", sourceWorkbook);


if ~isfile(inputMat)
    error( ...
        "Standardized dataset MAT not found:\n%s", ...
        inputMat);
end

if ~isfile(broadContextMat)
    error( ...
        "Canonical D0/Tp MAT not found:\n%s", ...
        broadContextMat);
end

S = load(inputMat, "event_region_system_records");
B = load(broadContextMat, "D0andTp");

if ~isfield(S, "event_region_system_records")
    error( ...
        "Input MAT does not contain event_region_system_records:\n%s", ...
        inputMat);
end

if ~isfield(B, "D0andTp")
    error( ...
        "Canonical D0/Tp MAT does not contain D0andTp:\n%s", ...
        broadContextMat);
end

records = S.event_region_system_records;
D0andTp = B.D0andTp;
D0andTp.record_ID = string(D0andTp.record_ID);
P = [0.10:0.10:0.90, 0.95, 1.00];

if isempty(records)
    error("event_region_system_records is empty.");
end

if numel(unique(D0andTp.record_ID)) ~= height(D0andTp)
    error("D0_Tp_broad_context.mat contains duplicate record_ID values.");
end

requiredBroadFields = ["record_ID", "D0", "D0_source"];
missingBroadFields = requiredBroadFields( ...
    ~ismember(requiredBroadFields, string(D0andTp.Properties.VariableNames)));
if ~isempty(missingBroadFields)
    error( ...
        "D0andTp is missing required field(s): %s", ...
        strjoin(missingBroadFields, ", "));
end

requiredFields = [ ...
    "record_ID", ...
    "restorationratioCurve", ...
    "PeakTime", ...
    "System", ...
    "nRestorationRatioPoints", ...
    "nUniqueRestorationValues"];

validateRequiredRecordFields(records, requiredFields);

unique3Cells = cell(numel(records), 1);
unique3Count = 0;


for k = 1:numel(records)

    rec = records(k);
    row = baseTauRow(rec);


    peakTime = getNumericField(rec, "PeakTime");
    row.PeakTime = peakTime;
    row.nRestorationRatioPoints = getNumericField(rec, "nRestorationRatioPoints");
    row.nUniqueRestorationValues = getNumericField(rec, "nUniqueRestorationValues");

    broadIdx = find(D0andTp.record_ID == row.record_ID, 1, "first");
    if isempty(broadIdx)
        error("record_ID %s was not found in D0_Tp_broad_context.mat.", row.record_ID);
    end

    row.D0 = double(D0andTp.D0(broadIdx));
    row.D0_source = string(D0andTp.D0_source(broadIdx));

    [tRR, rRR] = curveToArrays( ...
        keepExactCurvePoints( ...
            getFieldOrDefault( ...
                rec, ...
                "restorationratioCurve", ...
                repmat(makePoint(), 0, 1))));

    [Tvals, tauVals] = computeMilestonesFromObservedRR( ...
        tRR, ...
        rRR, ...
        peakTime, ...
        P);

    for i = 1:numel(P)
        label = pLabel(P(i));
        row.("T" + label) = Tvals(i);
        row.("tau" + label) = tauVals(i);
    end


    row = addPhaseRatios(row);

    if isfinite(row.nUniqueRestorationValues) && ...
            row.nUniqueRestorationValues >= 3

        curveRec = makeCurveRecord(rec);

        curveRec = addScalarAnalysisFields(curveRec, row);

        unique3Count = unique3Count + 1;
        unique3Cells{unique3Count, 1} = curveRec;
    end
end


unique3Records = cellStructArray(unique3Cells(1:unique3Count));

if isempty(unique3Records)

    wellResolvedRecords = struct([]);

else

    nPoints = [unique3Records.nRestorationRatioPoints];
    nUnique = [unique3Records.nUniqueRestorationValues];

    wellResolvedMask = ...
        isfinite(nPoints) & ...
        isfinite(nUnique) & ...
        nPoints >= 6 & ...
        nUnique >= 3;

    wellResolvedRecords = unique3Records(wellResolvedMask);
end


phaseRatioSummary = buildPhaseRatioSummary(wellResolvedRecords);
phaseRatioRecordLevel = buildPhaseRatioRecordLevel(wellResolvedRecords);
fig3aTrajectories = buildNormalizedTrajectorySourceData(wellResolvedRecords);

save( ...
    outMat, ...
    "unique3Records", ...
    "wellResolvedRecords", ...
    "phaseRatioSummary", ...
    "phaseRatioRecordLevel", ...
    "fig3aTrajectories");

writePhaseRatioSourceData( ...
    sourceWorkbook, ...
    fig3aTrajectories, ...
    phaseRatioRecordLevel, ...
    phaseRatioSummary);


fprintf("\nSaved emergent-simplicity dataset:\n%s\n", outMat);
fprintf("Saved phase-ratio source data:\n%s\n", sourceWorkbook);

fprintf( ...
    "Unique>=3 records: %d\n", ...
    numel(unique3Records));

fprintf( ...
    "Well-resolved records (>=6 points and >=3 unique values): %d\n", ...
    numel(wellResolvedRecords));

fprintf( ...
    "Phase-ratio summary rows: %d\n", ...
    height(phaseRatioSummary));

end



function row = baseTauRow(rec)

row = struct();

row.record_ID = getStringField(rec, "record_ID");
row.event_ISO_ID = getNumericField(rec, "event_ISO_ID");
row.event_ID = getNumericField(rec, "event_ID");
row.Country = getStringField(rec, "Country");
row.ISO = getStringField(rec, "ISO");
row.IncomeGroup = getStringField(rec, "IncomeGroup");
row.OccTime = getStringField(rec, "OccTime");
row.magnitude = getNumericField(rec, "magnitude");
row.regionIndex = getNumericField(rec, "regionIndex");
row.regionName = getStringField(rec, "regionName");
row.System = normalizeSystemName(getStringField(rec, "System"));
row.avg_PGA = getNumericField(rec, "avg_PGA");
row.PeakTime = getNumericField(rec, "PeakTime");

row.D0 = NaN;
row.D0_source = "none";

row.nRestorationRatioPoints = NaN;
row.nUniqueRestorationValues = NaN;

P = [0.10:0.10:0.90, 0.95, 1.00];

for i = 1:numel(P)
    label = pLabel(P(i));
    row.("T" + label) = NaN;
    row.("tau" + label) = NaN;
end

row.tau80_over_tau50 = NaN;
row.tau90_over_tau80 = NaN;
row.tau95_over_tau90 = NaN;

end


function row = addPhaseRatios(row)

row.tau80_over_tau50 = safePositiveRatio( ...
    row.tau80, ...
    row.tau50);

row.tau90_over_tau80 = safePositiveRatio( ...
    row.tau90, ...
    row.tau80);

row.tau95_over_tau90 = safePositiveRatio( ...
    row.tau95, ...
    row.tau90);

end

function v = safePositiveRatio(numVal, denVal)

v = NaN;

if isfinite(numVal) && ...
        isfinite(denVal) && ...
        denVal > 0 && ...
        numVal >= denVal

    v = numVal / denVal;
end

end


function T = buildPhaseRatioSummary(wellResolvedRecords)

metricNames = { ...
    'tau80_over_tau50', ...
    'tau90_over_tau80', ...
    'tau95_over_tau90'};

groupNames = { ...
    'All systems', ...
    'Power', ...
    'Water', ...
    'Gas'};

nRows = numel(metricNames) * numel(groupNames);

template = struct( ...
    'metric', "", ...
    'system', "", ...
    'n', 0, ...
    'median', NaN, ...
    'Q25', NaN, ...
    'Q75', NaN, ...
    'Q05', NaN, ...
    'Q95', NaN);

rows = repmat(template, nRows, 1);
rowIndex = 0;

if isempty(wellResolvedRecords)

    systems = strings(0, 1);

else

    systems = strings(numel(wellResolvedRecords), 1);

    for i = 1:numel(wellResolvedRecords)

        systems(i) = normalizeSystemName( ...
            getStringField(wellResolvedRecords(i), 'System'));
    end
end

for m = 1:numel(metricNames)

    metricName = metricNames{m};

    allValues = nan(numel(wellResolvedRecords), 1);

    for i = 1:numel(wellResolvedRecords)

        allValues(i) = getNumericField( ...
            wellResolvedRecords(i), ...
            metricName);
    end

    for g = 1:numel(groupNames)

        groupName = groupNames{g};
        rowIndex = rowIndex + 1;

        if strcmp(groupName, 'All systems')

            values = allValues;

        else

            values = allValues(systems == string(groupName));
        end

        values = double(values(:));
        values = values(isfinite(values));

        rows(rowIndex).metric = string(metricName);
        rows(rowIndex).system = displaySystemName(groupName);
        rows(rowIndex).n = numel(values);

        if ~isempty(values)

            q = quantile( ...
                values, ...
                [0.05, 0.25, 0.50, 0.75, 0.95]);

            rows(rowIndex).Q05 = q(1);
            rows(rowIndex).Q25 = q(2);
            rows(rowIndex).median = q(3);
            rows(rowIndex).Q75 = q(4);
            rows(rowIndex).Q95 = q(5);
        end
    end
end

T = struct2table(rows);

T = T(:, { ...
    'metric', ...
    'system', ...
    'n', ...
    'median', ...
    'Q25', ...
    'Q75', ...
    'Q05', ...
    'Q95'});

end

function T = buildPhaseRatioRecordLevel(wellResolvedRecords)

metrics = ["tau80_over_tau50", "tau90_over_tau80", "tau95_over_tau90"];
rows = cell(0, 1);

for i = 1:numel(wellResolvedRecords)

    rec = wellResolvedRecords(i);

    for j = 1:numel(metrics)

        value = getNumericField(rec, metrics(j));

        if isfinite(value)
            rows{end + 1, 1} = struct( ...
                "record_ID", getStringField(rec, "record_ID"), ...
                "system", normalizeSystemName(getStringField(rec, "System")), ...
                "metric", metrics(j), ...
                "value", value);
        end
    end
end

if isempty(rows)
    T = table( ...
        strings(0, 1), ...
        strings(0, 1), ...
        strings(0, 1), ...
        nan(0, 1), ...
        "VariableNames", {'record_ID', 'system', 'metric', 'value'});
else
    T = struct2table(vertcat(rows{:}));
end

end

function T = buildNormalizedTrajectorySourceData(wellResolvedRecords)

rows = cell(0, 1);

for i = 1:numel(wellResolvedRecords)
    rec = wellResolvedRecords(i);
    tau50 = getNumericField(rec, "tau50");
    peakTime = getNumericField(rec, "PeakTime");

    if ~(isfinite(tau50) && tau50 > 0 && isfinite(peakTime))
        continue;
    end

    curve = getFieldOrDefault( ...
        rec, ...
        "restorationratioCurve", ...
        repmat(makePoint(), 0, 1));

    for j = 1:numel(curve)
        dayVal = getNumericField(curve(j), "day");
        valueVal = getNumericField(curve(j), "value");
        if ~(isfinite(dayVal) && isfinite(valueVal))
            continue;
        end

        tau = dayVal - peakTime;
        if valueVal < -1e-8 || valueVal > 1 + 1e-8
            error( ...
                "restorationratioCurve contains value outside [0,1] for record_ID %s.", ...
                getStringField(rec, "record_ID"));
        end

        if tau < 0
            continue;
        end

        rows{end + 1, 1} = struct( ...
            "record_ID", getStringField(rec, "record_ID"), ...
            "system", normalizeSystemName(getStringField(rec, "System")), ...
            "day", dayVal, ...
            "PeakTime", peakTime, ...
            "tau", tau, ...
            "tau50", tau50, ...
            "tau_over_tau50", tau / tau50, ...
            "restoration_ratio", valueVal);
    end
end

if isempty(rows)
    T = table( ...
        strings(0, 1), ...
        strings(0, 1), ...
        nan(0, 1), ...
        nan(0, 1), ...
        nan(0, 1), ...
        nan(0, 1), ...
        nan(0, 1), ...
        nan(0, 1), ...
        "VariableNames", { ...
            'record_ID', ...
            'system', ...
            'day', ...
            'PeakTime', ...
            'tau', ...
            'tau50', ...
            'tau_over_tau50', ...
            'restoration_ratio'});
else
    T = struct2table(vertcat(rows{:}));
end

end

function writePhaseRatioSourceData(filename, fig3aTrajectories, recordLevel, summary)

if isfile(filename)
    delete(filename);
end

fig3aTrajectories = standardizeSystemNamesForOutput(fig3aTrajectories);
recordLevel = standardizeSystemNamesForOutput(recordLevel);
summary = standardizeSystemNamesForOutput(summary);

writetable(fig3aTrajectories, filename, "Sheet", "Fig3a_normalized_trajectories");
writetable(recordLevel, filename, "Sheet", "record_level");
writetable(summary, filename, "Sheet", "summary");

end

function T = standardizeSystemNamesForOutput(T)

vars = string(T.Properties.VariableNames);
for j = 1:numel(vars)
    name = lower(vars(j));
    if contains(name, "system") || contains(name, "lifeline") || name == "group"
        T.(char(vars(j))) = displaySystemName(T.(char(vars(j))));
    end
end

end

function x = displaySystemName(x)

if iscell(x)
    x = string(x);
end

if isstring(x) || iscategorical(x)
    s = string(x);
    s(s == "Power") = "Electric power";
    s(s == "Water") = "Water supply";
    s(s == "Gas") = "Natural gas";
    x = s;
end

end


function recOut = makeCurveRecord(rec)

recOut = rec;

recOut.functionalityCurve = keepExactCurvePoints( ...
    getFieldOrDefault( ...
        rec, ...
        "functionalityCurve", ...
        repmat(makePoint(), 0, 1)));

recOut.restorationratioCurve = keepExactCurvePoints( ...
    getFieldOrDefault( ...
        rec, ...
        "restorationratioCurve", ...
        repmat(makePoint(), 0, 1)));

end

function recOut = addScalarAnalysisFields(recOut, row)

names = fieldnames(row);

for i = 1:numel(names)
    recOut.(names{i}) = row.(names{i});
end

end


function curve = sortCurveByDay(curve)

if isempty(curve)
    return;
end

days = nan(numel(curve), 1);

for i = 1:numel(curve)
    days(i) = getNumericField(curve(i), "day");
end

[~, order] = sort(days);
curve = curve(order);

end

function curveOut = keepExactCurvePoints(curveIn)

curveOut = repmat(makePoint(), 0, 1);

if isempty(curveIn) || ~isstruct(curveIn)
    return;
end

keepMask = false(numel(curveIn), 1);

for i = 1:numel(curveIn)

    dayVal = getNumericField(curveIn(i), "day");
    valueVal = getNumericField(curveIn(i), "value");
    censorCode = getNumericField( ...
        curveIn(i), ...
        "day_censor_code");

    keepMask(i) = ...
        isfinite(dayVal) && ...
        isfinite(valueVal) && ...
        isfinite(censorCode) && ...
        censorCode == 0;
end

curveOut = sortCurveByDay(curveIn(keepMask));

end

function [t, r] = curveToArrays(curve)

t = [];
r = [];

if isempty(curve)
    return;
end

t = double([curve.day]');
r = double([curve.value]');
mask = isfinite(t) & isfinite(r);
t = t(mask);
r = r(mask);

if any(diff(t) == 0)
    error("Restoration-ratio curve contains duplicate exact-time day values.");
end

end

function [Tvals, tauVals] = computeMilestonesFromObservedRR(t, r, peakTime, P)

Tvals = nan(size(P));
tauVals = nan(size(P));

if isempty(t) || isempty(r) || ~isfinite(peakTime)
    return;
end

for i = 1:numel(P)
    Tvals(i) = firstCrossingTime(t, r, P(i));
    if isfinite(Tvals(i))
        tauVals(i) = Tvals(i) - peakTime;
    end
end

end

function tq = firstCrossingTime(t, r, p)

tol = 1e-12;
tq = NaN;

if isempty(t) || r(1) > p + tol
    return;
end

idx = find(r >= p - tol, 1, "first");
if isempty(idx)
    return;
end

if idx == 1
    if abs(r(idx) - p) <= tol
        tq = t(idx);
    end
elseif abs(r(idx) - p) <= tol || r(idx) <= r(idx - 1)
    tq = t(idx);
else
    tq = t(idx - 1) + ...
        (p - r(idx - 1)) * ...
        (t(idx) - t(idx - 1)) / ...
        (r(idx) - r(idx - 1));
end

end

function label = pLabel(p)

if abs(p - 0.95) < 1e-12

    label = "95";

elseif abs(p - 1.00) < 1e-12

    label = "100";

else

    label = string(round(p * 100));
end

end

function s = normalizeSystemName(s)

s = string(s);
sl = lower(strtrim(s));

if sl == "e" || ...
        sl == "power"

    s = "Power";

elseif sl == "w" || ...
        sl == "water"

    s = "Water";

elseif sl == "g" || ...
        sl == "gas"

    s = "Gas";
end

end


function S = cellStructArray(C)

if isempty(C)
    S = struct([]);
else
    S = vertcat(C{:});
end

end

function validateRequiredRecordFields(records, requiredFields)

availableFields = string(fieldnames(records));

missingFields = requiredFields( ...
    ~ismember(requiredFields, availableFields));

if ~isempty(missingFields)

    error( ...
        "The standardized records are missing required fields: %s", ...
        strjoin(missingFields, ", "));
end

end

function x = getNumericField(S, fieldName)

x = NaN;

if ~isstruct(S) || ~isfield(S, fieldName)
    return;
end

v = S.(fieldName);

if isnumeric(v) || islogical(v)

    if isscalar(v) && ...
            ~isempty(v) && ...
            isfinite(double(v))

        x = double(v);
    end

elseif ischar(v) || isstring(v)

    y = str2double(string(v));

    if isscalar(y) && isfinite(y)
        x = y;
    end
end

end

function s = getStringField(S, fieldName)

s = "";

if ~isstruct(S) || ~isfield(S, fieldName)
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

function v = getFieldOrDefault(S, fieldName, defaultValue)

if isstruct(S) && isfield(S, fieldName)
    v = S.(fieldName);
else
    v = defaultValue;
end

end

function p = makePoint()

p = struct( ...
    "dayRaw", [], ...
    "day", NaN, ...
    "value", NaN, ...
    "day_censor_code", 0);

end

function ensureDir(p)

if ~isfolder(p)
    mkdir(p);
end

end
