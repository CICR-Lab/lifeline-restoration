function run_01_make_cross_lifeline_timing_coupling_source_data()


moduleRoot = string(fileparts(mfilename("fullpath")));
emergentSimplicityRoot = string(fileparts(fileparts(moduleRoot)));
analysisRoot = string(fileparts(emergentSimplicityRoot));
intermediateRoot = fullfile( ...
    emergentSimplicityRoot, ...
    "outputs", ...
    "intermediate");

inputMat = fullfile( ...
    analysisRoot, ...
    "00_Dataset", ...
    "standardized_dataset.mat");

broadContextMat = fullfile( ...
    intermediateRoot, ...
    "01_restoration_shape_phase_ratio", ...
    "D0_Tp_broad_context.mat");

outputRoot = fullfile( ...
    intermediateRoot, ...
    "03_cross_lifeline_timing_coupling");

outMat = fullfile( ...
    outputRoot, ...
    "paired_timing_coupling_dataset.mat");
sourceFig3e = fullfile( ...
    outputRoot, ...
    "SourceData_Fig3e_cross_lifeline_timing_ratios.xlsx");

sourceFig3f = fullfile( ...
    outputRoot, ...
    "SourceData_Fig3f_cross_lifeline_coupling_Cmax90.xlsx");

ensureDir(outputRoot);


if ~isfile(inputMat)
    error( ...
        "Standardized dataset MAT not found:\n%s", ...
        inputMat);
end

S = load( ...
    inputMat, ...
    "event_region_system_records");

if ~isfield(S, "event_region_system_records")
    error( ...
        "Input MAT does not contain event_region_system_records:\n%s", ...
        inputMat);
end

records = S.event_region_system_records;

requiredFields = [ ...
    "record_ID", ...
    "event_ISO_ID", ...
    "event_ID", ...
    "regionIndex", ...
    "System"];

missingFields = requiredFields( ...
    ~isfield(records, cellstr(requiredFields)));

if ~isempty(missingFields)
    error( ...
        "event_region_system_records is missing required field(s): %s", ...
        strjoin(missingFields, ", "));
end


if ~isfile(broadContextMat)
    error( ...
        "Canonical D0/Tp MAT not found:\n%s\n" + ...
        "Run run_00_prepare_restoration_shape_phase_ratio_data.m first.", ...
        broadContextMat);
end

B = load( ...
    broadContextMat, ...
    "D0andTp");

if ~isfield(B, "D0andTp")
    error( ...
        "Canonical D0/Tp MAT does not contain D0andTp:\n%s", ...
        broadContextMat);
end

D0andTp = normalizeCanonicalTimingTable(B.D0andTp);

requiredTimingFields = [ ...
    "record_ID", ...
    "T50", ...
    "T80", ...
    "T90", ...
    "T95", ...
    "T100"];

missingTimingFields = requiredTimingFields( ...
    ~ismember(requiredTimingFields, string(D0andTp.Properties.VariableNames)));

if ~isempty(missingTimingFields)
    error( ...
        "D0andTp is missing required field(s): %s", ...
        strjoin(missingTimingFields, ", "));
end

if numel(unique(D0andTp.record_ID)) ~= height(D0andTp)
    error( ...
        "D0andTp contains duplicate record_ID values. " + ...
        "Canonical timing lookup must be one-to-one.");
end


keys = composeKey(records);
uniqueKeys = unique(keys);
pairRows = cell(0, 1);

for i = 1:numel(uniqueKeys)

    index = keys == uniqueKeys(i);
    groupRecords = records(index);

    recE = firstSystemRecord(groupRecords, "Power");
    recW = firstSystemRecord(groupRecords, "Water");
    recG = firstSystemRecord(groupRecords, "Gas");

    if ~isempty(recE) && ~isempty(recW)

        pairRows{end + 1, 1} = buildPairRow( ...
            recE, ...
            recW, ...
            "EW", ...
            D0andTp);
    end

    if ~isempty(recE) && ~isempty(recG)

        pairRows{end + 1, 1} = buildPairRow( ...
            recE, ...
            recG, ...
            "EG", ...
            D0andTp);
    end
end

if isempty(pairRows)
    error("No E->W or E->G paired records were generated.");
end

Paired = struct2table( ...
    vertcat(pairRows{:}));

save( ...
    outMat, ...
    "Paired");

if ~isfile(outMat)
    error( ...
        "Paired timing/coupling MAT was not created:\n%s", ...
        outMat);
end

fprintf( ...
    'Saved paired timing/coupling dataset:\n%s\n', ...
    char(outMat));

fprintf( ...
    'Canonical Tp source:\n%s\n', ...
    char(broadContextMat));

printDatasetSummary(Paired);


[TimingSource, TimingSummary] = buildTimingOutputs(Paired);
writeTwoSheetWorkbook( ...
    sourceFig3e, ...
    "record_level", ...
    TimingSource, ...
    "summary", ...
    TimingSummary);

[CouplingSource, CouplingSummary] = buildCouplingOutputs(Paired);
writeTwoSheetWorkbook( ...
    sourceFig3f, ...
    "record_level", ...
    CouplingSource, ...
    "summary", ...
    CouplingSummary);

fprintf( ...
    'Saved source data:\n%s\n%s\n', ...
    char(sourceFig3e), ...
    char(sourceFig3f));

end


function keys = composeKey(records)
keys = strings(numel(records), 1);
for k = 1:numel(records)
    keys(k) = string(records(k).event_ISO_ID) + "|" + string(records(k).regionIndex);
end
end

function rec = firstSystemRecord(G, systemName)
rec = [];
match = strcmpi(string({G.System}), string(systemName));
nMatch = nnz(match);
if nMatch == 0
    return;
end
if nMatch > 1
    records = G(match);
    recordIDs = string({records.record_ID});
    error( ...
        "Duplicate %s records for event_ISO_ID=%s, regionIndex=%s: %s", ...
        string(systemName), ...
        string(records(1).event_ISO_ID), ...
        string(records(1).regionIndex), ...
        strjoin(recordIDs, ", "));
end
rec = G(match);
end

function row = buildPairRow(recE, recB, pairType, D0andTp)
infoE = canonicalTimingInfo(recE, D0andTp);
infoB = canonicalTimingInfo(recB, D0andTp);

row = struct();
row.pair_ID = canonicalEventTag(recE.event_ISO_ID) + "_R" + string(recE.regionIndex) + "_" + string(pairType);
row.systemPair = string(pairType);
row.event_ISO_ID = recE.event_ISO_ID;
row.event_ID = getFieldOrNaN(recE, "event_ID");
row.regionIndex = recE.regionIndex;
row.regionName = string(getFieldOrDefault(recE, "regionName", ""));
row.record_ID_E = string(recE.record_ID);
row.record_ID_B = string(recB.record_ID);

row.Dmax_E = infoE.Dmax;
row.Dmax_B = infoB.Dmax;

milestones = [50 80 90 95 100];
for i = 1:numel(milestones)
    p = milestones(i);
    fieldName = "T" + string(p);
    TpE = infoE.(fieldName);
    TpB = infoB.(fieldName);

    row.(fieldName + "_E") = TpE;
    row.(fieldName + "_B") = TpB;
    row.("rho" + string(p) + "_B_over_E") = safeRatio(TpB, TpE);
    row.("deltaT" + string(p) + "_B_minus_E") = safeDifference(TpB, TpE);
end

row.CmaxAB90 = computeCmax(recB, infoE, infoB, infoE.T90, infoB.T90);
row.CmaxAB95 = computeCmax(recB, infoE, infoB, infoE.T95, infoB.T95);
row.CmaxAB100 = computeCmax(recB, infoE, infoB, infoE.T100, infoB.T100);
end

function info = canonicalTimingInfo(rec, D0andTp)
recordID = strtrim(string(rec.record_ID));
idx = find(D0andTp.record_ID == recordID);

if isempty(idx)
    error( ...
        "record_ID %s was not found in D0_Tp_broad_context.mat.", ...
        recordID);
end
if numel(idx) > 1
    error( ...
        "record_ID %s appears more than once in D0_Tp_broad_context.mat.", ...
        recordID);
end

info = struct();
info.Dmax = getFieldOrNaN(rec, "PeakOutageFraction");
requireUnitInterval(info.Dmax, "PeakOutageFraction", recordID);
info.PeakTime = getFieldOrNaN(rec, "PeakTime");

for p = [50 80 90 95 100]
    fieldName = "T" + string(p);
    col = D0andTp.(fieldName);
    info.(fieldName) = double(col(idx));
end

allTp = [info.T50 info.T80 info.T90 info.T95 info.T100];
zeroTol = 1e-12;

info.isUnaffected = ...
    (isfinite(info.Dmax) && abs(info.Dmax) <= zeroTol) || ...
    any(isfinite(allTp) & abs(allTp) <= zeroTol);

end

function T = normalizeCanonicalTimingTable(T)
if ~istable(T)
    error("D0andTp must be a table.");
end
if ~ismember("record_ID", string(T.Properties.VariableNames))
    error("D0andTp must contain record_ID.");
end
T.record_ID = strtrim(string(T.record_ID));
end

function val = safeRatio(num, den)
val = NaN;
if isfinite(num) && isfinite(den) && den > 0 && num > 0
    val = num ./ den;
end
end

function val = safeDifference(valueB, valueE)
val = NaN;
if isfinite(valueB) && isfinite(valueE)
    val = valueB - valueE;
end
end

function cmax = computeCmax(recB, infoE, infoB, TpE, TpB)
recoverTol = 1e-8;
zeroTol = 1e-12;
cmax = NaN;

if infoE.isUnaffected
    cmax = 0;
    return;
end

if (isfinite(infoB.Dmax) && abs(infoB.Dmax) <= zeroTol) || ...
        (isfinite(TpB) && abs(TpB) <= zeroTol)
    cmax = 0;
    return;
end

if countDistinctExactFunctionalityValues(recB) < 3
    return;
end

if ~(isfinite(infoB.Dmax) && infoB.Dmax > 0 && ...
        isfinite(TpE) && TpE >= 0 && ...
        isfinite(TpB) && TpB > 0)
    return;
end

if TpE >= TpB - recoverTol
    return;
end

[t, loss, okCurve] = preparedBLossCurve(recB, infoB, recoverTol);
if ~okCurve || max(t) < TpB - recoverTol
    return;
end

[lossAtTpE, okAt] = valueAtOrInterpolate(t, loss, TpE, recoverTol);
if ~okAt
    return;
end

[RLB_to_Tp, okRL] = integrateOnWindow(t, loss, 0, TpB, recoverTol);
if ~(okRL && isfinite(RLB_to_Tp) && RLB_to_Tp > recoverTol)
    return;
end

[x, lossUse, okWindow] = curveOnWindow(t, loss, 0, TpE, recoverTol);
if ~okWindow
    return;
end

numerator = trapz(x, max(lossUse - lossAtTpE, 0));
cmax = numerator ./ RLB_to_Tp;
if ~isfinite(cmax)
    cmax = NaN;
end
end

function [t, loss, ok] = preparedBLossCurve(rec, infoB, recoverTol)
ok = false;
t = [];
loss = [];

Dmax = infoB.Dmax;
peakTime = infoB.PeakTime;
if ~(isfinite(Dmax) && Dmax >= 0 && isfinite(peakTime) && peakTime >= 0)
    return;
end

C = getFieldOrDefault(rec, "functionalityCurve", []);
[t, functionality] = exactCurve(C);
if isempty(t)
    return;
end

loss = 1 - functionality;
requireUnitInterval(loss, "functionality-derived loss", string(rec.record_ID));

keepMask = t > peakTime + recoverTol;
tAfter = t(keepMask);
lossAfter = loss(keepMask);

if peakTime > recoverTol
    t = [0; peakTime; tAfter];
    loss = [Dmax; Dmax; lossAfter];
else
    t = [0; tAfter];
    loss = [Dmax; lossAfter];
end

if any(diff(t) == 0)
    error("Functionality curve contains duplicate day values after peak-time alignment.");
end

ok = numel(t) >= 2 && all(isfinite(t)) && all(isfinite(loss));
end

function n = countDistinctExactFunctionalityValues(rec)
[~, service] = exactCurve(getFieldOrDefault(rec, "functionalityCurve", repmat(makePoint(), 0, 1)));
service = service(isfinite(service));
if isempty(service)
    n = 0;
    return;
end
service = round(service(:), 12);
n = numel(unique(service));
end

function curve = keepExactFiniteCurve(curve)
if isempty(curve)
    curve = repmat(makePoint(), 0, 1);
    return;
end
keep = false(numel(curve), 1);
for i = 1:numel(curve)
    if ~isfield(curve(i), "day_censor_code") || ~isfield(curve(i), "day") || ~isfield(curve(i), "value")
        continue;
    end
    keep(i) = isfinite(curve(i).day) && isfinite(curve(i).value) && double(curve(i).day_censor_code) == 0;
end
curve = curve(keep);
end

function [t, v] = exactCurve(curve)
t = [];
v = [];
curve = keepExactFiniteCurve(curve);
if isempty(curve)
    return;
end
t = [curve.day]';
v = [curve.value]';
[t, ord] = sort(t);
v = v(ord);
if any(diff(t) == 0)
    error("Exact-time curve contains duplicate day values.");
end
end

function [yq, ok] = valueAtOrInterpolate(t, y, xq, tol)
yq = NaN;
ok = false;
if isempty(t) || xq < min(t) - tol || xq > max(t) + tol
    return;
end
idxEq = find(abs(t - xq) <= tol, 1, "first");
if ~isempty(idxEq)
    yq = y(idxEq);
    ok = true;
    return;
end
idx = find(t < xq, 1, "last");
if isempty(idx) || idx >= numel(t)
    return;
end
t1 = t(idx);
t2 = t(idx+1);
y1 = y(idx);
y2 = y(idx+1);
if abs(t2 - t1) <= tol
    return;
end
w = (xq - t1) / (t2 - t1);
yq = y1 + w * (y2 - y1);
ok = isfinite(yq);
end

function [x, y, ok] = curveOnWindow(t, loss, x0, x1, tol)
ok = false;
x = [];
y = [];
if ~(isfinite(x0) && isfinite(x1) && x1 >= x0)
    return;
end
[y0, ok0] = valueAtOrInterpolate(t, loss, x0, tol);
[y1, ok1] = valueAtOrInterpolate(t, loss, x1, tol);
if ~(ok0 && ok1)
    return;
end
mask = t > x0 + tol & t < x1 - tol;
x = [x0; t(mask); x1];
y = [y0; loss(mask); y1];
ok = true;
end

function [area, ok] = integrateOnWindow(t, y, x0, x1, tol)
[xw, yw, ok] = curveOnWindow(t, y, x0, x1, tol);
if ~ok
    area = NaN;
    return;
end
area = trapz(xw, yw);
ok = isfinite(area);
end

function pt = makePoint()
pt = struct("dayRaw", NaN, "day", NaN, "value", NaN, "day_censor_code", NaN);
end

function printDatasetSummary(Paired)
if isempty(Paired)
    fprintf('No E->W/E->G paired records were generated.\n');
    return;
end

fprintf('Paired records generated: %d\n', height(Paired));
for pair = ["EW", "EG"]
    index = Paired.systemPair == pair;
    fprintf( ...
        '%s: total=%d, rho50=%d, rho80=%d, rho90=%d, rho95=%d, rho100=%d, deltaT50=%d, deltaT80=%d, deltaT90=%d, deltaT95=%d, deltaT100=%d, Cmax90=%d, Cmax95=%d, Cmax100=%d\n', ...
        char(pair), ...
        nnz(index), ...
        nnz(index & isfinite(Paired.rho50_B_over_E)), ...
        nnz(index & isfinite(Paired.rho80_B_over_E)), ...
        nnz(index & isfinite(Paired.rho90_B_over_E)), ...
        nnz(index & isfinite(Paired.rho95_B_over_E)), ...
        nnz(index & isfinite(Paired.rho100_B_over_E)), ...
        nnz(index & isfinite(Paired.deltaT50_B_minus_E)), ...
        nnz(index & isfinite(Paired.deltaT80_B_minus_E)), ...
        nnz(index & isfinite(Paired.deltaT90_B_minus_E)), ...
        nnz(index & isfinite(Paired.deltaT95_B_minus_E)), ...
        nnz(index & isfinite(Paired.deltaT100_B_minus_E)), ...
        nnz(index & isfinite(Paired.CmaxAB90)), ...
        nnz(index & isfinite(Paired.CmaxAB95)), ...
        nnz(index & isfinite(Paired.CmaxAB100)));
end
end

function tag = canonicalEventTag(eventISOID)

tag = strtrim(string(eventISOID));

if strlength(tag) == 0 || ismissing(tag)
    tag = "";
    return;
end

if ~startsWith(upper(tag), "E")
    tag = "E" + tag;
end

end

function requireUnitInterval(x, fieldName, recordID)
tol = 1e-8;
x = double(x(:));
x = x(isfinite(x));
if any(x < -tol | x > 1 + tol)
    error( ...
        "%s contains value(s) outside [0,1] for record_ID %s.", ...
        fieldName, ...
        recordID);
end
end

function v = getFieldOrNaN(s, name)
if isfield(s, name)
    v = double(s.(name));
    if isempty(v)
        v = NaN;
    end
else
    v = NaN;
end
end

function v = getFieldOrDefault(s, name, defaultValue)
if isfield(s, name)
    v = s.(name);
    if isempty(v)
        v = defaultValue;
    end
else
    v = defaultValue;
end
end




function [Source, Summary] = buildTimingOutputs(T)

pairs = ["EW", "EG"];
pairLabels = ["E->W", "E->G"];
metrics = ["rho90_B_over_E", "rho95_B_over_E"];
milestones = ["T90", "T95"];

sourceRows = cell(0, 1);
summaryRows = cell(0, 1);

for m = 1:numel(metrics)

    for p = 1:numel(pairs)

        values = double( ...
            T.(char(metrics(m))));

        index = ...
            T.systemPair == pairs(p) & ...
            isfinite(values) & ...
            values > 0;

        nRows = nnz(index);

        if nRows > 0

            pair_ID = string(T.pair_ID(index));
            system_pair = repmat(pairs(p), nRows, 1);
            pair = repmat(pairLabels(p), nRows, 1);
            record_ID_E = string(T.record_ID_E(index));
            record_ID_B = string(T.record_ID_B(index));
            event_ISO_ID = string(T.event_ISO_ID(index));
            event_ID = double(T.event_ID(index));
            regionIndex = double(T.regionIndex(index));
            regionName = string(T.regionName(index));
            milestone = repmat(milestones(m), nRows, 1);
            ratio = values(index);

            sourceRows{end + 1, 1} = table( ...
                pair_ID, ...
                system_pair, ...
                pair, ...
                record_ID_E, ...
                record_ID_B, ...
                event_ISO_ID, ...
                event_ID, ...
                regionIndex, ...
                regionName, ...
                milestone, ...
                ratio);
        end

        q = quantileLocal( ...
            values(index), ...
            [0.05 0.25 0.50 0.75 0.95]);

        milestone = milestones(m);
        system_pair = pairs(p);
        pair = pairLabels(p);
        n = nRows;
        Q05 = q(1);
        Q25 = q(2);
        median = q(3);
        Q75 = q(4);
        Q95 = q(5);

        summaryRows{end + 1, 1} = table( ...
            milestone, ...
            system_pair, ...
            pair, ...
            n, ...
            median, ...
            Q25, ...
            Q75, ...
            Q05, ...
            Q95);
    end
end

if isempty(sourceRows)

    Source = table( ...
        strings(0, 1), ...
        strings(0, 1), ...
        strings(0, 1), ...
        strings(0, 1), ...
        strings(0, 1), ...
        strings(0, 1), ...
        zeros(0, 1), ...
        zeros(0, 1), ...
        strings(0, 1), ...
        strings(0, 1), ...
        nan(0, 1), ...
        'VariableNames', { ...
            'pair_ID', ...
            'system_pair', ...
            'pair', ...
            'record_ID_E', ...
            'record_ID_B', ...
            'event_ISO_ID', ...
            'event_ID', ...
            'regionIndex', ...
            'regionName', ...
            'milestone', ...
            'ratio'});

else

    Source = vertcat(sourceRows{:});
end

Summary = vertcat(summaryRows{:});

end


function [Source, Summary] = buildCouplingOutputs(T)

pairs = ["EW", "EG"];
pairLabels = ["E->W", "E->G"];

sourceRows = cell(0, 1);
summaryRows = cell(0, 1);
valuesAll = double(T.CmaxAB90);

for i = 1:numel(pairs)

    index = ...
        T.systemPair == pairs(i) & ...
        isfinite(valuesAll);

    values = max(0, valuesAll(index));
    nRows = nnz(index);

    if nRows > 0

        pair_ID = string(T.pair_ID(index));
        system_pair = repmat(pairs(i), nRows, 1);
        pair = repmat(pairLabels(i), nRows, 1);
        record_ID_E = string(T.record_ID_E(index));
        record_ID_B = string(T.record_ID_B(index));
        event_ISO_ID = double(T.event_ISO_ID(index));
        event_ID = double(T.event_ID(index));
        regionIndex = double(T.regionIndex(index));
        regionName = string(T.regionName(index));
        CmaxAB90 = values;

        sourceRows{end + 1, 1} = table( ...
            pair_ID, ...
            system_pair, ...
            pair, ...
            record_ID_E, ...
            record_ID_B, ...
            event_ISO_ID, ...
            event_ID, ...
            regionIndex, ...
            regionName, ...
            CmaxAB90);
    end

    q = quantileLocal( ...
        values, ...
        [0.05 0.25 0.50 0.75 0.95]);

    system_pair = pairs(i);
    pair = pairLabels(i);
    n = numel(values);
    median = q(3);
    Q25 = q(2);
    Q75 = q(4);
    Q05 = q(1);
    Q95 = q(5);
    fraction_le_5pct = safeFraction(nnz(values <= 0.05), n);
    fraction_le_10pct = safeFraction(nnz(values <= 0.10), n);

    summaryRows{end + 1, 1} = table( ...
        system_pair, ...
        pair, ...
        n, ...
        median, ...
        Q25, ...
        Q75, ...
        Q05, ...
        Q95, ...
        fraction_le_5pct, ...
        fraction_le_10pct);
end

if isempty(sourceRows)

    Source = table( ...
        strings(0, 1), ...
        strings(0, 1), ...
        strings(0, 1), ...
        strings(0, 1), ...
        strings(0, 1), ...
        zeros(0, 1), ...
        zeros(0, 1), ...
        zeros(0, 1), ...
        strings(0, 1), ...
        nan(0, 1), ...
        'VariableNames', { ...
            'pair_ID', ...
            'system_pair', ...
            'pair', ...
            'record_ID_E', ...
            'record_ID_B', ...
            'event_ISO_ID', ...
            'event_ID', ...
            'regionIndex', ...
            'regionName', ...
            'CmaxAB90'});

else

    Source = vertcat(sourceRows{:});
end

Summary = vertcat(summaryRows{:});

end

function q = quantileLocal(values, probabilities)

values = sort(values(:));
values = values(isfinite(values));

if isempty(values)
    q = nan(size(probabilities));
    return;
end

q = nan(size(probabilities));

for i = 1:numel(probabilities)
    probability = probabilities(i);

    if numel(values) == 1
        q(i) = values(1);
    else
        position = 1 + (numel(values) - 1) * probability;
        lowerIndex = floor(position);
        upperIndex = ceil(position);

        if lowerIndex == upperIndex
            q(i) = values(lowerIndex);
        else
            q(i) = ...
                values(lowerIndex) + ...
                (values(upperIndex) - values(lowerIndex)) * ...
                (position - lowerIndex);
        end
    end
end

end

function value = safeFraction(numerator, denominator)

if denominator > 0
    value = numerator / denominator;
else
    value = NaN;
end

end



function ensureDir(directoryPath)

if exist(char(directoryPath), 'dir') ~= 7
    mkdir(char(directoryPath));
end

end

function writeTwoSheetWorkbook( ...
        filename, ...
        sheet1, ...
        table1, ...
        sheet2, ...
        table2)

if isfile(filename)
    delete(filename);
end

writeTableSheetCompat( ...
    table1, ...
    filename, ...
    sheet1);

writeTableSheetCompat( ...
    table2, ...
    filename, ...
    sheet2);

if ~isfile(filename)
    error( ...
        "Source-data workbook was not created:\n%s", ...
        filename);
end

end

function writeTableSheetCompat( ...
        T, ...
        filename, ...
        sheetName)

C = tableToCellWithHeader(T);
C = sanitizeCell(C);

try

    writecell( ...
        C, ...
        filename, ...
        "Sheet", ...
        char(sheetName), ...
        "Range", ...
        "A1");

catch firstError

    try

        xlswrite( ...
            filename, ...
            C, ...
            char(sheetName), ...
            "A1");

    catch secondError

        error( ...
            [ ...
            'Unable to write source-data workbook: %s\n', ...
            'writecell: %s\n', ...
            'xlswrite: %s'], ...
            char(filename), ...
            firstError.message, ...
            secondError.message);
    end
end

end

function C = tableToCellWithHeader(T)

variableNames = T.Properties.VariableNames;

C = cell( ...
    height(T) + 1, ...
    width(T));

C(1, :) = variableNames;

for c = 1:width(T)

    column = T.(variableNames{c});

    for r = 1:height(T)

        if iscell(column)
            value = column{r};
        else
            value = column(r, :);
        end

        C{r + 1, c} = value;
    end
end

end

function C = sanitizeCell(C)

for r = 1:size(C, 1)

    for c = 1:size(C, 2)

        value = C{r, c};

        if isstring(value)

            value(ismissing(value)) = "";

            if numel(value) > 1
                value = strjoin(value(:).', " ");
            end

            C{r, c} = char(value);

        elseif isnumeric(value)

            if ~isscalar(value)
                C{r, c} = char(strjoin(string(value(:).'), " "));
            elseif ~isfinite(value)
                C{r, c} = [];
            end

        elseif islogical(value)

            C{r, c} = double(value);

        elseif iscategorical(value)

            C{r, c} = char(string(value));
        end
    end
end

end

