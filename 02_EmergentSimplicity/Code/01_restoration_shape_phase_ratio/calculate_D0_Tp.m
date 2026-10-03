function calculate_D0_Tp()

moduleRoot = string(fileparts(mfilename("fullpath")));
emergentSimplicityRoot = string(fileparts(fileparts(moduleRoot)));
analysisRoot = string(fileparts(emergentSimplicityRoot));
inMat = fullfile(analysisRoot, "00_Dataset", "standardized_dataset.mat");
outRoot = fullfile( ...
    emergentSimplicityRoot, ...
    "outputs", ...
    "intermediate", ...
    "01_restoration_shape_phase_ratio");
outMat = fullfile(outRoot, "D0_Tp_broad_context.mat");

P = [0.10:0.10:0.90, 0.95, 1.00];

if ~isfile(inMat)
    error("Input MAT not found: %s", inMat);
end

ensureDir(outRoot);

S = load(inMat);
if ~isfield(S, "event_ISO") || ~isfield(S, "event_region_system_records")
    error("standardized_dataset.mat must contain event_ISO and event_region_system_records.");
end

event_ISO = S.event_ISO;
event_region_system_records = S.event_region_system_records;

D0andTp = buildBroadContextTable(event_region_system_records, event_ISO, P);
save(outMat, "D0andTp", "inMat");

fprintf("Rows: %d\n", height(D0andTp));
fprintf("Rows with finite D0: %d\n", nnz(isfinite(D0andTp.D0)));
for i = 1:numel(P)
    lab = pLabel(P(i));
    fprintf("Rows with finite T%s: %d\n", lab, nnz(isfinite(D0andTp.("T" + lab))));
end
fprintf("Saved:\n%s\n", outMat);
end

function T = buildBroadContextTable(records, event_ISO, P)
rows = cell(numel(records), 1);
for i = 1:numel(records)
    rows{i} = buildOneRow(records(i), event_ISO, P);
end
T = struct2table([rows{:}]');
end

function row = buildOneRow(rec, event_ISO, P)
row = baseRow(rec, event_ISO);
[row.D0, row.D0_source] = computeD0(rec);

[rrCurve, peakTimeUsed] = buildBroadContextRR(rec);

[Tvals, ~] = computeMilestones(rrCurve, peakTimeUsed, P);

if isfinite(row.D0) && abs(row.D0) <= 1e-12 && isSingleDay0NoLoss(rec)
    Tvals(:) = 0;
end
for i = 1:numel(P)
    lab = pLabel(P(i));
    row.("T" + lab) = Tvals(i);
end
end

function tf = isSingleDay0NoLoss(rec)
tf = false;
curveSpecs = { ...
    {"functionalityCurve", 1}, ...
    {"outagefractionCurve", 0}, ...
    {"outagecountCurve", 0}, ...
    {"restorationratioCurve", 1}};
hasAnyCurve = false;

for k = 1:numel(curveSpecs)
    fieldName = curveSpecs{k}{1};
    neutralValue = curveSpecs{k}{2};
    curve = keepExactFiniteCurve(getFieldOrDefault(rec, fieldName, repmat(makePoint(), 0, 1)));
    curve = sortCurveByDay(curve);
    if isempty(curve)
        continue;
    end
    hasAnyCurve = true;
    if numel(curve) ~= 1 || abs(curve(1).day) > 1e-12 || ...
            abs(curve(1).value - neutralValue) > 1e-12
        return;
    end
end

if hasAnyCurve
    peakO = getFieldOrDefault(rec, "PeakOutageCount", NaN);
    peakL = getFieldOrDefault(rec, "PeakOutageFraction", NaN);
    if (isfinite(peakO) && peakO > 1e-12) || (isfinite(peakL) && peakL > 1e-12)

        return;
    end
    tf = true;
    return;
end

peakO = getFieldOrDefault(rec, "PeakOutageCount", NaN);
peakL = getFieldOrDefault(rec, "PeakOutageFraction", NaN);
tf = (isfinite(peakO) && abs(peakO) <= 1e-12) || ...
     (isfinite(peakL) && abs(peakL) <= 1e-12);
end

function row = baseRow(rec, event_ISO)
row = struct();
row.record_ID = string(getFieldOrDefault(rec, "record_ID", ""));
row.event_ISO_ID = getFieldOrDefault(rec, "event_ISO_ID", NaN);
row.event_ID = getFieldOrDefault(rec, "event_ID", NaN);
row.Country = string(getFieldOrDefault(rec, "Country", ""));
row.ISO = string(getFieldOrDefault(rec, "ISO", ""));
row.IncomeGroup = string(getFieldOrDefault(rec, "IncomeGroup", ""));
row.OccTime = string(getFieldOrDefault(rec, "OccTime", ""));
row.magnitude = getFieldOrDefault(rec, "magnitude", NaN);
if ~isfinite(row.magnitude)
    row.magnitude = getFieldOrDefault(rec, "Magnitude", NaN);
end
row.regionIndex = getRecordRegionIndex(rec);
row.regionName = string(getFieldOrDefault(rec, "regionName", ""));
row.System = string(getFieldOrDefault(rec, "System", ""));

row.Depth = NaN;
row.Epicenter = "";
row.EMDATDisNo = "";
row.usgsEventTitle = "";

row.avg_PGA = NaN;
row.population = NaN;
row.popDensity_km2 = NaN;
row.popWeightedPGA_g = NaN;
row.popMeanPGA_g = NaN;
row.pgaMedian_g = NaN;
row.pgaMax_g = NaN;

eid = row.event_ISO_ID;
rid = row.regionIndex;
if isfinite(eid)
    evIdx = findEventIsoIndex(event_ISO, eid);
    if isempty(evIdx)
        error("event_ISO_ID %g from record %s was not found in event_ISO.", eid, char(row.record_ID));
    end
    ev = event_ISO(evIdx);
    if ~isfinite(row.event_ID)
        row.event_ID = getFieldOrDefault(ev, "event_ID", NaN);
    end
    row.Depth = getFieldOrDefault(ev, "Depth", NaN);
    row.Epicenter = string(getFieldOrDefault(ev, "Epicenter", ""));
    row.EMDATDisNo = string(getFieldOrDefault(ev, "EMDATDisNo", ""));
    row.usgsEventTitle = string(getFieldOrDefault(ev, "usgsEventTitle", ""));
    if ~isfinite(row.magnitude)
        row.magnitude = getFieldOrDefault(ev, "magnitude", NaN);
    end
    if ~isfield(ev, "Region") || isempty(ev.Region) || ~isfinite(rid) || rid < 1 || rid > numel(ev.Region)
        return;
    end
    region = ev.Region(rid);
    info = getFieldOrDefault(region, "regionInfo", struct());
    row.avg_PGA = getFieldOrDefault(info, "avg_PGA", NaN);
    row.population = getFieldOrDefault(info, "population", NaN);
    row.popDensity_km2 = getFieldOrDefault(info, "popDensity_km2", NaN);
    row.popWeightedPGA_g = getFieldOrDefault(info, "popWeightedPGA_g", NaN);
    row.popMeanPGA_g = getFieldOrDefault(info, "popMeanPGA_g", NaN);
    row.pgaMedian_g = getFieldOrDefault(info, "pgaMedian_g", NaN);
    row.pgaMax_g = getFieldOrDefault(info, "pgaMax_g", NaN);
end
end

function idx = findEventIsoIndex(event_ISO, eventISOID)
idx = [];
for k = 1:numel(event_ISO)
    candidate = getFieldOrDefault(event_ISO(k), "event_ISO_ID", NaN);
    if isfinite(candidate) && candidate == eventISOID
        idx = k;
        return;
    end
end
end

function [D0, source] = computeD0(rec)
D0 = NaN;
source = "none";

funCurve = sortCurveByDay(getFieldOrDefault(rec, "functionalityCurve", repmat(makePoint(), 0, 1)));
funCurve = keepExactFiniteCurve(funCurve);
if ~isempty(funCurve)
    if isscalar(funCurve)
        d = funCurve(1).day;
        v = funCurve(1).value;
        if abs(d) <= 1e-12
            D0 = 1 - v;
            source = "functionalityCurve_day0_single_point";
            return;
        elseif abs(v - 1) <= 1e-12
        else
            D0 = 1 - v;
            source = "functionalityCurve_single_postevent_point";
            return;
        end
    else
        idx = firstFiniteValueIndex(funCurve);
        if ~isempty(idx)
            D0 = 1 - funCurve(idx).value;
            source = "functionalityCurve_first_point";
            return;
        end
    end
end

peakLoss = getFieldOrDefault(rec, "maxOutageFraction", NaN);
if isfinite(peakLoss)
    D0 = peakLoss;
    source = "maxOutageFraction";
end
end

function [rrCurve, peakTimeUsed] = buildBroadContextRR(rec)
recordID = string(getFieldOrDefault(rec, "record_ID", ""));

rrCurve = keepExactFiniteCurve(getFieldOrDefault(rec, "restorationratioCurve", repmat(makePoint(), 0, 1)));
rrCurve = sortCurveByDay(rrCurve);

if ~isempty(rrCurve)
    rrCurve = requireCurveUnitInterval(rrCurve, recordID, "restorationratioCurve");
    if ~hasValue(rrCurve, 0)
        rrCurve = prependPoint(rrCurve, 0, 0);
        peakTimeUsed = 0;
    else
        peakTimeUsed = firstValueTime(rrCurve, 0);
    end
    rrCurve = truncateAfterFirstOne(rrCurve);
    return;
end

peakO = getFieldOrDefault(rec, "PeakOutageCount", NaN);
peakL = getFieldOrDefault(rec, "PeakOutageFraction", NaN);

orgCurve = keepExactFiniteCurve(getFieldOrDefault(rec, "outagecountCurve", repmat(makePoint(), 0, 1)));
orgCurve = sortCurveByDay(orgCurve);
if isscalar(orgCurve) && isfinite(peakO) && peakO > 0
    currentLoss = requireUnitInterval( ...
        orgCurve(1).value ./ peakO, ...
        recordID, ...
        "outagecountCurve/PeakOutageCount");
    rrVal = requireUnitInterval( ...
        1 - currentLoss, ...
        recordID, ...
        "outagecountCurve-derived restoration ratio");
    rrCurve = [makeScalarPoint(0, 0); makeScalarPoint(orgCurve(1).day, rrVal)];
    rrCurve = truncateAfterFirstOne(sortCurveByDay(rrCurve));
    peakTimeUsed = 0;
    return;
end

funCurve = keepExactFiniteCurve(getFieldOrDefault(rec, "functionalityCurve", repmat(makePoint(), 0, 1)));
funCurve = sortCurveByDay(funCurve);
if isscalar(funCurve) && isfinite(peakL) && peakL > 0
    currentLoss = requireUnitInterval( ...
        1 - funCurve(1).value, ...
        recordID, ...
        "functionalityCurve-derived loss");
    rrVal = requireUnitInterval( ...
        1 - currentLoss ./ peakL, ...
        recordID, ...
        "functionalityCurve/PeakOutageFraction");
    rrCurve = [makeScalarPoint(0, 0); makeScalarPoint(funCurve(1).day, rrVal)];
    rrCurve = truncateAfterFirstOne(sortCurveByDay(rrCurve));
    peakTimeUsed = 0;
    return;
end

if isfinite(peakO) && abs(peakO) <= 1e-12
    rrCurve = makeScalarPoint(0, 1);
    peakTimeUsed = 0;
    return;
end
if isfinite(peakL) && abs(peakL) <= 1e-12
    rrCurve = makeScalarPoint(0, 1);
    peakTimeUsed = 0;
    return;
end

rrCurve = repmat(makePoint(), 0, 1);
peakTimeUsed = NaN;
end

function curve = requireCurveUnitInterval(curve, recordID, sourceName)
for i = 1:numel(curve)
    curve(i).value = requireUnitInterval(curve(i).value, recordID, sourceName);
end
end

function x = requireUnitInterval(x, recordID, sourceName)
tol = 1e-10;
if ~isfinite(x) || x < -tol || x > 1 + tol
    error( ...
        "%s produced a value outside [0,1] for %s: %.15g", ...
        sourceName, ...
        recordID, ...
        x);
end
if x < 0
    x = 0;
elseif x > 1
    x = 1;
end
end
function [Tvals, tauVals] = computeMilestones(rrCurve, peakTimeUsed, P)
Tvals = nan(size(P));
tauVals = nan(size(P));
if isempty(rrCurve)
    return;
end

rrCurve = sortCurveByDay(rrCurve);
times = [rrCurve.day]';
vals = [rrCurve.value]';

for i = 1:numel(P)
    p = P(i);
    if abs(p - 1.0) <= 1e-12
        t = firstValueTimeAtLeast(rrCurve, 1 - 1e-12);
    else
        t = firstCrossingTimeLinear(times, vals, p);
    end
    Tvals(i) = t;
    if isfinite(t) && isfinite(peakTimeUsed)
        tauVals(i) = t - peakTimeUsed;
    end
end
end

function t = firstCrossingTimeLinear(times, vals, p)
t = NaN;
if isscalar(times)
    if vals(1) >= p - 1e-12
        t = times(1);
    end
    return;
end

for i = 1:numel(times)-1
    t1 = times(i);  t2 = times(i+1);
    y1 = vals(i);   y2 = vals(i+1);
    if ~(isfinite(t1) && isfinite(t2) && isfinite(y1) && isfinite(y2))
        continue;
    end
    if y1 >= p - 1e-12
        t = t1;
        return;
    end
    lo = min(y1, y2);
    hi = max(y1, y2);
    if p < lo - 1e-12 || p > hi + 1e-12
        continue;
    end
    if abs(y2 - y1) <= 1e-12
        continue;
    end
    alpha = (p - y1) ./ (y2 - y1);
    if all(alpha >= -1e-12 & alpha <= 1 + 1e-12)
        t = t1 + alpha .* (t2 - t1);
        return;
    end
end

if vals(end) >= p - 1e-12
    t = times(end);
end
end

function tf = hasValue(curve, target)
tf = false;
if isempty(curve)
    return;
end
vals = [curve.value]';
tf = any(isfinite(vals) & abs(vals - target) <= 1e-12);
end

function t = firstValueTime(curve, target)
t = NaN;
if isempty(curve)
    return;
end
curve = sortCurveByDay(curve);
for i = 1:numel(curve)
    if isfinite(curve(i).value) && abs(curve(i).value - target) <= 1e-12
        t = curve(i).day;
        return;
    end
end
end

function t = firstValueTimeAtLeast(curve, thr)
t = NaN;
if isempty(curve)
    return;
end
curve = sortCurveByDay(curve);
for i = 1:numel(curve)
    if isfinite(curve(i).value) && curve(i).value >= thr
        t = curve(i).day;
        return;
    end
end
end

function curve = prependPoint(curve, day, value)
if isempty(curve)
    curve = makeScalarPoint(day, value);
    return;
end
curve = [makeScalarPoint(day, value); curve];
curve = sortCurveByDay(curve);
curve = uniqueByDayKeepFirst(curve);
end

function curveOut = truncateAfterFirstOne(curveIn)
curveOut = curveIn;
if isempty(curveIn)
    return;
end
vals = [curveIn.value]';
idx = find(isfinite(vals) & abs(vals - 1) <= 1e-12, 1, "first");
if ~isempty(idx)
    curveOut = curveIn(1:idx);
end
end

function curve = keepExactFiniteCurve(curveIn)
curve = repmat(makePoint(), 0, 1);
if isempty(curveIn)
    return;
end
buf = repmat(makePoint(), 0, 1);
for i = 1:numel(curveIn)
    d = getFieldOrDefault(curveIn(i), "day", NaN);
    v = getFieldOrDefault(curveIn(i), "value", NaN);
    c = getFieldOrDefault(curveIn(i), "day_censor_code", 0);
    if c == 0 && isfinite(d) && isfinite(v)
        p = curveIn(i);
        p.day = d;
        p.value = v;
        p.day_censor_code = 0;
        buf(end+1,1) = p;
    end
end
curve = uniqueByDayKeepLast(sortCurveByDay(buf));
end

function curve = sortCurveByDay(curve)
if isempty(curve)
    return;
end
[~, idx] = sort([curve.day]);
curve = curve(idx);
end

function curve = uniqueByDayKeepLast(curve)
if isempty(curve)
    return;
end
days = [curve.day]';
keep = true(size(days));
for i = 1:numel(days)-1
    if isfinite(days(i)) && isfinite(days(i+1)) && abs(days(i)-days(i+1)) <= 1e-12
        keep(i) = false;
    end
end
curve = curve(keep);
end

function curve = uniqueByDayKeepFirst(curve)
if isempty(curve)
    return;
end
days = [curve.day]';
keep = true(size(days));
for i = 2:numel(days)
    if isfinite(days(i)) && isfinite(days(i-1)) && abs(days(i)-days(i-1)) <= 1e-12
        keep(i) = false;
    end
end
curve = curve(keep);
end

function idx = firstFiniteValueIndex(curve)
idx = [];
if isempty(curve)
    return;
end
curve = sortCurveByDay(curve);
for i = 1:numel(curve)
    if isfinite(getFieldOrDefault(curve(i), "value", NaN))
        idx = i;
        return;
    end
end
end

function p = makeScalarPoint(day, value)
p = makePoint();
p.dayRaw = day;
p.day = day;
p.value = value;
p.day_censor_code = 0;
end

function p = makePoint()
p = struct('dayRaw', [], 'day', NaN, 'value', NaN, 'day_censor_code', 0);
end

function lab = pLabel(p)
if abs(p - 0.95) <= 1e-12
    lab = "95";
elseif abs(p - 1.00) <= 1e-12
    lab = "100";
else
    lab = string(round(p * 100));
end
end

function v = getRecordRegionIndex(rec)
v = getFieldOrDefault(rec, "regionIndex", NaN);
if isfinite(v)
    return;
end
rid = regexp(char(string(getFieldOrDefault(rec, "record_ID", ""))), '_R(\d+)_', 'tokens', 'once');
if ~isempty(rid)
    v = str2double(rid{1});
end
end

function v = getFieldOrDefault(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName)
    v = s.(fieldName);
    if isempty(v)
        v = defaultValue;
    end
else
    v = defaultValue;
end
end


function ensureDir(folderPath)
if ~isfolder(folderPath)
    mkdir(folderPath);
end
end
