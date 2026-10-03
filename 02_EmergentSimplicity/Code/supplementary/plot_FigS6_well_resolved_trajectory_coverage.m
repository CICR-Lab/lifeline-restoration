
function plot_FigS6_well_resolved_trajectory_coverage(supplementaryRoot, figureDataRoot)

if nargin < 1 || isempty(supplementaryRoot)
    scriptDir = string(fileparts(mfilename("fullpath")));
    moduleRoot = scriptDir;
    emergentRoot = string(fileparts(fileparts(moduleRoot)));
    supplementaryRoot = fullfile(emergentRoot, "outputs", "supplementary");
end
if nargin >= 2 && strlength(string(figureDataRoot)) > 0
    renderFigS6FromSourceData(supplementaryRoot, figureDataRoot);
    return;
end
scriptDir = string(fileparts(mfilename("fullpath")));
moduleRoot = scriptDir;
emergentSimplicityRoot = string(fileparts(fileparts(moduleRoot)));
analysisRoot = string(fileparts(emergentSimplicityRoot));

addpath(scriptDir);

inputMat = fullfile( ...
    analysisRoot, ...
    "00_Dataset", ...
    "standardized_dataset.mat");

suppOutputRoot = supplementaryRoot;

figureDir = fullfile(suppOutputRoot, "figures");
sourceDataDir = fullfile(suppOutputRoot, "source_data");

ensureDir(suppOutputRoot);
ensureDir(figureDir);
ensureDir(sourceDataDir);

if ~isfile(inputMat)
    error("Input MAT not found: %s", inputMat);
end

S = load(inputMat, "event_region_system_records");

if ~isfield(S, "event_region_system_records")
    error( ...
        "Input MAT must contain event_region_system_records: %s", ...
        inputMat);
end

T = scalarTableFromRecords(S.event_region_system_records);
T = normalizeInputTable(T);

requiredVars = [ ...
    "System", ...
    "IncomeGroup", ...
    "avg_PGA", ...
    "magnitude", ...
    "nRestorationRatioPoints", ...
    "nUniqueRestorationValues"];

missingVars = requiredVars( ...
    ~ismember(requiredVars, string(T.Properties.VariableNames)));

if ~isempty(missingVars)
    error( ...
        "event_region_system_records is missing required scalar field(s): %s", ...
        strjoin(missingVars, ", "));
end

T.isWellResolved = ...
    T.nRestorationRatioPoints >= 6 & ...
    T.nUniqueRestorationValues >= 3;

summary = summarizePrimaryVsNonprimary(T);

sourceWorkbook = fullfile( ...
    sourceDataDir, ...
    "SourceData_FigS6.xlsx");

if isfile(sourceWorkbook)
    delete(sourceWorkbook);
end

sourceVars = intersect( ...
    [ ...
        "record_ID", ...
        "event_ID", ...
        "System", ...
        "Country", ...
        "IncomeGroup", ...
        "avg_PGA", ...
        "avg_PGA_plot", ...
        "magnitude", ...
        "eventYear", ...
        "nRestorationRatioPoints", ...
        "nUniqueRestorationValues", ...
        "isWellResolved"], ...
    string(T.Properties.VariableNames), ...
    "stable");

writeTableCompatPortable( ...
    T(:, cellstr(sourceVars)), ...
    sourceWorkbook, ...
    "FigS6_records");

writeTableCompatPortable( ...
    summary, ...
    sourceWorkbook, ...
    "FigS6_summary");

plotFigS6( ...
    T, ...
    fullfile(figureDir, "FigS6_well_resolved_trajectory_coverage"));

reportFigS6PGAComparison(T);

fprintf( ...
    "Saved Fig S6 outputs under:\n%s\n", ...
    suppOutputRoot);


end

function renderFigS6FromSourceData(figureDir, figureDataRoot)

sourceWorkbook = fullfile(figureDataRoot, "SourceData_FigS6.xlsx");
if ~isfile(sourceWorkbook)
    error("Frozen Figure S6 source workbook not found: %s", sourceWorkbook);
end
ensureDir(figureDir);

T = readtable(sourceWorkbook, "Sheet", "FigS6_records", ...
    "TextType", "string", "VariableNamingRule", "preserve");
T = normalizeInputTable(T);
if ~ismember("isWellResolved", string(T.Properties.VariableNames))
    error("Frozen Figure S6 source data do not contain isWellResolved.");
end
T.isWellResolved = logical(double(T.isWellResolved));

plotFigS6(T, fullfile(figureDir, "FigS6_well_resolved_trajectory_coverage"));
fprintf("Saved Fig S6 from frozen figure data under:\n%s\n", figureDir);

end

function T = normalizeInputTable(T)

if ~istable(T)
    error("Input event-region-system records must be convertible to a table.");
end

T = renameVarIfPresent(T, "recordID", "record_ID");

vars = string(T.Properties.VariableNames);

if ismember("System", vars)
    T.System = canonicalSystemName(T.System);
end

for name = ["IncomeGroup", "Country", "ISO"]
    if ismember(name, vars)
        value = string(T.(char(name)));
        value(ismissing(value)) = "";
        T.(char(name)) = value;
    else
        T.(char(name)) = strings(height(T), 1);
    end
end

if ismember("avg_PGA", vars)
    T.avg_PGA = double(T.avg_PGA);
    T.avg_PGA_plot = adjustedPGAForLogScale(T.avg_PGA);
    T.logPGA = safeLog(T.avg_PGA_plot);
else
    T.avg_PGA = nan(height(T), 1);
    T.avg_PGA_plot = nan(height(T), 1);
    T.logPGA = nan(height(T), 1);
end

if ismember("popDensity_km2", vars)
    T.logPopDensity = safeLog(double(T.popDensity_km2));
else
    T.popDensity_km2 = nan(height(T), 1);
    T.logPopDensity = nan(height(T), 1);
end

if ismember("magnitude", vars)
    T.magnitude = double(T.magnitude);
else
    T.magnitude = nan(height(T), 1);
end

T.eventYear = deriveEventYear(T);

finiteYears = T.eventYear(isfinite(T.eventYear));

if isempty(finiteYears)
    yearCenter = 0;
else
    yearCenter = median(finiteYears);
end

finiteMagnitude = T.magnitude(isfinite(T.magnitude));

if isempty(finiteMagnitude)
    magnitudeCenter = 0;
else
    magnitudeCenter = median(finiteMagnitude);
end

T.year_c = T.eventYear - yearCenter;
T.magnitude_c = T.magnitude - magnitudeCenter;

T.sys = categorical( ...
    T.System, ...
    ["Power", "Water", "Gas"]);

T.eqid = eventGroupingFromInput(T);
T.highIncomeCat = highIncomeCategorical(T.IncomeGroup);

end

function y = deriveEventYear(T)

vars = string(T.Properties.VariableNames);
y = nan(height(T), 1);

if ismember("eventYear", vars)

    raw = T.eventYear;

    if isnumeric(raw)
        y = double(raw(:));
        return;
    end
end

if ~ismember("OccTime", vars)
    return;
end

raw = T.OccTime;

if isdatetime(raw)

    y = year(raw);
    return;
end

if isnumeric(raw)

    x = double(raw(:));

    directYear = ...
        isfinite(x) & ...
        x >= 1800 & ...
        x <= 2200;

    y(directYear) = x(directYear);

    other = isfinite(x) & ~directYear;

    if any(other)

        try

            dt = datetime(x(other), "ConvertFrom", "datenum");
            y(other) = year(dt);

        catch
        end
    end

    return;
end

x = string(raw);

for i = 1:numel(x)

    if ismissing(x(i)) || strlength(strtrim(x(i))) == 0
        continue;
    end

    try

        dt = datetime(x(i));

        if ~isnat(dt)
            y(i) = year(dt);
        end

    catch
    end
end

end

function T = scalarTableFromRecords(records)
if isempty(records)
    T = table();
    return;
end
names = fieldnames(records);
keep = false(size(names));
for i = 1:numel(names)
    v = records(1).(names{i});
    keep(i) = isscalar(v) || ischar(v) || isstring(v) || islogical(v) || ...
        (isnumeric(v) && isscalar(v));
end
T = struct2table(rmfield(records, names(~keep)));
end

function y = yearFromOccTime(x)
x = string(x);
y = nan(size(x));
for i = 1:numel(x)
    try
        dt = datetime(x(i), "InputFormat", "yyyy-MM-dd HH:mm:ss");
        if isnat(dt)
            dt = datetime(x(i));
        end
        if ~isnat(dt)
            y(i) = year(dt);
        end
    catch
    end
end
end

function Summary = summarizePrimaryVsNonprimary(T)
vars = ["avg_PGA","magnitude","eventYear"];
rows = {};
for i = 1:numel(vars)
    a = T.(vars(i))(T.isWellResolved == 1);
    b = T.(vars(i))(T.isWellResolved == 0);
    a = a(isfinite(a));
    b = b(isfinite(b));
    qa = quantileOrNaN(a, [0.25 0.50 0.75]);
    qb = quantileOrNaN(b, [0.25 0.50 0.75]);
    rows{end+1,1} = table(vars(i), numel(a), numel(b), ...
        qa(1), qa(2), qa(3), qb(1), qb(2), qb(3), ...
        'VariableNames', {'variable','n_well_resolved','n_lower_resolution', ...
        'Q1_well_resolved','median_well_resolved','Q3_well_resolved', ...
        'Q1_lower_resolution','median_lower_resolution','Q3_lower_resolution'});
end
Summary = vertcat(rows{:});
end

function reportFigS6PGAComparison(T)
hi = T.avg_PGA_plot(T.isWellResolved == 1);
lo = T.avg_PGA_plot(T.isWellResolved == 0);
hi = hi(isfinite(hi) & hi > 0);
lo = lo(isfinite(lo) & lo > 0);

fprintf("\nFig S6 regional PGA comparison\n");
fprintf("--------------------------------------------------\n");
fprintf("Well-resolved cohort:   n = %d | mean = %.4f g | median = %.4f g\n", ...
    numel(hi), mean(hi, "omitnan"), median(hi, "omitnan"));
fprintf("Remaining records cohort: n = %d | mean = %.4f g | median = %.4f g\n", ...
    numel(lo), mean(lo, "omitnan"), median(lo, "omitnan"));
if ~isempty(hi) && ~isempty(lo) && mean(lo, "omitnan") > 0
    fprintf("Mean ratio (well-resolved / remaining records): %.4f\n", ...
        mean(hi, "omitnan") / mean(lo, "omitnan"));
end
fprintf("--------------------------------------------------\n\n");
end

function plotFigS6(T, outBase)
fig = figure("Color", "w", "Units", "pixels", "Position", [80 80 1180 900], "Visible", "off");
set(fig, "DefaultAxesFontName", "Arial", "DefaultTextFontName", "Arial");

positions = [
    0.080 0.565 0.380 0.355
    0.575 0.565 0.380 0.355
    0.080 0.105 0.380 0.355
    0.575 0.105 0.380 0.355];

ax1 = axes(fig, "Position", positions(1,:)); hold(ax1, "on");
drawTwoGroupBoxScatter(ax1, T.avg_PGA_plot, T.isWellResolved, true, "general");
ylabel(ax1, "Service-region mean PGA (g)");
text(ax1, -0.13, 1.06, "a", "Units", "normalized", "FontWeight", "bold", "FontSize", 18);

ax2 = axes(fig, "Position", positions(2,:)); hold(ax2, "on");
drawTwoGroupBoxScatter(ax2, T.magnitude, T.isWellResolved, false, "general");
ylabel(ax2, "Magnitude");
text(ax2, -0.13, 1.06, "b", "Units", "normalized", "FontWeight", "bold", "FontSize", 18);

ax3 = axes(fig, "Position", positions(3,:)); hold(ax3, "on");
drawTwoGroupBoxScatter(ax3, T.eventYear, T.isWellResolved, false, "year");
ylabel(ax3, "Earthquake year");
text(ax3, -0.13, 1.06, "c", "Units", "normalized", "FontWeight", "bold", "FontSize", 18);

ax4 = axes(fig, "Position", positions(4,:)); hold(ax4, "on");
drawIncomeComposition(ax4, T.IncomeGroup, T.isWellResolved);
ylabel(ax4, "Share");
text(ax4, -0.08, 1.07, "d", "Units", "normalized", "FontWeight", "bold", "FontSize", 18);

set(findall(fig, "-property", "FontName"), "FontName", "Arial");
set(findall(fig, "-property", "FontSize"), "FontSize", 9);

finishFigure(fig, outBase, 300);
close(fig);
end

function drawTwoGroupBoxScatter(ax, y, isPrimary, useLog, valueFormat)
groups = {y(isPrimary), y(~isPrimary)};
cols = {[0.17 0.46 0.72], [0.38 0.38 0.38]};
labels = ["Well-resolved", "Remaining records"];
for g = 1:2
    vals = groups{g};
    vals = vals(isfinite(vals));
    if useLog
        vals = vals(vals > 0);
    end
    if isempty(vals), continue; end
    drawManualBox(ax, g, vals, cols{g});
    jitter = 0.075 * sin((1:numel(vals))' * 2.399963);
    scatter(ax, g + jitter, vals, 14, cols{g}, "filled", ...
        "MarkerFaceAlpha", 0.22, "MarkerEdgeColor", "none");
    q = quantile(vals, [0.25 0.50 0.75]);
    text(ax, 0.25 + 0.50*(g-1), 1.025, formatStatLabel(numel(vals), q, valueFormat), ...
        "Units", "normalized", "HorizontalAlignment", "center", "VerticalAlignment", "top", ...
        "FontSize", 8, "FontWeight", "bold", "Color", cols{g}, ...
        "FontName", "Arial", "Interpreter", "tex", "Clipping", "off");
end
set(ax, "XTick", 1:2, "XTickLabel", labels, "FontSize", 8, "TickLength", [0 0]);
if useLog
    set(ax, "YScale", "log");
end
xlim(ax, [0.5 2.5]);
box(ax, "off");
grid(ax, "off");
end

function label = formatStatLabel(n, q, valueFormat)
if valueFormat == "year"
    label = sprintf("{\\it n} = %d\n%.0f [%.0f-%.0f]", n, q(2), q(1), q(3));
else
    label = sprintf("{\\it n} = %d\n%.3g [%.3g-%.3g]", n, q(2), q(1), q(3));
end
end

function drawManualBox(ax, x, vals, col)
q = quantile(vals, [0.05 0.25 0.50 0.75 0.95]);
w = 0.52;
patch(ax, [x-w/2 x+w/2 x+w/2 x-w/2], [q(2) q(2) q(4) q(4)], col, ...
    "FaceAlpha", 0.16, "EdgeColor", col, "LineWidth", 1.2);
plot(ax, [x-w/2 x+w/2], [q(3) q(3)], "-", "Color", col, "LineWidth", 1.6);
plot(ax, [x x], [q(1) q(2)], "-", "Color", [0.20 0.20 0.20], "LineWidth", 0.9);
plot(ax, [x x], [q(4) q(5)], "-", "Color", [0.20 0.20 0.20], "LineWidth", 0.9);
plot(ax, [x-w/5 x+w/5], [q(1) q(1)], "-", "Color", [0.20 0.20 0.20], "LineWidth", 0.9);
plot(ax, [x-w/5 x+w/5], [q(5) q(5)], "-", "Color", [0.20 0.20 0.20], "LineWidth", 0.9);
end

function drawIncomeComposition(ax, incomeGroup, isPrimary)
levels = ["Low income", "Lower middle income", "Upper middle income", "High income"];
shortLabels = ["Low", "Lower middle", "Upper middle", "High"];
cols = [0.38 0.38 0.38; 0.17 0.46 0.72];
shares = nan(2, numel(levels));
nGroup = zeros(2,1);
for g = 1:2
    idx = isPrimary == (g == 1);
    vals = string(incomeGroup(idx));
    vals = vals(strlength(strtrim(vals)) > 0);
    nGroup(g) = numel(vals);
    for j = 1:numel(levels)
        shares(g,j) = mean(vals == levels(j), "omitnan");
    end
end
bh = bar(ax, shares', "grouped", "EdgeColor", "none");
for g = 1:2
    bh(g).FaceColor = cols(g,:);
    bh(g).FaceAlpha = 0.82;
end
ylim(ax, [0 1]);
set(ax, "XTick", 1:numel(levels), "XTickLabel", shortLabels, "FontSize", 8, "TickLength", [0 0]);
legend(ax, ["Well-resolved {\it n} = " + nGroup(1), "Remaining records {\it n} = " + nGroup(2)], ...
    "Location", "northwest", "Box", "off", "FontSize", 8, ...
    "FontName", "Arial", "Interpreter", "tex");
box(ax, "off");
grid(ax, "off");
end

function finishFigure(fig, outBase, resolution)
set(findall(fig, "-property", "FontName"), "FontName", "Arial");
set(findall(fig, "-property", "FontSize"), "FontSize", 9);
set(findall(fig, "Type", "axes"), "LineWidth", 0.8, "Layer", "top");
set(fig, "Visible", "on");
drawnow;
saveas(fig, char(outBase + ".fig"), "fig");
exportgraphics(fig, char(outBase + ".png"), "Resolution", resolution);

try
    exportgraphics(fig, char(outBase + ".pdf"), "ContentType", "vector");
catch
    print(fig, char(outBase + ".pdf"), "-dpdf", "-painters");
end
end

function g = highIncomeCategorical(x)
x = string(x);
g = strings(size(x));
g(x == "High income") = "High income";
g(x ~= "High income" & strlength(x) > 0) = "Non-high income";
g = categorical(g, ["High income","Non-high income"]);
end

function eqid = eventGroupingFromInput(T)
n = height(T);
raw = strings(n,1);
candidateNames = ["event_ID","earthquake_ID","event_ISO_ID","earthquake_ISO_ID"];
for name = candidateNames
    if ~any(strcmp(T.Properties.VariableNames, name))
        continue;
    end
    v = T.(name);
    if isnumeric(v)
        s = strings(size(v));
        good = isfinite(v);
        s(good) = string(v(good));
    else
        s = string(v);
    end
    good = strlength(strtrim(s)) > 0 & ~ismissing(s) & s ~= "<missing>";
    take = (ismissing(raw) | strlength(raw) == 0) & good;
    raw(take) = s(take);
end
recordVar = "";
if any(strcmp(T.Properties.VariableNames, "record_ID"))
    recordVar = "record_ID";
elseif any(strcmp(T.Properties.VariableNames, "recordID"))
    recordVar = "recordID";
end
if strlength(recordVar) > 0
    rid = string(T.(char(recordVar)));
    for i = 1:n
        if (ismissing(raw(i)) || strlength(raw(i)) == 0) && ...
                ~ismissing(rid(i)) && strlength(rid(i)) > 0
            pos = strfind(char(rid(i)), '_');
            if ~isempty(pos)
                raw(i) = extractBefore(rid(i), pos(1));
            else
                raw(i) = rid(i);
            end
        end
    end
end
eqid = categorical(raw);
end

function sys = canonicalSystemName(x)
x = string(x);
sys = strings(size(x));
for i = 1:numel(x)
    switch lower(strtrim(x(i)))
        case {"e","power"}
            sys(i) = "Power";
        case {"w","water"}
            sys(i) = "Water";
        case {"g","gas"}
            sys(i) = "Gas";
        otherwise
            sys(i) = string(x(i));
    end
end
end

function q = quantileOrNaN(x, p)
if isempty(x)
    q = nan(size(p));
else
    q = quantile(x, p);
end
end

function y = safeLog(x)
y = nan(size(x));
mask = isfinite(x) & x > 0;
y(mask) = log(x(mask));
end

function y = adjustedPGAForLogScale(x)
y = double(x);
zeroMask = isfinite(y) & abs(y) <= 1e-12;
y(zeroMask) = 1e-6;
end

function ensureDir(p)

if exist(p, "dir") ~= 7
    mkdir(p);
end

end

function T = renameVarIfPresent(T, oldName, newName)
vars = string(T.Properties.VariableNames);
if any(vars == oldName) && ~any(vars == newName)
    idx = find(vars == oldName, 1, "first");
    T.Properties.VariableNames{idx} = char(newName);
end
end
