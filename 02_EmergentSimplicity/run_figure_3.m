function run_figure_3(varargin)
% Render Figure 3 from packaged or generated source-data workbooks.

analysisRoot = string(fileparts(mfilename("fullpath")));
[dataSource, outputDir, dataRoot] = parseFigureOptions(analysisRoot, varargin{:});
outDir = outputDir;
outBase = fullfile(outDir, "Figure3");

global FONT_AXIS_LABEL FONT_AXIS_TICK FONT_LEGEND FONT_COUNT FONT_PANEL FONT_ANNOT FONT_INSET
BASE_FONT_SIZE = 25.0;
FONT_AXIS_LABEL = 27.0;
FONT_AXIS_TICK = BASE_FONT_SIZE;
FONT_LEGEND = BASE_FONT_SIZE;
FONT_COUNT = BASE_FONT_SIZE;
FONT_PANEL = 27.0;
FONT_ANNOT = BASE_FONT_SIZE;
FONT_INSET = BASE_FONT_SIZE;

if strlength(dataRoot) == 0
    if dataSource == "packaged"
        dataRoot = fullfile(analysisRoot, "figure_data", "main");
    else
        dataRoot = fullfile(analysisRoot, "outputs", "figure_data", "main");
    end
end
sourceFigure3 = fullfile(dataRoot, "SourceData_Figure3.xlsx");

source3b = sourceFigure3;
source3c = sourceFigure3;
source3d = sourceFigure3;
source3e = sourceFigure3;
source3f = sourceFigure3;

requiredFiles = [ ...
    sourceFigure3];

missingFiles = ...
    requiredFiles( ...
        ~isfile(requiredFiles));

if ~isempty(missingFiles)

    error("EmergentFigure:MissingData", ...
        "Required Figure 3 data are missing for DataSource='%s':\n%s", ...
        dataSource, strjoin(missingFiles, newline));
end
ensureDir(outDir);
fprintf("Figure 3 | DataSource=%s | DataRoot=%s\n", dataSource, dataRoot);


sysOrderLong = ["Power","Water","Gas"];
sysNames = ["Electric power","Water supply","Natural gas"];

phaseSysColors = [ ...
    0.70 0.14 0.16
    0.17 0.46 0.72
    0.78 0.55 0.18];

lossSysColors = [ ...
    0.75 0.22 0.25
    0.20 0.50 0.78
    0.86 0.55 0.16];

TrajSource = readSourceSheetContaining( ...
    sourceFigure3, ...
    ["record_ID", "system", "tau_over_tau50", "restoration_ratio"], ...
    ["Fig3a_normalized_trajectories"]);

Phase = readSourceSheetContaining( ...
    source3b, ["system", "metric", "value"], ["Fig3b_record_level"]);
PhaseSummary = readSourceSheetContaining( ...
    source3b, ["metric", "system", "n", "median", "Q25", "Q75", "Q05", "Q95"], ["Fig3b_summary"]);
Phase = normalizeFigure3PhaseSource(Phase);
PhaseSummary = normalizeFigure3PhaseSummary(PhaseSummary);

LossSource = readSourceSheetContaining( ...
    source3c, ["system", "Dmax_tau50", "Lrest", "Lrest_over_Dmax_tau50"], ["Fig3c_record_level"]);
LossSlope = readSourceSheetContaining( ...
    source3c, ["system", "slope_b", "R2"], ["Fig3c_fit_summary"]);
LossRatioSummary = readSourceSheetContaining( ...
    source3c, ["system", "n", "median", "Q25", "Q75", "Q05", "Q95"], ...
    ["Fig3c_ratio_summary"]);
LossSource = normalizeFigure3LossSource(LossSource);
LossSlope = normalizeFigure3LossFit(LossSlope);
LossRatioSummary = normalizeFigure3LossCoefficientSummary(LossRatioSummary);

Tail = readSourceSheetContaining( ...
    source3d, ["system", "eta90"], ["Fig3d_record_level"]);
TailSummary = readSourceSheetContaining( ...
    source3d, ["system", "n", "median", "Q25", "Q75", "Q05", "Q95"], ["Fig3d_summary"]);
Tail = normalizeFigure3TailSource(Tail);
TailSummary = normalizeFigure3TailSummary(TailSummary);

TimingSource = readSourceSheetContaining( ...
    source3e, ["milestone", "system_pair", "pair", "ratio"], ["Fig3e_record_level"]);
TimingSummary = readSourceSheetContaining( ...
    source3e, ["milestone", "system_pair", "pair", "n", "median", "Q25", "Q75", "Q05", "Q95"], ["Fig3e_summary"]);
TimingSource = normalizeFigure3TimingSource(TimingSource);
TimingSummary = normalizeFigure3TimingSummary(TimingSummary);

CouplingSource = readSourceSheetContaining( ...
    source3f, ["system_pair", "pair", "CmaxAB90"], ["Fig3f_record_level"]);
CouplingSummary = readSourceSheetContaining( ...
    source3f, ["system_pair", "pair", "n", "fraction_le_5pct"], ["Fig3f_summary"]);
CouplingSource = normalizeFigure3CouplingSource(CouplingSource);
CouplingSummary = normalizeFigure3CouplingSummary(CouplingSummary);

exportFigure3Main3Row( ...
    outBase, TrajSource, sysOrderLong, sysNames, phaseSysColors, lossSysColors, ...
    Phase, PhaseSummary, LossSource, LossSlope, LossRatioSummary, Tail, TailSummary, ...
    TimingSource, TimingSummary, CouplingSource, CouplingSummary);
fprintf("Saved Figure 3 main figure (curve-a 3x3 layout):\n%s.fig\n%s.png\n%s.pdf\n", ...
    outBase, outBase, outBase);

end

function T = readSourceSheetContaining( ...
        filename, ...
        requiredColumns, ...
        preferredSheets)

filename = string(filename);
requiredColumns = string(requiredColumns);
preferredSheets = string(preferredSheets);

if ~isfile(filename)
    error( ...
        "SourceData workbook not found:\n%s", ...
        filename);
end

try

    sheets = ...
        string(sheetnames(filename));

catch

    [~, sheetsCell] = ...
        xlsfinfo(char(filename));

    sheets = ...
        string(sheetsCell);
end

if isempty(sheets)

    error( ...
        "No worksheets were found in:\n%s", ...
        filename);
end

orderedSheets = strings(0,1);

for i = 1:numel(preferredSheets)

    idx = ...
        find( ...
            strcmpi( ...
                strtrim(sheets), ...
                strtrim(preferredSheets(i))), ...
            1, ...
            "first");

    if ~isempty(idx)

        orderedSheets(end+1,1) = ...
            sheets(idx);
    end
end

for i = 1:numel(sheets)

    if ~any(strcmpi(orderedSheets, sheets(i)))

        orderedSheets(end+1,1) = ...
            sheets(i);
    end
end

for i = 1:numel(orderedSheets)

    try

        candidate = ...
            readtable( ...
                filename, ...
                "Sheet", ...
                char(orderedSheets(i)), ...
                "TextType", ...
                "string", ...
                "VariableNamingRule", ...
                "preserve");

    catch

        continue;
    end

    available = ...
        string(candidate.Properties.VariableNames);

    if all( ...
            ismember( ...
                requiredColumns, ...
                available))

        T = candidate;
        return;
    end
end

error( ...
    "No worksheet in %s contains required column(s): %s", ...
    filename, ...
    strjoin(requiredColumns, ", "));

end

function T = normalizeFigure3PhaseSource(T)
T = normalizeStrings(T);
T.system = normalizeSystemName(T.system);
T.metric = normalizePhaseMetric(T.metric);
T.value = double(T.value);
if ismember("record_ID", string(T.Properties.VariableNames))
    T.curveID = string(T.record_ID);
elseif ~ismember("curveID", string(T.Properties.VariableNames))
    T.curveID = "curve_" + string((1:height(T))');
end
end

function T = normalizeFigure3PhaseSummary(T)
T = normalizeStrings(T);
T.metric = normalizePhaseMetric(T.metric);
system = string(T.system);
for i = 1:numel(system)
    if any(lower(strtrim(system(i))) == ["all systems", "pooled", "all"])
        system(i) = "All systems";
    else
        system(i) = normalizeSystemName(system(i));
    end
end
T.system = system;
T.n = double(T.n); T.median = double(T.median);
T.Q25 = double(T.Q25); T.Q75 = double(T.Q75);
T.Q05 = double(T.Q05); T.Q95 = double(T.Q95);
end

function metric = normalizePhaseMetric(metric)
metric = lower(strtrim(string(metric)));
metric(metric == "tau80_over_tau50" | metric == "tau80/tau50") = "tau80/tau50";
metric(metric == "tau90_over_tau80" | metric == "tau90/tau80") = "tau90/tau80";
metric(metric == "tau95_over_tau90" | metric == "tau95/tau90") = "tau95/tau90";
end

function T = normalizeFigure3LossSource(T)
T = normalizeStrings(T);
T.system = normalizeSystemName(T.system);
T.Dmax_tau50 = double(T.Dmax_tau50);
T.Lrest = double(T.Lrest);
T.Lrest_over_Dmax_tau50 = double(T.Lrest_over_Dmax_tau50);
end

function T = normalizeFigure3LossFit(T)

T = normalizeStrings(T);
system = string(T.system);
sl = lower(strtrim(system));
system(sl == "pooled" | sl == "all" | sl == "all systems") = "All systems";
idx = system ~= "All systems";
system(idx) = normalizeSystemName(system(idx));
T.system = system;

T.slope = ...
    double(T.slope_b);

vars = ...
    string(T.Properties.VariableNames);

if ismember( ...
        "bootstrap_slope_ci95_low", ...
        vars)

    T.slope_ci95_low = ...
        double(T.bootstrap_slope_ci95_low);

elseif ismember( ...
        "slope_ci95_low", ...
        vars)

    T.slope_ci95_low = ...
        double(T.slope_ci95_low);

else

    error( ...
        "Fig. 3c fit_summary is missing a lower 95%% slope CI.");
end

if ismember( ...
        "bootstrap_slope_ci95_high", ...
        vars)

    T.slope_ci95_high = ...
        double(T.bootstrap_slope_ci95_high);

elseif ismember( ...
        "slope_ci95_high", ...
        vars)

    T.slope_ci95_high = ...
        double(T.slope_ci95_high);

else

    error( ...
        "Fig. 3c fit_summary is missing an upper 95%% slope CI.");
end

T.R2 = ...
    double(T.R2);

end

function T = normalizeFigure3LossCoefficientSummary(T)
T = normalizeStrings(T);
s = string(T.system);
sl = lower(strtrim(s));
s(sl == "pooled" | sl == "all" | sl == "all systems") = "All systems";
idx = s ~= "All systems";
if any(idx)
    s(idx) = normalizeSystemName(s(idx));
end
T.system = s;
T.n = double(T.n);
T.median = double(T.median);
T.Q25 = double(T.Q25);
T.Q75 = double(T.Q75);
T.Q05 = double(T.Q05);
T.Q95 = double(T.Q95);
end

function T = normalizeFigure3TailSource(T)
T = normalizeStrings(T);
T.system = normalizeSystemName(T.system);
T.value = double(T.eta90);
end

function T = normalizeFigure3TailSummary(T)
T = normalizeStrings(T);
T.system = normalizeSystemName(T.system);
T.n = double(T.n); T.median = double(T.median);
T.Q25 = double(T.Q25); T.Q75 = double(T.Q75);
T.Q05 = double(T.Q05); T.Q95 = double(T.Q95);
end

function T = normalizeFigure3TimingSource(T)
T = normalizeStrings(T);
T.system_pair = upper(strtrim(string(T.system_pair)));
T.pair = normalizePairLabel(T.pair);
T.milestone = upper(strtrim(string(T.milestone)));
T.ratio = double(T.ratio);
end

function T = normalizeFigure3TimingSummary(T)
T = normalizeStrings(T);
T.system_pair = upper(strtrim(string(T.system_pair)));
T.pair = normalizePairLabel(T.pair);
T.milestone = upper(strtrim(string(T.milestone)));
T.n = double(T.n); T.median = double(T.median);
T.Q25 = double(T.Q25); T.Q75 = double(T.Q75);
T.Q05 = double(T.Q05); T.Q95 = double(T.Q95);
end

function T = normalizeFigure3CouplingSource(T)
T = normalizeStrings(T);
T.system_pair = upper(strtrim(string(T.system_pair)));
T.pair = normalizePairLabel(T.pair);
T.CmaxAB90 = double(T.CmaxAB90);
end

function T = normalizeFigure3CouplingSummary(T)
T = normalizeStrings(T);
T.system_pair = upper(strtrim(string(T.system_pair)));
T.pair = normalizePairLabel(T.pair);
T.n = double(T.n);
T.fraction_le_5pct = double(T.fraction_le_5pct);
if ismember("fraction_le_10pct", string(T.Properties.VariableNames))
    T.fraction_le_10pct = double(T.fraction_le_10pct);
end
end

function pair = normalizePairLabel(pair)
pair = strtrim(string(pair));
pair(pair == "W/E") = "E->W";
pair(pair == "G/E") = "E->G";
end


function exportFigure3Main3Row(outBase, TrajSource, sysOrderLong, sysNames, phaseSysColors, lossSysColors, ...
    Phase, PhaseSummary, LossSource, LossSlope, LossRatioSummary, Tail, TailSummary, ...
    TimingSource, TimingSummary, CouplingSource, CouplingSummary)

global FONT_AXIS_LABEL FONT_AXIS_TICK FONT_LEGEND FONT_COUNT FONT_PANEL FONT_ANNOT FONT_INSET
fontBackup = [FONT_AXIS_LABEL, FONT_AXIS_TICK, FONT_LEGEND, FONT_COUNT, FONT_PANEL, FONT_ANNOT, FONT_INSET];
baseFontSize = 25.0;
axisLabelFontSize = 27.0;
panelFontSize = 27.0;
FONT_AXIS_LABEL = axisLabelFontSize; FONT_AXIS_TICK = baseFontSize;
FONT_COUNT = baseFontSize; FONT_LEGEND = baseFontSize;
FONT_PANEL = panelFontSize; FONT_ANNOT = baseFontSize; FONT_INSET = baseFontSize;

figWidth = 2140; figHeight = 2180;
fig = figure("Color", "w", "Units", "pixels", "Position", [35 25 figWidth figHeight], ...
    "Renderer", "painters", "Visible", "off", "Resize", "off");
setFigureDefaults(fig);

leftX = 0.060; leftW = 0.275;
rightX1 = 0.400; rightW = 0.230; rightGap = 0.090;
rightX2 = rightX1 + rightW + rightGap; rightW2 = 0.230;
rowH = 0.230; yTop = 0.705; yMid = 0.385; yBot = 0.065;

axA1 = axes(fig, "Position", [leftX yTop leftW rowH]); hold(axA1, "on");
axA2 = axes(fig, "Position", [leftX yMid leftW rowH]); hold(axA2, "on");
axA3 = axes(fig, "Position", [leftX yBot leftW rowH]); hold(axA3, "on");

axB = axes(fig, "Position", [rightX1 yTop rightW + rightGap + rightW2 rowH]); hold(axB, "on");
axC = axes(fig, "Position", [rightX1 yMid rightW rowH]); hold(axC, "on");
axD = axes(fig, "Position", [rightX2 yMid rightW2 rowH]); hold(axD, "on");

yBotRight = yBot;
axE = axes(fig, "Position", [rightX1 yBotRight rightW rowH]); hold(axE, "on");
axF = axes(fig, "Position", [rightX2 yBotRight rightW2 rowH]); hold(axF, "on");

drawNormalizedTrajectoriesThreeSystems( ...
    [axA1 axA2 axA3], ...
    TrajSource, ...
    sysOrderLong, ...
    sysNames, ...
    phaseSysColors);


for ax = [axA1 axA2 axA3]
    set(ax, "YTick", 0:0.2:1, "YTickLabel", compose("%.1f", 0:0.2:1));
end
xlabA1 = xlabel(axA1, "$\tau/\tau_{50}$", "FontSize", FONT_AXIS_LABEL, "Interpreter", "latex", "FontName", "Arial");
xlabA2 = xlabel(axA2, "$\tau/\tau_{50}$", "FontSize", FONT_AXIS_LABEL, "Interpreter", "latex", "FontName", "Arial");
xlabA3 = xlabel(axA3, "$\tau/\tau_{50}$", "FontSize", FONT_AXIS_LABEL, "Interpreter", "latex", "FontName", "Arial");
set([xlabA1 xlabA2 xlabA3], "Units", "normalized", "VerticalAlignment", "top");

posA1 = get(xlabA1, "Position"); posA1(2) = -0.08; set(xlabA1, "Position", posA1);
posA2 = get(xlabA2, "Position"); posA2(2) = -0.08; set(xlabA2, "Position", posA2);
posA3 = get(xlabA3, "Position"); posA3(2) = -0.08; set(xlabA3, "Position", posA3);
ylabel(axA1, "$R(\tau)$", ...
    "FontSize", FONT_AXIS_LABEL, "Interpreter", "latex");
ylabel(axA2, "$R(\tau)$", ...
    "FontSize", FONT_AXIS_LABEL, "Interpreter", "latex");
ylabel(axA3, "$R(\tau)$", ...
    "FontSize", FONT_AXIS_LABEL, "Interpreter", "latex");

drawPhaseRatiosPanel(axB, Phase, PhaseSummary, sysOrderLong, phaseSysColors);
drawLossCompressionPanelFromSource(fig, axC, LossSource, LossSlope, LossRatioSummary, sysOrderLong, lossSysColors);
drawTailLossPanel(axD, Tail, TailSummary, sysOrderLong, sysNames, lossSysColors);
drawTimingRatioPanel(axE, TimingSource, TimingSummary, true);
drawCmaxDistributionPanel(axF, CouplingSource, CouplingSummary);

unifyXAxisFont([axA1 axA2 axA3 axB axC axD axE axF], FONT_AXIS_TICK);

addPanelLetter(fig, axA1, "a");
addPanelLetter(fig, axB, "b");
addPanelLetter(fig, axC, "c", [], 0.010);
addPanelLetter(fig, axD, "d", [], 0.010);
addPanelLetter(fig, axE, "e", [], 0.012);
addPanelLetter(fig, axF, "f", [], 0.012);

turnAllGridsOff(fig);
hideAxesToolbars(fig);
applyUniformFigureTypography(fig, baseFontSize, axisLabelFontSize, panelFontSize);
drawnow;

exportFixedCanvasPng(fig, outBase + ".png", figWidth, figHeight);
try
    exportgraphics(fig, outBase + ".pdf", "ContentType", "vector", "BackgroundColor", "white");
catch
    print(fig, char(outBase + ".pdf"), "-dpdf", "-painters");
end
prepareFigureForVisibleFigSave(fig, figWidth, figHeight);
savefig(fig, char(outBase + ".fig"), "-v7");
close(fig);

FONT_AXIS_LABEL = fontBackup(1); FONT_AXIS_TICK = fontBackup(2); FONT_LEGEND = fontBackup(3);
FONT_COUNT = fontBackup(4); FONT_PANEL = fontBackup(5); FONT_ANNOT = fontBackup(6); FONT_INSET = fontBackup(7);

end

function drawNormalizedTrajectoriesThreeSystems(axList, T, sysOrder, sysNames, sysColors)

global FONT_AXIS_LABEL FONT_AXIS_TICK FONT_COUNT

requiredVars = ["record_ID", "system", "tau_over_tau50", "restoration_ratio"];
vars = string(T.Properties.VariableNames);
missing = requiredVars(~ismember(requiredVars, vars));
if ~isempty(missing)
    error("Fig. 3a source data is missing required column(s): %s", strjoin(missing, ", "));
end

recordID = string(T.record_ID);
system = normalizeSystemName(T.system);
uAll = double(T.tau_over_tau50);
rAll = double(T.restoration_ratio);
gridU = linspace(0, 4, 161);

for s = 1:numel(sysOrder)
    ax = axList(s);
    hold(ax, "on");

    idxSys = system == sysOrder(s) & isfinite(uAll) & isfinite(rAll) & ...
        uAll >= 0 & rAll >= 0 & rAll <= 1 & strlength(recordID) > 0;

    ids = unique(recordID(idxSys), "stable");
    Y = nan(numel(ids), numel(gridU));
    curveCount = 0;

    lineColor = sysColors(s,:);
    lightColor = blendWithWhite(lineColor, 0.65);

    for i = 1:numel(ids)
        idx = idxSys & recordID == ids(i);
        u = uAll(idx);
        r = rAll(idx);
        [u, order] = sort(u(:));
        r = r(order);
        [u, keepUnique] = unique(u, "stable");
        r = r(keepUnique);

        if numel(u) < 2
            continue;
        end

        plot(ax, u, r, "-", "Color", lightColor, "LineWidth", 0.42, ...
            "HandleVisibility", "off");

        validGrid = gridU >= min(u) & gridU <= max(u);
        if any(validGrid)
            Y(i, validGrid) = interp1(u, r, gridU(validGrid), "linear");
        end
        curveCount = curveCount + 1;
    end

    medY = median(Y, 1, "omitnan");
    q10 = pointwiseQuantile(Y, 0.10);
    q90 = pointwiseQuantile(Y, 0.90);
    finiteBand = isfinite(q10) & isfinite(q90);

    if any(finiteBand)
        xBand = gridU(finiteBand);
        fill(ax, [xBand fliplr(xBand)], [q10(finiteBand) fliplr(q90(finiteBand))], ...
            lineColor, "FaceAlpha", 0.15, "EdgeColor", "none", "HandleVisibility", "off");
    end

    plot(ax, gridU, medY, "-", "Color", lineColor, "LineWidth", 2.6, ...
        "HandleVisibility", "off");
    xline(ax, 1, ":", "Color", [0.38 0.38 0.38], "LineWidth", 1.15);
    yline(ax, 0.5, ":", "Color", [0.38 0.38 0.38], "LineWidth", 1.15);

    xlim(ax, [0 4]);
    ylim(ax, [0 1]);
    set(ax, "FontSize", FONT_AXIS_TICK, "TickLength", [0 0]);
    xlabel(ax, "$\tau/\tau_{50}$", "FontSize", FONT_AXIS_LABEL, ...
        "Interpreter", "latex", "FontName", "Arial");

    if s == 1
        ylabel(ax, "$R(\tau)$", "FontSize", FONT_AXIS_LABEL, "Interpreter", "latex");
    else
        set(ax, "YTickLabel", []);
    end

    text(ax, 0.500, 1.025, char(sysNames(s)), "Units", "normalized", ...
        "HorizontalAlignment", "center", "VerticalAlignment", "bottom", ...
        "FontSize", FONT_COUNT, "Color", [0.10 0.10 0.10], ...
        "FontName", "Arial", "Interpreter", "none", "Clipping", "off");

    sysSubscripts = ["E", "W", "G"];
    text(ax, 0.965, 0.035, sprintf("$n_{%s}=%d$", sysSubscripts(s), curveCount), ...
        "Units", "normalized", "HorizontalAlignment", "right", ...
        "VerticalAlignment", "bottom", "FontSize", FONT_COUNT, ...
        "Color", [0.10 0.10 0.10], "FontName", "Arial", "Interpreter", "latex");

    box(ax, "off");
    grid(ax, "off");
end

linkaxes(axList, "xy");

end
function c = blendWithWhite(baseColor, amount)
amount = max(0, min(1, amount));
c = amount * [1 1 1] + (1 - amount) * baseColor;
end

function q = pointwiseQuantile(Y, p)
q = nan(1, size(Y, 2));
for j = 1:size(Y, 2)
    x = Y(:, j);
    x = x(isfinite(x));
    if ~isempty(x)
        q(j) = quantile(x, p);
    end
end
end

function T = normalizeStrings(T)
vars = string(T.Properties.VariableNames);
for i = 1:numel(vars)
    if iscellstr(T.(vars(i))) || isstring(T.(vars(i))) || ischarColumn(T.(vars(i)))
        T.(vars(i)) = string(T.(vars(i)));
    end
end
end

function tf = ischarColumn(v)
tf = iscell(v) && (isempty(v) || all(cellfun(@(x) ischar(x) || isstring(x), v)));
end

function s = normalizeSystemName(s)
s = string(s);
sl = lower(strtrim(s));

s(sl == "e" | ...
    sl == "electric power" | ...
  sl == "power") = "Power";

s(sl == "w" | ...
  sl == "water" | ...
  sl == "water supply") = "Water";

s(sl == "g" | ...
  sl == "gas" | ...
  sl == "natural gas") = "Gas";
end

function drawPhaseRatiosPanel(ax,R,Summary,sysOrder,sysColors)
global FONT_AXIS_LABEL FONT_AXIS_TICK FONT_COUNT
metrics=["tau80/tau50","tau90/tau80","tau95/tau90"];
ratioLabels=["$\tau_{80}/\tau_{50}$","$\tau_{90}/\tau_{80}$","$\tau_{95}/\tau_{90}$"];

systemLabels=["All","E","W","G"];
pooledColor = [0.42 0.42 0.42];
off=[-0.30 -0.10 0.10 0.30];

phaseN=zeros(numel(metrics),numel(sysOrder)+1); hold(ax,"on");

for m=1:numel(metrics)

    pooledVals = R.value(ismember(R.system,sysOrder) & R.metric==metrics(m));
    pooledVals = pooledVals(isfinite(pooledVals));
    [qPooled,nPooled]=getPhaseSummaryValues(Summary,metrics(m),"All systems");
    phaseN(m,1) = nPooled;
    if ~isempty(pooledVals)
        x0 = m + off(1);
      
        for s = 1:numel(sysOrder)
            valsSys = R.value(R.system==sysOrder(s) & R.metric==metrics(m));
            valsSys = valsSys(isfinite(valsSys));
            if isempty(valsSys), continue; end
            jitterSys = 0.060*sin((1:numel(valsSys))'*(2.399963 + 0.17*s));
            scatter(ax,x0+jitterSys,valsSys,18,sysColors(s,:),"filled", ...
                "MarkerFaceAlpha",0.82,"MarkerEdgeColor","none");
        end
        if all(isfinite(qPooled)), drawBoxDot(ax,x0,qPooled,pooledColor,0.085); end
    end


    for s=1:numel(sysOrder)
        vals=R.value(R.system==sysOrder(s) & R.metric==metrics(m));
        vals=vals(isfinite(vals));
        [q,n]=getPhaseSummaryValues(Summary,metrics(m),sysOrder(s));
        phaseN(m,s+1)=n;
        if isempty(vals) || any(~isfinite(q)), continue; end
        x0=m+off(s+1);
        jitter=0.060*sin((1:numel(vals))'*2.399963);
        scatter(ax,x0+jitter,vals,18,sysColors(s,:),"filled", ...
            "MarkerFaceAlpha",0.82,"MarkerEdgeColor","none");
        drawBoxDot(ax,x0,q,sysColors(s,:),0.085);
    end


    nY1 = 5.10;
    nY2 = 4.78;
    eqX1 = m - 0.25;
    eqX2 = m + 0.25;
    drawAlignedNEntry(ax,eqX1,nY1,'All',phaseN(m,1),FONT_COUNT);
    drawAlignedNEntry(ax,eqX2,nY1,'W',  phaseN(m,3),FONT_COUNT);
    drawAlignedNEntry(ax,eqX1,nY2,'E',  phaseN(m,2),FONT_COUNT);
    drawAlignedNEntry(ax,eqX2,nY2,'G',  phaseN(m,4),FONT_COUNT);
end


allTicks=[];
for m=1:numel(metrics)
    allTicks=[allTicks, m+off];
end
set(ax,"XTick",allTicks,"XTickLabel",repmat({''},1,numel(allTicks)), ...
    "TickLabelInterpreter","none","FontSize",FONT_AXIS_TICK,"TickLength",[0 0]);

ySystem = -0.01;
yRatio  = -0.08;
for m=1:numel(metrics)
    for s=1:numel(systemLabels)
        xNormSystem = (m + off(s) - 0.52) / (3.48 - 0.52);
        text(ax,xNormSystem,ySystem,systemLabels(s), ...
            "Units","normalized", ...
            "HorizontalAlignment","center","VerticalAlignment","top", ...
            "FontSize",FONT_AXIS_TICK,"Interpreter","none","Clipping","off");
    end
    xNorm = (m - 0.52) / (3.48 - 0.52);
    text(ax,xNorm,yRatio,ratioLabels(m), ...
        "Units","normalized", ...
        "HorizontalAlignment","center","VerticalAlignment","top", ...
        "FontSize",FONT_AXIS_TICK,"Interpreter","latex","Clipping","off", ...
        "Tag","PanelBRatioLabel");
end

ylabel(ax,"Milestone ratio","FontSize",FONT_AXIS_LABEL);
ylim(ax,[0.75 4.85]);
xlim(ax,[0.52 3.48]);
box(ax,"off");
end

function drawAlignedNEntry(ax,xEq,y,label,nValue,fontSize)
if strcmpi(string(label), "All")
    prefix = '$n$';
    prefixRightShift = -0.040;
else
    prefix = sprintf('$n_{%s}$', char(string(label)));
    prefixRightShift = 0;
end

valueText = sprintf('$%d$', round(nValue));

gapLeft = 0.045;
gapRight = 0.055;
text(ax,xEq-gapLeft+prefixRightShift,y,prefix, ...
    "HorizontalAlignment","right","VerticalAlignment","bottom", ...
    "FontSize",fontSize,"Interpreter","latex","Clipping","off");
text(ax,xEq,y,'$=$', ...
    "HorizontalAlignment","center","VerticalAlignment","bottom", ...
    "FontSize",fontSize,"Interpreter","latex","Clipping","off");
text(ax,xEq+gapRight,y,valueText, ...
    "HorizontalAlignment","left","VerticalAlignment","bottom", ...
    "FontSize",fontSize,"Interpreter","latex","Clipping","off");
end

function [q,n]=getPhaseSummaryValues(Summary,metric,system)
idx=Summary.metric==metric & Summary.system==system;
if nnz(idx)~=1, error("Expected one Fig. 3b summary row for metric '%s' and system '%s'; found %d.",metric,system,nnz(idx)); end
row=Summary(idx,:); q=[double(row.Q05),double(row.Q25),double(row.median),double(row.Q75),double(row.Q95)]; n=double(row.n);
end

function drawLossCompressionPanelFromSource(fig, ax, D, Slope, RatioSummary, sysOrder, sysColors)
global FONT_AXIS_LABEL FONT_AXIS_TICK FONT_ANNOT FONT_INSET FONT_COUNT
hold(ax, "on");
nBySys = zeros(numel(sysOrder), 1);
for s = 1:numel(sysOrder)
    idx = D.system == sysOrder(s);
    nBySys(s) = nnz(idx);
    scatter(ax, double(D.Dmax_tau50(idx)), double(D.Lrest(idx)), 34, sysColors(s,:), "filled", ...
        "MarkerFaceAlpha", 0.85, "MarkerEdgeColor", "none");
end

allRow = Slope(string(Slope.system) == "All systems", :);
if height(allRow) ~= 1
    error("Expected exactly one All systems Fig. 3c fit row; found %d.", height(allRow));
end
    b = double(allRow.slope(1));
    blo = double(allRow.slope_ci95_low(1));
    bhi = double(allRow.slope_ci95_high(1));
    r2 = double(allRow.R2(1));
    x = double(D.Dmax_tau50);
    x = x(isfinite(x) & x > 0);
    if ~isempty(x) && all(isfinite([b blo bhi]))
        xx = logspace(log10(min(x)), log10(max(x)), 260)';
        fill(ax, [xx; flipud(xx)], [blo*xx; flipud(bhi*xx)], [0 0 0], ...
            "FaceAlpha", 0.14, "EdgeColor", "none");
        plot(ax, xx, b*xx, "k-", "LineWidth", 2.0);
        annotationText = sprintf('$b = %.2f,\\; R^2 = %.3f$', b, r2);
        text(ax, 0.055, 0.955, annotationText, "Units", "normalized", ...
            "HorizontalAlignment", "left", "VerticalAlignment", "top", "FontSize", FONT_COUNT, ...
            "FontWeight", "normal", "Color", [0.10 0.10 0.10], "BackgroundColor", "w", "Margin", 1.0, ...
            "Interpreter", "latex");
    end

set(ax, "XScale", "log", "YScale", "log", "FontSize", FONT_AXIS_TICK);
xlimLogData(ax, double(D.Dmax_tau50));
ylimLogData(ax, double(D.Lrest));
setLogTicksForVisibleRange(ax, "x", 3);
setLogTicksForVisibleRange(ax, "y", 4);
xlabel(ax, "$D_{\mathrm{max}}\tau_{50}$", "FontSize", FONT_AXIS_LABEL, "Interpreter", "latex", "FontName", "Arial");
ylabel(ax, "$L_{\mathrm{rest}}$", "FontSize", FONT_AXIS_LABEL, "Interpreter", "latex");

countX = linspace(0.065,0.86,4);
countStrings = {sprintf('$n=%d$',height(D)), ...
    sprintf('$n_E=%d$',nBySys(1)), ...
    sprintf('$n_W=%d$',nBySys(2)), ...
    sprintf('$n_G=%d$',nBySys(3))};
for ii = 1:4
    text(ax,countX(ii),1.018,countStrings{ii},"Units","normalized", ...
        "HorizontalAlignment","center","VerticalAlignment","bottom", ...
        "FontSize",FONT_COUNT,"Interpreter","latex","Clipping","off");
end
legendX = 0.080; legendY = [0.74 0.665 0.59];
text(ax,legendX,legendY(1),char(9679),"Units","normalized","HorizontalAlignment","center", ...
    "VerticalAlignment","middle","FontSize",FONT_COUNT+1.0,"Color",sysColors(1,:),"Clipping","off");
text(ax,legendX,legendY(2),char(9679),"Units","normalized","HorizontalAlignment","center", ...
    "VerticalAlignment","middle","FontSize",FONT_COUNT+1.0,"Color",sysColors(2,:),"Clipping","off");
text(ax,legendX,legendY(3),char(9679),"Units","normalized","HorizontalAlignment","center", ...
    "VerticalAlignment","middle","FontSize",FONT_COUNT+1.0,"Color",sysColors(3,:),"Clipping","off");
text(ax,legendX+0.050,legendY(1),'E',"Units","normalized","HorizontalAlignment","left", ...
    "VerticalAlignment","middle","FontSize",FONT_COUNT,"Interpreter","none");
text(ax,legendX+0.050,legendY(2),'W',"Units","normalized","HorizontalAlignment","left", ...
    "VerticalAlignment","middle","FontSize",FONT_COUNT,"Interpreter","none");
text(ax,legendX+0.050,legendY(3),'G',"Units","normalized","HorizontalAlignment","left", ...
    "VerticalAlignment","middle","FontSize",FONT_COUNT,"Interpreter","none");
box(ax, "off");
addBiBoxInsetFromSource(fig, ax, D, RatioSummary, sysOrder, sysColors);
end

function addBiBoxInsetFromSource(fig, axMain, D, RatioSummary, sysOrder, sysColors)
global FONT_INSET
mainPos = get(axMain, "Position");
axInset = axes(fig, "Position", [mainPos(1)+mainPos(3)*0.685, ...
    mainPos(2)+mainPos(4)*0.075, mainPos(3)*0.270, mainPos(4)*0.54]);
hold(axInset, "on");
box(axInset, "on");

allRatio = [];
for s = 1:numel(sysOrder)
    idx = D.system == sysOrder(s);
    ratio = double(D.Lrest_over_Dmax_tau50(idx));
    ratio = ratio(isfinite(ratio) & ratio > 0);
    allRatio = [allRatio; ratio];
    scatter(axInset, 1 + 0.080*sin((1:numel(ratio))'*(3.1 + s)), ratio, ...
        12, sysColors(s,:), "filled", "MarkerFaceAlpha", 0.78, "MarkerEdgeColor", "none");
end


pooled = RatioSummary(string(RatioSummary.system) == "All systems", :);
if height(pooled) ~= 1
    error("Fig. 3c ratio_summary must contain exactly one All systems row; found %d.", height(pooled));
end

q = [ ...
    double(pooled.Q05(1)), ...
    double(pooled.Q25(1)), ...
    double(pooled.median(1)), ...
    double(pooled.Q75(1)), ...
    double(pooled.Q95(1))];

if any(~isfinite(q)) || any(q <= 0)
    error("Fig. 3c All systems ratio_summary contains invalid quantiles.");
end

if ~isempty(allRatio)
    set(axInset, "YScale", "log");
    drawMiniBoxFromSummary(axInset, 1, q, [0.45 0.45 0.45]);
    yline(axInset, 1, "--", "Color", [0.45 0.45 0.45], "LineWidth", 0.8);
    lo = min([allRatio; q(:)]);
    hi = max([allRatio; q(:)]);
    pad = 10 ^ (0.06 * max(1, log10(hi / lo)));
    ylim(axInset, [lo / pad, hi * pad]);
end
xlim(axInset, [0.55 1.45]);
setLogTicksForVisibleRange(axInset, "y");
set(axInset, "XTick", [], "XTickLabel", [], ...
    "FontName", "Arial", "FontSize", FONT_INSET, "LineWidth", 0.7, ...
    "TickDir", "out", "Layer", "top");
hY = ylabel(axInset, "$\kappa$", "FontName", "Arial", "FontSize", FONT_INSET, "Interpreter", "latex");
set(hY, "Units", "normalized");
ylPos = get(hY, "Position");
ylPos(1) = -0.16;
set(hY, "Position", ylPos);
end

function setLogTicksForVisibleRange(ax, axisName, maxTicks)
if nargin < 3
    maxTicks = Inf;
end
limits = xlim(ax);
if axisName == "y"
    limits = ylim(ax);
end
ticks = 10.^(ceil(log10(limits(1))):floor(log10(limits(2))));
ticks = ticks(ticks >= limits(1) & ticks <= limits(2));
if numel(ticks) < 2
    ticks = limits;
end
if isfinite(maxTicks) && numel(ticks) > maxTicks
    idx = unique(round(linspace(1, numel(ticks), maxTicks)));
    ticks = ticks(idx);
end
if axisName == "y"
    set(ax, "YTick", ticks);
else
    set(ax, "XTick", ticks);
end
end

function drawMiniBoxFromSummary(ax, x, q, color)
q = double(q(:)');
if numel(q) ~= 5 || any(~isfinite(q)) || any(q <= 0)
    return;
end
boxW = 0.22;
patch(ax, [x-boxW x+boxW x+boxW x-boxW], [q(2) q(2) q(4) q(4)], color, ...
    "FaceAlpha", 0.20, "EdgeColor", darkenColor(color, 0.35), "LineWidth", 0.85);
plot(ax, [x x], [q(1) q(2)], "-", "Color", darkenColor(color, 0.35), "LineWidth", 0.75);
plot(ax, [x x], [q(4) q(5)], "-", "Color", darkenColor(color, 0.35), "LineWidth", 0.75);
plot(ax, [x-boxW*0.55 x+boxW*0.55], [q(1) q(1)], "-", "Color", darkenColor(color, 0.35), "LineWidth", 0.75);
plot(ax, [x-boxW*0.55 x+boxW*0.55], [q(5) q(5)], "-", "Color", darkenColor(color, 0.35), "LineWidth", 0.75);
plot(ax, [x-boxW x+boxW], [q(3) q(3)], "k-", "LineWidth", 1.0);
plot(ax, x, q(3), "o", "MarkerFaceColor", "w", "MarkerEdgeColor", [0.05 0.05 0.05], "MarkerSize", 4.8, "LineWidth", 0.8);
end

function drawTailLossPanel(ax,D,Summary,sysOrder,sysNames,sysColors)
global FONT_AXIS_LABEL FONT_AXIS_TICK FONT_COUNT
hold(ax,"on"); tailN=zeros(numel(sysOrder)+1,1);
pooledColor=[0.42 0.42 0.42];


pooledVals=D.value(ismember(D.system,sysOrder));
pooledVals=pooledVals(isfinite(pooledVals)&pooledVals>=0);
pooledSummary=Summary(Summary.system=="All systems",:);
if height(pooledSummary)~=1
    error("Expected one Fig. 3d summary row for All systems; found %d.",height(pooledSummary));
end
qPooled=[double(pooledSummary.Q05),double(pooledSummary.Q25),double(pooledSummary.median),double(pooledSummary.Q75),double(pooledSummary.Q95)];
tailN(1)=double(pooledSummary.n);
if ~isempty(pooledVals)
    drawHalfViolin(ax,1,pooledVals,0.18,pooledColor,"right",false);
    for s = 1:numel(sysOrder)
        valsSys = D.value(D.system==sysOrder(s));
        valsSys = valsSys(isfinite(valsSys)&valsSys>=0);
        if isempty(valsSys), continue; end
        jitterSys = 0.055*sin((1:numel(valsSys))'*(2.7 + 0.21*s));
        scatter(ax,1+jitterSys,valsSys,17,sysColors(s,:),"filled", ...
            "MarkerFaceAlpha",0.82,"MarkerEdgeColor","none");
    end
    if all(isfinite(qPooled)), drawBoxDot(ax,1,qPooled,pooledColor,0.28); end
end


for s=1:numel(sysOrder)
    xPos=s+1;
    vals=D.value(D.system==sysOrder(s)); vals=vals(isfinite(vals)&vals>=0);
    idxSummary=Summary.system==sysOrder(s);
    if nnz(idxSummary)~=1, error("Expected one Fig. 3d summary row for %s; found %d.",sysOrder(s),nnz(idxSummary)); end
    row=Summary(idxSummary,:);
    q=[double(row.Q05),double(row.Q25),double(row.median),double(row.Q75),double(row.Q95)];
    tailN(s+1)=double(row.n);
    drawHalfViolin(ax,xPos,vals,0.18,sysColors(s,:),"right",false,darkenColor(sysColors(s,:),0.22));
    if all(isfinite(q)), drawBoxDot(ax,xPos,q,sysColors(s,:),0.28); end
end

if isempty(pooledVals)
    yUpper = 0.05;
else
    yMax = max(pooledVals);
    yUpper = max(0.05, ceil((1.02*yMax)/0.05)*0.05);
end


countX = linspace(0.065,0.86,4);
countStrings = {sprintf('$n=%d$',tailN(1)), ...
    sprintf('$n_E=%d$',tailN(2)), ...
    sprintf('$n_W=%d$',tailN(3)), ...
    sprintf('$n_G=%d$',tailN(4))};
for ii=1:4
    text(ax,countX(ii),1.018,countStrings{ii},"Units","normalized", ...
        "HorizontalAlignment","center","VerticalAlignment","bottom", ...
        "FontSize",FONT_COUNT,"Interpreter","latex","Clipping","off");
end

set(ax,"XTick",1:4,"XTickLabel",{"","","",""},"TickLabelInterpreter","none","FontSize",FONT_AXIS_TICK,"TickLength",[0 0]);
ylabel(ax,"\eta_{90} (%)","FontSize",FONT_AXIS_LABEL,"Interpreter","tex");
ylim(ax,[0 yUpper]); xlim(ax,[0.45 4.55]);
if yUpper <= 0.35
    yStep = 0.05;
elseif yUpper <= 0.65
    yStep = 0.10;
else
    yStep = 0.20;
end
yTicks = 0:yStep:yUpper;
yticks(ax,yTicks); set(ax,"YTickLabel",compose("%d",round(100*yTicks)));

ySystem = -0.01;
labels = {'All','E','W','G'};
for ii=1:4
    xNorm=(ii-0.45)/(4.55-0.45);
    text(ax,xNorm,ySystem,labels{ii},"Units","normalized", ...
        "HorizontalAlignment","center","VerticalAlignment","top", ...
        "FontSize",FONT_AXIS_TICK,"Interpreter","none","Clipping","off");
end
box(ax,"off");
end

function drawTimingRatioPanel(ax,Source,Summary,compactN)
global FONT_AXIS_LABEL FONT_AXIS_TICK FONT_COUNT

blue=[117 153 95]./255;
blueLight=[196 222 151]./255;
gold=[122 90 151]./255;
goldLight=[149 120 176]./255;

milestones=["T90","T95"];
metricLabels=["$p=90$","$p=95$"];

allVals=[];
qWAll=nan(2,5);
qGAll=nan(2,5);
nW=zeros(2,1);
nG=zeros(2,1);

hold(ax,"on");

for m=1:2

    valsW = timingSourceValues(Source,milestones(m),"E->W");
    valsG = timingSourceValues(Source,milestones(m),"E->G");

    allVals=[allVals;valsW;valsG];

    drawHalfViolin(ax,m,valsW,0.27,blueLight,"left",true,darkenColor(blue,0.22));
    drawHalfViolin(ax,m,valsG,0.27,goldLight,"right",true,darkenColor(gold,0.22));

    [qW,nW(m)] = timingSummaryValues(Summary,milestones(m),"E->W");
    [qG,nG(m)] = timingSummaryValues(Summary,milestones(m),"E->G");

    qWAll(m,:) = qW;
    qGAll(m,:) = qG;

    if all(isfinite(qW))
        drawBoxDot(ax,m-0.07,qW,blue,0.055);
    end
    if all(isfinite(qG))
        drawBoxDot(ax,m+0.07,qG,gold,0.055);
    end
end

set(ax, ...
    "YScale","log", ...
    "XTick",1:2, ...
    "XTickLabel",metricLabels, ...
    "TickLabelInterpreter","latex", ...
    "TickLength",[0 0], ...
    "FontSize",FONT_AXIS_TICK);

yline(ax,1,"--","Color",[0.25 0.25 0.25],"LineWidth",1.0);
xlim(ax,[0.48 2.62]);

if isempty(allVals)
    ylim(ax,[0.04 180]);
else
    ylim(ax,[max(0.01,min(allVals)/1.8), min(250,max(180,max(allVals)*1.6))]);
end

ylabel(ax,"$\rho_p^{B/E}$", ...
    "Interpreter","latex", ...
    "FontSize",FONT_AXIS_LABEL);

h1 = plot(ax,NaN,NaN,"o", ...
    "MarkerFaceColor",blueLight, ...
    "MarkerEdgeColor","none", ...
    "MarkerSize",8.5, ...
    "LineStyle","none");

h2 = plot(ax,NaN,NaN,"o", ...
    "MarkerFaceColor",goldLight, ...
    "MarkerEdgeColor","none", ...
    "MarkerSize",8.5, ...
    "LineStyle","none");

lgd = legend(ax,[h1 h2],["W/E","G/E"],"Location","southeast","Box","off");
styleLegend(lgd,FONT_COUNT);

if compactN
    drawAlignedNEntryNormalized(ax,0.255,1.075,'E,W',nW(1),FONT_COUNT);
    drawAlignedNEntryNormalized(ax,0.255,1.015,'E,G',nG(1),FONT_COUNT);
    drawAlignedNEntryNormalized(ax,0.735,1.075,'E,W',nW(2),FONT_COUNT);
    drawAlignedNEntryNormalized(ax,0.735,1.015,'E,G',nG(2),FONT_COUNT);
else
    text(ax,0.275,1.055, ...
        sprintf('$n_{E,W}=%d; n_{E,G}=%d$',nW(1),nG(1)), ...
        "Units","normalized", ...
        "HorizontalAlignment","center", ...
        "VerticalAlignment","bottom", ...
        "FontSize",FONT_COUNT, ...
        "Interpreter","latex", ...
        "Clipping","off");

    text(ax,0.725,1.055, ...
        sprintf('$n_{E,W}=%d; n_{E,G}=%d$',nW(2),nG(2)), ...
        "Units","normalized", ...
        "HorizontalAlignment","center", ...
        "VerticalAlignment","bottom", ...
        "FontSize",FONT_COUNT, ...
        "Interpreter","latex", ...
        "Clipping","off");
end



medianLabelColor = [0 0 0];
addMedianOnly(ax,0.70,qWAll(1,3) * 1.62,medianLabelColor,"center",qWAll(1,3),-1);
addMedianOnly(ax,1.44,qGAll(1,3) * 1.25,medianLabelColor,"center",qGAll(1,3));
addMedianOnly(ax,1.60,qWAll(2,3),medianLabelColor,"center");
addMedianOnly(ax,2.46,qGAll(2,3),medianLabelColor,"center");

box(ax,"off");
end

function values=timingSourceValues(Source,milestone,pair)
idx=Source.milestone==milestone & Source.pair==pair; values=double(Source.ratio(idx)); values=values(isfinite(values)&values>0);
end
function [q,n]=timingSummaryValues(Summary,milestone,pair)
idx=Summary.milestone==milestone & Summary.pair==pair;
if nnz(idx)~=1, error("Expected one Fig. 3e summary row for %s / %s; found %d.",milestone,pair,nnz(idx)); end
row=Summary(idx,:); q=[double(row.Q05),double(row.Q25),double(row.median),double(row.Q75),double(row.Q95)]; n=double(row.n);
end

function drawAlignedNEntryNormalized(ax,xEq,y,label,nValue,fontSize)

prefix = sprintf('$n_{%s}$', char(string(label)));
valueText = sprintf('$%d$', round(nValue));
gapLeft = 0.038;
gapRight = 0.038;
if string(label) == "E,G"
    gapLeft = 0.052;
end
text(ax,xEq-gapLeft,y,prefix, ...
    "Units","normalized", ...
    "HorizontalAlignment","right","VerticalAlignment","bottom", ...
    "FontSize",fontSize,"Interpreter","latex","Clipping","off");
text(ax,xEq,y,'$=$', ...
    "Units","normalized", ...
    "HorizontalAlignment","center","VerticalAlignment","bottom", ...
    "FontSize",fontSize,"Interpreter","latex","Clipping","off");
text(ax,xEq+gapRight,y,valueText, ...
    "Units","normalized", ...
    "HorizontalAlignment","left","VerticalAlignment","bottom", ...
    "FontSize",fontSize,"Interpreter","latex","Clipping","off");
end

function drawCmaxDistributionPanel(ax,Source,Summary)
global FONT_AXIS_LABEL FONT_AXIS_TICK FONT_COUNT FONT_ANNOT
pairs=["EW","EG"]; colors=[117 153 95;122 90 151]./255; hold(ax,"on");
legendHandles=gobjects(2,1); legendLabels=strings(2,1); fractions=nan(2,1);
for i=1:2
    idx=Source.system_pair==pairs(i)&isfinite(Source.CmaxAB90); values=double(Source.CmaxAB90(idx));
    if any(values < -1e-12)
        error("Fig. 3f CmaxAB90 contains negative value(s) for %s.",pairs(i));
    end
    values=sort(values);
    if isempty(values), continue; end
    sidx=Summary.system_pair==pairs(i); if nnz(sidx)~=1, error("Expected one Fig. 3f summary row for %s; found %d.",pairs(i),nnz(sidx)); end
    row=Summary(sidx,:); n=double(row.n); fractions(i)=double(row.fraction_le_5pct);
    y=100.*(1:numel(values))'./numel(values);
    if pairs(i)=="EW"
        legendLabels(i)="$n_{E,W}="+string(n)+"$";
    else
        legendLabels(i)="$n_{E,G}="+string(n)+"$";
    end
    legendHandles(i)=stairs(ax,100.*values,y,"LineWidth",2.2,"Color",colors(i,:),"DisplayName",legendLabels(i));
    scatter(ax,100.*values,y,14,colors(i,:),"filled","MarkerFaceAlpha",0.35,"MarkerEdgeColor","none","HandleVisibility","off");
end
xline(ax,5,"--","Color",[0.20 0.20 0.20],"LineWidth",1.2,"HandleVisibility","off");
xlabF = xlabel(ax,"$C_{E\rightarrow B}^{\max}(90)\;(\%)$","FontSize",FONT_AXIS_LABEL,"Interpreter","latex");
set(xlabF,"Units","normalized");
posF = get(xlabF,"Position"); posF(2) = -0.08; set(xlabF,"Position",posF);
ylabel(ax,"ECDF (%)","FontSize",FONT_AXIS_LABEL);
set(ax,"FontSize",FONT_AXIS_TICK,"LineWidth",0.9,"Box","off","Layer","top"); xlim(ax,[0 100]); ylim(ax,[0 102]);
xticks(ax,[5 20 40 60 80 100]);
valid=isgraphics(legendHandles); lgd=legend(ax,legendHandles(valid),legendLabels(valid),"Location","southeast","Box","off"); lgd.Interpreter="latex"; styleLegend(lgd,FONT_COUNT);
for i=1:numel(pairs)
    if ~isfinite(fractions(i)), continue; end
    yText=100*fractions(i);
    if i==1, yText=yText-2.0; else, yText=yText-4.0; end
    text(ax,7.5,yText,sprintf("%.0f%%",100*fractions(i)), ...
        "HorizontalAlignment","left","VerticalAlignment","middle", ...
        "FontSize",FONT_COUNT,"Color",darkenColor(colors(i,:),0.12));
end
grid(ax,"off");
end

function lbl = buildCmaxLegendLabel(pair, n)
pair = upper(strtrim(string(pair)));
switch pair
    case "EW"
        lbl = sprintf("E\rightarrowW (n_{E,W} = %d)", n);
    case "EG"
        lbl = sprintf("E\rightarrowG (n_{E,G} = %d)", n);
    otherwise
        lbl = sprintf("%s (n = %d)", char(pair), n);
end
end

function tf = toLogicalColumn(v)
if islogical(v)
    tf = v;
    return;
end
if isnumeric(v)
    tf = isfinite(v) & v ~= 0;
    return;
end
vs = lower(strtrim(string(v)));
tf = vs == "true" | vs == "1" | vs == "yes";
end

function drawHalfViolin(ax, x, vals, width, color, side, useLog, scatterColor)
vals = vals(:);
vals = vals(isfinite(vals));
if nargin < 8 || isempty(scatterColor)
    scatterColor = color;
end
if useLog
    vals = vals(vals > 0);
    z = log10(vals);
else
    z = vals;
end
if isempty(vals), return; end
if numel(vals) < 3 || range(z) <= 1e-12
    y = median(vals, "omitnan");
    scatter(ax, x, y, 32, scatterColor, "filled", "MarkerFaceAlpha", 0.85, ...
        "MarkerEdgeColor", "none", "HandleVisibility", "off");
    return;
end
try
    [f, zi] = ksdensity(z, "NumPoints", 120);
catch
    zi = linspace(min(z), max(z), 120);
    edges = linspace(min(z), max(z), 22);
    counts = histcounts(z, edges, "Normalization", "pdf");
    mids = (edges(1:end-1) + edges(2:end)) ./ 2;
    f = interp1(mids, counts, zi, "linear", 0);
end
if max(f) <= 0, return; end
f = f ./ max(f) .* width;
if useLog, yi = 10.^zi; else, yi = zi; end
if side == "left"
    xx = [x; x - f(:); x];
else
    xx = [x; x + f(:); x];
end
yy = [yi(1); yi(:); yi(end)];
patch(ax, xx, yy, color, "FaceAlpha", 0.24, "EdgeColor", "none", "HandleVisibility", "off");
jit = 0.028 * sin((1:numel(vals))' * 2.399963);
if side == "left", jit = -abs(jit); else, jit = abs(jit); end
scatter(ax, x + jit, vals, 16, scatterColor, "filled", "MarkerFaceAlpha", 0.82, ...
    "MarkerEdgeColor", "none", "HandleVisibility", "off");
end

function drawBoxDot(ax, x, q, color, halfWidth)
if any(~isfinite(q)), return; end
patch(ax, [x-halfWidth x+halfWidth x+halfWidth x-halfWidth], [q(2) q(2) q(4) q(4)], ...
    color, "FaceAlpha", 0.26, "EdgeColor", darkenColor(color, 0.35), "LineWidth", 0.8, "HandleVisibility", "off");
plot(ax, [x x], [q(1) q(5)], "-", "Color", darkenColor(color,0.35), "LineWidth", 1.0, "HandleVisibility", "off");
plot(ax, [x-halfWidth*0.55 x+halfWidth*0.55], [q(1) q(1)], "-", ...
    "Color", darkenColor(color,0.35), "LineWidth", 0.8, "HandleVisibility", "off");
plot(ax, [x-halfWidth*0.55 x+halfWidth*0.55], [q(5) q(5)], "-", ...
    "Color", darkenColor(color,0.35), "LineWidth", 0.8, "HandleVisibility", "off");
plot(ax, [x-halfWidth x+halfWidth], [q(3) q(3)], "w-", "LineWidth", 2.0, "HandleVisibility", "off");
plot(ax, x, q(3), "o", "MarkerFaceColor", "w", ...
    "MarkerEdgeColor", darkenColor(color,0.20), "MarkerSize", 5.8, ...
    "LineWidth", 0.9, "HandleVisibility", "off");
end

function addMedianOnly(ax, x, y, color, align, displayValue, fontDelta)
global FONT_ANNOT
if ~isfinite(y), return; end
if nargin < 6 || isempty(displayValue)
    displayValue = y;
end
if nargin < 7 || isempty(fontDelta)
    fontDelta = 0;
end
text(ax, x, y, sprintf("%.3f", displayValue), ...
    "Color", [0.22 0.22 0.22], ...
    "FontWeight", "normal", ...
    "FontName", "Arial", ...
    "Interpreter", "none", ...
    "HorizontalAlignment", align, ...
    "VerticalAlignment", "middle", ...
    "FontSize", FONT_ANNOT - 2 + fontDelta, ...
    "Tag", "TimingMedianLabel");
end

function c = sysColorFromCode(code)
code = upper(strtrim(string(code)));
switch code
    case "E"
        c = [0.70 0.14 0.16];
    case "W"
        c = [0.17 0.46 0.72];
    case "G"
        c = [0.78 0.55 0.18];
    otherwise
        c = [0.45 0.45 0.45];
end
end

function c = darkenColor(c, f)
c = max(0, c .* (1 - f));
end

function addCustomNColumns(ax, labels, values, x0, y0, dy, fontSize)
labels = string(labels(:));
values = values(:);
xLabel = x0; xEq = x0 + 0.095; xVal = x0 + 0.128;
for ii = 1:numel(values)
    y = y0 - (ii-1) * dy;
    text(ax, xLabel, y, "$" + labels(ii) + "$", "Units", "normalized", "HorizontalAlignment", "left", ...
        "VerticalAlignment", "top", "FontSize", fontSize, "FontWeight", "normal", ...
        "Interpreter", "latex", "Color", [0.10 0.10 0.10], "Margin", 0.5);
    text(ax, xEq, y, "=", "Units", "normalized", "HorizontalAlignment", "center", ...
        "VerticalAlignment", "top", "FontSize", fontSize, "FontWeight", "normal", ...
        "Interpreter", "none", "Color", [0.10 0.10 0.10], "Margin", 0.5);
    text(ax, xVal, y, sprintf("%d", values(ii)), "Units", "normalized", "HorizontalAlignment", "left", ...
        "VerticalAlignment", "top", "FontSize", fontSize, "FontWeight", "normal", ...
        "Interpreter", "none", "Color", [0.10 0.10 0.10], "Margin", 0.5);
end

function name = sysOrderDisplayName(code)
code = upper(strtrim(string(code)));
switch code
    case "E"
        name = "Electric power";
    case "W"
        name = "Water supply";
    case "G"
        name = "Natural gas";
    otherwise
        name = char(code);
end
end
end

function xlimLogData(ax, x)
x = x(isfinite(x) & x > 0);
if isempty(x)
    xlim(ax, [1e-4 1]);
    return;
end
xmin = min(x); xmax = max(x);
pad = 10^(0.06 * max(1, log10(xmax / xmin)));
xlim(ax, [xmin / pad, xmax * pad]);
end

function ylimLogData(ax, y)
y = y(isfinite(y) & y > 0);
if isempty(y)
    ylim(ax, [1e-4 1]);
    return;
end
ymin = min(y); ymax = max(y);
pad = 10^(0.06 * max(1, log10(ymax / ymin)));
ylim(ax, [ymin / pad, ymax * pad]);
end

function alignYLabels(axList)
if isempty(axList)
    return;
end
labelList = gobjects(0);
axPosList = [];
figXPos = [];
for i = 1:numel(axList)
    if ~isgraphics(axList(i)), continue; end
    hY = get(axList(i), "YLabel");
    if ~isgraphics(hY), continue; end
    set(hY, "Units", "normalized");
    pos = get(hY, "Position");
    axPos = get(axList(i), "Position");
    labelList(end+1) = hY;
    axPosList(end+1,:) = axPos;
    figXPos(end+1) = axPos(1) + pos(1) * axPos(3);
end
if isempty(labelList)
    return;
end
targetFigX = min(figXPos);
for i = 1:numel(labelList)
    pos = get(labelList(i), "Position");
    axPos = axPosList(i,:);
    pos(1) = (targetFigX - axPos(1)) / axPos(3);
    set(labelList(i), "Position", pos);
end
end

function unifyXAxisFont(axList, fontSize)
for i = 1:numel(axList)
    if ~isgraphics(axList(i)), continue; end
    try
        axList(i).XAxis.FontSize = fontSize;
    catch
        set(axList(i), "FontSize", fontSize);
    end
end
end

function addPanelLetter(fig, ax, letter, xAnchorAx, yExtra)
global FONT_PANEL
pos = get(ax, "Position");
letterDx = 0.047;
if nargin < 5 || isempty(yExtra)
    yExtra = 0;
end
if nargin >= 4 && ~isempty(xAnchorAx) && isgraphics(xAnchorAx)
    anchorPos = get(xAnchorAx, "Position");
    x = anchorPos(1) - letterDx;
else
    x = pos(1) - letterDx;
end
y = pos(2) + pos(4) + 0.010 + yExtra;
annotation(fig, "textbox", [x y 0.03 0.03], "String", letter, ...
    "EdgeColor", "none", "HorizontalAlignment", "left", ...
    "VerticalAlignment", "bottom", "FontWeight", "bold", ...
    "FontSize", FONT_PANEL, "FontName", "Arial", "Color", [0.05 0.05 0.05], ...
    "Tag", "FigurePanelLetter");
end

function styleLegend(lgd, fontSize)
if isempty(lgd) || ~isgraphics(lgd), return; end
lgd.Box = "off";
lgd.FontSize = fontSize;
lgd.ItemTokenSize = [16 12];
end

function applyUniformFigureTypography(fig, regularFontSize, axisLabelFontSize, panelFontSize)

objs = findall(fig, "-property", "FontSize");
for i = 1:numel(objs)
    isPanelLetter = false;
    isTimingMedian = false;
    try
        if isprop(objs(i), "Tag")
            isPanelLetter = strcmp(string(objs(i).Tag), "FigurePanelLetter");
            isTimingMedian = strcmp(string(objs(i).Tag), "TimingMedianLabel");
        end
    catch
    end

    try
        if isPanelLetter
            objs(i).FontSize = panelFontSize;
        elseif isTimingMedian
            objs(i).FontSize = regularFontSize - 1.5;
        else
            objs(i).FontSize = regularFontSize;
        end
    catch
    end

    try
        if isprop(objs(i), "FontWeight")
            if isPanelLetter
                objs(i).FontWeight = "bold";
            else
                objs(i).FontWeight = "normal";
            end
        end
    catch
    end
end

axList = findall(fig, "Type", "axes");
for i = 1:numel(axList)

    try
        axList(i).FontSize = regularFontSize;
    catch
    end

    labelHandles = [
        axList(i).XLabel
        axList(i).YLabel
        axList(i).ZLabel];

    for j = 1:numel(labelHandles)
        if isgraphics(labelHandles(j))
            try
                labelHandles(j).FontSize = axisLabelFontSize;
                labelHandles(j).FontWeight = "normal";
            catch
            end
        end
    end
end
end

function setFigureDefaults(fig)
global FONT_AXIS_TICK
set(fig, "WindowStyle", "normal", "DefaultAxesFontName", "Arial", ...
    "DefaultTextFontName", "Arial", "DefaultAxesFontSize", FONT_AXIS_TICK, "DefaultTextFontSize", FONT_AXIS_TICK);
end

function hideAxesToolbars(fig)
axList = findall(fig, "Type", "axes");
for i = 1:numel(axList)
    try
        if isprop(axList(i), "Toolbar") && ~isempty(axList(i).Toolbar)
            axList(i).Toolbar.Visible = "off";
        end
    catch
    end
end
end

function turnAllGridsOff(fig)
axList = findall(fig, "Type", "axes");
for i = 1:numel(axList)
    try, grid(axList(i), "off"); catch, end
end
end

function exportFixedCanvasPng(fig, filename, figWidth, figHeight)
dpi = 100;
set(fig, "PaperUnits", "inches", "PaperPosition", [0 0 figWidth/dpi figHeight/dpi], ...
    "PaperSize", [figWidth/dpi figHeight/dpi], "PaperPositionMode", "manual", "InvertHardcopy", "off");
print(fig, char(filename), "-dpng", "-r" + string(dpi));
end

function prepareFigureForVisibleFigSave(fig, figWidth, figHeight)
screen = get(groot, "ScreenSize");
maxW = 0.70 * screen(3);
maxH = 0.78 * screen(4);
displayW = min(figWidth, maxW);
displayH = displayW * figHeight / figWidth;
if displayH > maxH
    displayH = maxH;
    displayW = displayH * figWidth / figHeight;
end
figScale = displayW / figWidth;
scaleFigureGraphicProperties(fig, figScale);
setVisibleFigAxisLabelFontSize(fig, 27.0 * figScale);
setTaggedTextFontSize(fig, "PanelBRatioLabel", 27.0 * figScale);
left = max(40, screen(1) + 0.04 * screen(3));
bottom = max(40, screen(2) + 0.08 * screen(4));
set(fig, "Visible", "on", "WindowState", "normal", "Units", "pixels", ...
    "Position", [left bottom displayW displayH]);
drawnow;
end

function setVisibleFigAxisLabelFontSize(fig, labelFontSize)
axList = findall(fig, "Type", "axes");
for i = 1:numel(axList)
    labelHandles = [axList(i).XLabel; axList(i).YLabel; axList(i).ZLabel];
    for j = 1:numel(labelHandles)
        if isgraphics(labelHandles(j))
            try
                labelHandles(j).FontSize = labelFontSize;
            catch
            end
        end
    end
end
end

function setTaggedTextFontSize(fig, tagName, fontSize)
textList = findall(fig, "Type", "text", "Tag", tagName);
for i = 1:numel(textList)
    try
        textList(i).FontSize = fontSize;
    catch
    end
end
end

function scaleFigureGraphicProperties(fig, scale)
if ~isfinite(scale) || scale <= 0 || abs(scale - 1) < 1e-6
    return;
end
objects = findall(fig);
for k = 1:numel(objects)
    h = objects(k);
    if isprop(h, "FontSize")
        try
            h.FontSize = max(5, h.FontSize * scale);
        catch
        end
    end
    if isprop(h, "LineWidth")
        try
            h.LineWidth = max(0.35, h.LineWidth * scale);
        catch
        end
    end
    if isprop(h, "MarkerSize")
        try
            h.MarkerSize = max(2, h.MarkerSize * scale);
        catch
        end
    end
    if isprop(h, "SizeData")
        try
            h.SizeData = max(4, h.SizeData * scale.^2);
        catch
        end
    end
end
end

function ensureDir(d)
if ~exist(d, "dir")
    mkdir(d);
end
end

function [dataSource, outputDir, dataRoot] = parseFigureOptions(moduleRoot, varargin)
% Preserve the historical one-argument output-directory call.
if numel(varargin) == 1 && isempty(varargin{1})
    varargin = {};
elseif numel(varargin) == 1 && isTextScalar(varargin{1}) && ...
        ~any(strcmpi(char(varargin{1}), {'DataSource', 'OutputDir', 'DataRoot'}))
    varargin = {"OutputDir", varargin{1}};
end

p = inputParser;
p.addParameter("DataSource", "packaged", @(x) isTextScalar(x));
p.addParameter("OutputDir", fullfile(moduleRoot, "outputs", "figures", "main"), ...
    @(x) isTextScalar(x));
p.addParameter("DataRoot", "", @(x) isTextScalar(x));
p.parse(varargin{:});

dataSource = lower(char(string(p.Results.DataSource)));
if ~ismember(dataSource, {'packaged', 'generated'})
    error("EmergentFigure:InvalidDataSource", ...
        "DataSource must be 'packaged' or 'generated'.");
end
outputDir = string(p.Results.OutputDir);
if strlength(outputDir) == 0
    outputDir = fullfile(moduleRoot, "outputs", "figures", "main");
end
dataRoot = string(p.Results.DataRoot);
if strlength(dataRoot) > 0 && ~isfile(fullfile(dataRoot, "SourceData_Figure3.xlsx"))
    error("EmergentFigure:MissingData", ...
        "Figure 3 source workbook was not found under DataRoot: %s", dataRoot);
end
end

function tf = isTextScalar(value)
tf = (ischar(value) && (isrow(value) || isempty(value))) || ...
    (isstring(value) && isscalar(value));
end
