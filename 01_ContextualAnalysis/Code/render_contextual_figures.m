function render_contextual_figures(projectRoot, runDir, dataDir, figDir) %#ok<INUSD>
% Render Figure 2 from CSV inputs, without model fits or historical FIG files.
if nargin < 2 || isempty(runDir) || strlength(string(runDir)) == 0
    runDir = fullfile(ctx_module_root(), "outputs", "analysis");
end
if nargin < 3 || isempty(dataDir) || strlength(string(dataDir)) == 0
    dataDir = fullfile(runDir, "figure_data");
end
if nargin < 4 || isempty(figDir) || strlength(string(figDir)) == 0
    figDir = fullfile(runDir, "figures");
end
% projectRoot is retained for compatibility with existing callers.
dataDir = string(dataDir);
figDir = string(figDir);
required = ["figure2a_D0_records.csv", "figure2b_T90_records.csv", ...
    "figure2c_fitted_r2.csv", "figure2d_adjusted_T90_factors.csv"];
if any(~isfile(fullfile(dataDir, required)))
    error("ContextualFigure:MissingData", "Figure 2 data files are missing from %s.", dataDir);
end
A = readtable(fullfile(dataDir, required(1)), "TextType", "string");
B = readtable(fullfile(dataDir, required(2)), "TextType", "string");
C = readtable(fullfile(dataDir, required(3)), "TextType", "string");
D = readtable(fullfile(dataDir, required(4)), "TextType", "string");
if ~isfolder(figDir), mkdir(figDir); end

fig = figure("Name", "Figure 2 - systems and contexts", "NumberTitle", "off", ...
    "Color", "w", "Units", "centimeters", "Position", [1 1 18.31 17.75], ...
    "Visible", "on", "Tag", "ContextualFigure2", ...
    "PaperUnits", "centimeters", "PaperSize", [26 28], ...
    "PaperPosition", [0 0 26 28], "PaperPositionMode", "manual");
% Window creation can clamp the requested height to the current display.
drawnow;
fig.Position = [1 1 18.31 17.75];
fig.UserData = struct("dataDirectory", dataDir, "sourceFiles", required);
colors = [.75 .22 .25; .20 .50 .78; .86 .55 .16];
axA = axes(fig, "Position", [.090 .725 .370 .225], "Tag", "PanelA");
drawD0(axA, A, colors);
axB = axes(fig, "Position", [.090 .422 .370 .235], "Tag", "PanelB");
drawT90(axB, B, colors);
axC = axes(fig, "Position", [.090 .085 .370 .230], "Tag", "PanelC");
drawR2(axC, C);
axD = axes(fig, "Position", [.640 .105 .244 .835], "Tag", "PanelD");
drawFactors(axD, D);
panelLetter(fig, [.035 .953 .025 .025], "a");
panelLetter(fig, [.035 .664 .025 .025], "b");
panelLetter(fig, [.035 .318 .025 .025], "c");
panelLetter(fig, [.505 .953 .025 .025], "d");

drawnow;
fig.Position = [1 1 18.31 17.75];
drawnow;
name = "Figure_2_contextual_analysis";
savefig(fig, fullfile(figDir, name + ".fig"));
exportgraphics(fig, fullfile(figDir, name + ".png"), "Resolution", 300);
reopened = openfig(fullfile(figDir, name + ".fig"), "visible");
drawnow;
exportgraphics(reopened, fullfile(figDir, name + "_reopened.png"), "Resolution", 180);
close(reopened);
fprintf("Figure 2 rendered from CSV data: %s\n", fullfile(figDir, name));
end

function drawD0(ax, T, colors)
systems = ["E" "W" "G"];
names = ["Electric power" "Water supply" "Natural gas"];
counts = zeros(3, 3);
for s = 1:3
    values = T.outage0(T.sys == systems(s));
    counts(:, s) = [sum(values == 0); sum(values > 0 & values < 1); sum(values == 1)];
end
hold(ax, "on");
x = [1 2.25 3.5];
bars = bar(ax, x, counts, "stacked", "BarWidth", .46);
for s = 1:3
    bars(s).FaceColor = colors(s, :);
    bars(s).EdgeColor = "none";
    bars(s).Tag = "D0Counts_" + systems(s);
end
totals = sum(counts, 2);
for k = 1:3
    text(ax, x(k), totals(k) + 16, sprintf("%d", totals(k)), ...
        "HorizontalAlignment", "center", "FontName", "Arial", "FontSize", 7.5, ...
        "Color", [.18 .18 .18], "Tag", "D0Total_" + k);
end
ax.XLim = [.52 3.98];
ax.YLim = [0 max(820, max(totals) * 1.10)];
ax.YTick = 0:200:800;
ax.XTick = x;
ax.XTickLabel = ["\itD_{0}\rm = 0", "0 < \itD_{0}\rm < 1", "\itD_{0}\rm = 1"];
ylabel(ax, "Number of records");
styleAxes(ax, 8);
ax.UserData = struct("counts", counts, "systems", systems);
lgd = legend(ax, bars, names, "Orientation", "vertical", "NumColumns", 1, ...
    "Box", "off", "FontName", "Arial", "FontSize", 7.2, "AutoUpdate", "off");
lgd.Units = "normalized";
lgd.Position = [.307481 .859854 .162717 .057471];
lgd.Tag = "LifelineLegend";

inset = axes(ancestor(ax, "figure"), "Position", [.13748191 .825 .09351809 .082], ...
    "Color", "w", "Tag", "D0ECDF");
hold(inset, "on");
for s = 1:3
    values = sort(T.outage0(T.sys == systems(s) & T.outage0 > 0 & T.outage0 < 1));
    plot(inset, values, (1:numel(values))' / numel(values), ...
        "Color", colors(s, :), "LineWidth", 1.2, "Tag", "ECDF_" + systems(s));
end
inset.XLim = [0 1];
inset.YLim = [0 1];
inset.XTick = [0 .5 1];
inset.YTick = [0 .5 1];
xlabel(inset, "\itD_{0}\rm (0 < \itD_{0}\rm < 1)");
ylabel(inset, "ECDF");
styleAxes(inset, 6.5);
end

function drawT90(ax, T, colors)
systems = ["E" "W" "G"];
hold(ax, "on");
oldRng = rng;
restoreRng = onCleanup(@() rng(oldRng));
rng(19, "twister");
axisColor = [.18 .18 .18];
summary = zeros(3, 4);
annotationY = max(T.T90) * 1.72;
for s = 1:3
    values = T.T90(T.sys == systems(s));
    z = log(values);
    support = linspace(min(z), max(z), 256);
    density = ksdensity(z, support);
    width = .16 * density / max(density);
    patch(ax, [s + .07 + width, s + .07 + zeros(size(width))], ...
        [exp(support), fliplr(exp(support))], colors(s, :), ...
        "FaceAlpha", .22, "EdgeColor", "none", "Tag", "T90Violin_" + systems(s));
    scatter(ax, s + (rand(numel(values), 1) - .5) * .11, values, 8, ...
        colors(s, :), "filled", "MarkerFaceAlpha", .24, "MarkerEdgeColor", "none", ...
        "Tag", "T90Records_" + systems(s));
    q = prctile(values, [25 50 75]);
    spread = q(3) - q(1);
    low = min(values(values >= q(1) - 1.5 * spread));
    high = max(values(values <= q(3) + 1.5 * spread));
    line(ax, [s s], [low q(1)], "Color", axisColor, "LineWidth", .75);
    line(ax, [s s], [q(3) high], "Color", axisColor, "LineWidth", .75);
    line(ax, s + [-.035 .035], [low low], "Color", axisColor, "LineWidth", .75);
    line(ax, s + [-.035 .035], [high high], "Color", axisColor, "LineWidth", .75);
    patch(ax, s + [-.055 .055 .055 -.055], [q(1) q(1) q(3) q(3)], ...
        "w", "FaceColor", "none", "EdgeColor", axisColor, "LineWidth", .85, ...
        "Tag", "T90Box_" + systems(s));
    line(ax, s + [-.055 .055], [q(2) q(2)], "Color", axisColor, "LineWidth", 1.25);
    text(ax, s, annotationY, ...
        sprintf("\\itn_{%s}\\rm = %d\nmedian = %.3f\nIQR = %.3f-%.3f", ...
        systems(s), numel(values), q(2), q(1), q(3)), ...
        "Interpreter", "tex", "HorizontalAlignment", "center", ...
        "VerticalAlignment", "bottom", "FontName", "Arial", "FontSize", 7, ...
        "Color", axisColor, "Tag", "T90Summary_" + systems(s));
    summary(s, :) = [numel(values), q];
end
ax.YScale = "log";
ax.XLim = [.55 3.5];
ax.YLim = [min(T.T90) / 1.5, max(T.T90) * 8];
ax.YTick = [0.01 1 100];
ax.XTick = 1:3;
ax.XTickLabel = ["Electric power" "Water supply" "Natural gas"];
ylabel(ax, "\itT_{90}\rm (days)");
styleAxes(ax, 8);
ax.YMinorTick = "on";
ax.UserData = struct("summaryColumns", ["n" "q25" "median" "q75"], ...
    "summary", summary, "systems", systems);
end

function drawR2(ax, T)
models = ["system_and_context_full_D0_record_set"; ...
    "system_and_context_full_T90_record_set"; ...
    "system_and_context_D0_aligned_T90_record_set"; ...
    "system_and_context_plus_D0_D0_aligned_T90_record_set"];
T = orderedRows(T, "model", models);
hold(ax, "on");
axisColor = [.18 .18 .18];
contextColor = [.70 .75 .76];
plusColor = [.29 .43 .46];
x = [1 2.25 3.55 4.20];
colors = [contextColor; contextColor; contextColor; plusColor];
patch(ax, [.58 1.42 1.42 .58], [0 0 .78 .78], [.96 .94 .88], ...
    "EdgeColor", "none", "FaceAlpha", .75, "HandleVisibility", "off");
patch(ax, [1.72 4.55 4.55 1.72], [0 0 .78 .78], [.92 .96 .96], ...
    "EdgeColor", "none", "FaceAlpha", .75, "HandleVisibility", "off");
bars = gobjects(4, 1);
for k = 1:4
    bars(k) = bar(ax, x(k), T.fitted_r2(k), .43, "FaceColor", colors(k, :), ...
        "EdgeColor", axisColor, "LineWidth", .55, "Tag", "R2Bar_" + models(k));
    errorbar(ax, x(k), T.fitted_r2(k), T.fitted_r2(k) - T.r2_ci95_low(k), ...
        T.r2_ci95_high(k) - T.fitted_r2(k), "o", "Color", axisColor, ...
        "MarkerFaceColor", colors(k, :), "MarkerEdgeColor", axisColor, ...
        "MarkerSize", 4.8, "LineWidth", 1.1, "CapSize", 6, ...
        "Tag", "R2CI_" + models(k));
    if k <= 2
        labelX = x(k) + .18;
        align = "left";
    else
        labelX = x(k) - .18;
        align = "right";
    end
    text(ax, labelX, T.fitted_r2(k) + .028, sprintf("%.3f", T.fitted_r2(k)), ...
        "FontName", "Arial", "FontSize", 7, "Color", axisColor, ...
        "HorizontalAlignment", align, "VerticalAlignment", "middle");
end
bracketY = .660;
plot(ax, [x(3) x(3) x(4) x(4)], bracketY + [-.018 0 0 -.018], ...
    "Color", axisColor, "LineWidth", .8, "Tag", "DeltaR2Bracket");
text(ax, mean(x(3:4)), bracketY + .028, ...
    sprintf("\\Delta\\itR\\rm^{2} = %.3f\n(%.3f-%.3f)", ...
    T.delta_fitted_r2(4), T.delta_ci95_low(4), T.delta_ci95_high(4)), ...
    "Interpreter", "tex", "FontName", "Arial", "FontSize", 7, ...
    "HorizontalAlignment", "center", "VerticalAlignment", "bottom", ...
    "Color", axisColor, "Tag", "DeltaR2Summary");
text(ax, .64, .755, "\itD_{0}\rm", "FontName", "Arial", "FontSize", 7.5);
text(ax, 1.80, .755, "\itT_{90}\rm", "FontName", "Arial", "FontSize", 7.5);
ax.XLim = [.55 4.60];
ax.YLim = [0 .80];
ax.YTick = 0:.2:.8;
ax.XTick = x;
ax.XTickLabel = ["\itD_{0}\rm records"; "Full \itT_{90}\rm record set"; ...
    "\itD_{0}\rm-aligned"; "\itD_{0}\rm-aligned"];
ax.XTickLabelRotation = 20;
ylabel(ax, "Fitted \itR\rm^{2}");
styleAxes(ax, 7.2);
ax.UserData = T;
lgd = legend(ax, bars([1 4]), ...
    ["System-and-context model", "System-and-context-plus-\itD_{0}\rm model"], ...
    "Interpreter", "tex", "Orientation", "vertical", "NumColumns", 1, ...
    "Box", "off", "FontName", "Arial", "FontSize", 7, "AutoUpdate", "off");
lgd.Units = "normalized";
lgd.Position = [.0916774 .321052 .291618 .0412419];
lgd.Tag = "ModelLegend";
end

function drawFactors(ax, T)
order = ["water_vs_electric"; "gas_vs_electric"; "pga_0p1g"; ...
    "earthquake_year_10y"; "income_LowIncome"; "income_LowerMiddleIncome"; ...
    "income_UpperMiddleIncome"; "population_doubling"; "D0_doubling"];
T = orderedRows(T, "effect_id", order);
rowY = [12.4 11.4 9.4 8.4 6.4 5.4 4.4 2.4 .8];
labels = ["Water supply"; "Natural gas"; "Mean PGA"; ...
    "Earthquake year"; "Low income"; "Lower-middle income"; ...
    "Upper-middle income"; "Matched population"; "\itD_{0}\rm"];
hold(ax, "on");
ax.XScale = "log";
ax.XLim = [.30 20];
ax.YLim = [.1 13];
xline(ax, 1, "--", "Color", [.55 .55 .55], "LineWidth", .75);
for k = 1:height(T)
    row = T(k, :);
    offScale = row.ci95_low < ax.XLim(1) || row.ci95_high > ax.XLim(2);
    if ~offScale
        plot(ax, [row.ci95_low row.ci95_high], [rowY(k) rowY(k)], ...
            "Color", "k", "LineWidth", 1.1, "Tag", "FactorCI_" + order(k));
    end
    plot(ax, row.estimate, rowY(k), "ko", "MarkerFaceColor", "k", ...
        "MarkerSize", 5.2, "Tag", "FactorEstimate_" + order(k));
    if row.p_value < .001
        pText = "\itP\rm < 0.001";
    else
        pText = sprintf("\\itP\\rm = %.3f", row.p_value);
    end
    effectLabel = sprintf("%.3f\n%s", row.estimate, pText);
    if offScale
        effectLabel = sprintf("%.3f\n95%% CI off scale\n%s", row.estimate, pText);
    end
    % A separate annotation column avoids drawing text across long intervals.
    labelY = (rowY(k) - ax.YLim(1)) / diff(ax.YLim);
    text(ax, 1.035, labelY, effectLabel, "Units", "normalized", ...
        "FontName", "Arial", "FontSize", 7, "Color", [.18 .18 .18], ...
        "Interpreter", "tex", "HorizontalAlignment", "left", ...
        "VerticalAlignment", "middle", "Clipping", "off", ...
        "Tag", "FactorSummary_" + order(k));
end
ax.XTick = [.5 1 2 5 10 20];
ax.XTickLabel = ["0.5" "1" "2" "5" "10" "20"];
ax.YTick = fliplr(rowY);
ax.YTickLabel = flipud(labels);
xlabel(ax, "Adjusted \itT_{90}\rm factor");
styleAxes(ax, 7.2);
ax.YAxis.TickLength = [0 0];
ax.UserData = T;
end

function T = orderedRows(T, field, order)
indices = zeros(numel(order), 1);
for k = 1:numel(order)
    hit = find(T.(field) == order(k));
    if numel(hit) ~= 1
        error("ContextualFigure:ExpectedRow", "Expected one %s row for %s.", field, order(k));
    end
    indices(k) = hit;
end
T = T(indices, :);
end

function styleAxes(ax, fontSize)
ax.FontName = "Arial";
ax.TickLabelInterpreter = "tex";
ax.FontSize = fontSize;
ax.LineWidth = .75;
ax.XColor = [.18 .18 .18];
ax.YColor = [.18 .18 .18];
ax.TickDir = "out";
ax.TickLength = [.012 .012];
ax.Box = "off";
ax.XGrid = "off";
ax.YGrid = "off";
ax.XMinorGrid = "off";
ax.YMinorGrid = "off";
ax.Layer = "top";
end

function panelLetter(fig, position, label)
annotation(fig, "textbox", position, "String", label, "LineStyle", "none", ...
    "FontName", "Arial", "FontSize", 10, "FontWeight", "bold", ...
    "Margin", 0, "Color", [.18 .18 .18]);
end
