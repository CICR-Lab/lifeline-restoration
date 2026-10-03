function run_figure_1()
% Render Figure 1 from the canonical standardized dataset.

close all force;

scriptDir = string(fileparts(mfilename("fullpath")));
fprintf("Running Figure 1 script:\n  %s\n", mfilename("fullpath"));

standardizedMatPath = fullfile(scriptDir, ...
    "standardized_dataset.mat");

countryShp = findNaturalEarthCountryShapefile( ...
    scriptDir);

outDir = fullfile(scriptDir, "outputs");
outBase = fullfile(outDir, "Figure1");
if ~isfolder(outDir)
    mkdir(outDir);
end

if ~isfile(standardizedMatPath)
    error( ...
        "Cannot find standardized dataset:\n%s\n" + ...
        "Request access at https://doi.org/10.5281/zenodo.23113551. " + ...
        "After author approval, download standardized_dataset.mat and " + ...
        "place it in 00_Dataset/ at the path above; see README.md.", ...
        standardizedMatPath);
end

if strlength(countryShp) == 0
    error( ...
        "Cannot find a complete Natural Earth country shapefile under:\n%s\n\n" + ...
        "Required files must be extracted and kept together:\n" + ...
        "  ne_10m_admin_0_countries.shp\n" + ...
        "  ne_10m_admin_0_countries.dbf\n" + ...
        "  ne_10m_admin_0_countries.shx", ...
        scriptDir);
end

fprintf( ...
    "Using Natural Earth country shapefile:\n  %s\n", ...
    countryShp);

Sstd = load( ...
    standardizedMatPath, ...
    "event_ISO", ...
    "event_region_system_records");

if ~isfield(Sstd, "event_ISO")
    error( ...
        "standardized_dataset.mat does not contain event_ISO:\n%s", ...
        standardizedMatPath);
end

stats = buildFigure1StatsFromStandardizedEventISO( ...
    Sstd.event_ISO);

if isfield(Sstd, "event_region_system_records")
    currentRecords = Sstd.event_region_system_records;
else
    currentRecords = struct([]);
end

[rawRecordBins, rawRecordTotals] = ...
    countObservationBinsFromStandardizedEventISO( ...
        Sstd.event_ISO, ...
        ["E", "W", "G"]);

countries = shaperead( ...
    countryShp, ...
    "UseGeoCoords", ...
    true);

countryCountMap = buildCountryCountMap(stats);
overlapCounts = countOverlapOnly(stats);

printFigure1DatasetTotals( ...
    stats, ...
    currentRecords);

fig = figure("Color", "w", "Units", "pixels", ...
    "Position", [45 25 1600 1320], "Resize", "off", ...
    "Renderer", "opengl", "Visible", "on");
set(fig, "WindowStyle", "normal");
set(fig, "DefaultAxesFontName", "Arial");
set(fig, "DefaultTextFontName", "Arial");
set(fig, "DefaultAxesFontSize", 16);
set(fig, "DefaultTextFontSize", 16);
set(fig, "DefaultAxesFontWeight", "normal");
set(fig, "DefaultTextFontWeight", "normal");

%% Panel a
panelAShiftX = -0.014;
panelAShiftY = 0.015;
fixedPanelAColorbarPosition = [0.7619 0.5697 0.0133 0.3500];

axMap = axes(fig, "Position", [0.002 + panelAShiftX 0.570 + panelAShiftY 0.795 0.350]);
hold(axMap, "on");
drawWorldCountries(axMap, countries, countryCountMap);
drawJapanBox(axMap);
drawEpicenters(axMap, stats);
axis(axMap, "equal");
xlim(axMap, [-180 180]);
ylim(axMap, [-58 84]);
axis(axMap, "off");
addPanelLabel(axMap, "a");

cb = colorbar(axMap, "eastoutside");
cb.Label.String = "";
cb.FontSize = 16;
cb.FontWeight = "normal";
cb.Position = fixedPanelAColorbarPosition;
maxCountryEvents = max(cell2mat(values(countryCountMap)));
clim(axMap, [-0.5 maxCountryEvents]);
if maxCountryEvents >= 10
    tickVals = unique(round([0, 10, (10 + maxCountryEvents)/2, maxCountryEvents]));
else
    tickVals = unique(round([0, maxCountryEvents]));
end
cb.Ticks = tickVals;
cb.TickLabels = string(tickVals);
colormap(axMap, datasetMapColormap(maxCountryEvents));
drawCountryRecordColorbarLabel(fig, cb);
drawMagnitudeSizeLegend(axMap, stats);
drawJapanInset(fig, countries, stats, countryCountMap, panelAShiftX, panelAShiftY);

%% Panel b
mapPlotLeft = mapVisibleLeft(axMap, fig);
panelLeft = max(mapPlotLeft + 0.004, 0.040);
panelRight = 0.804;
axB = axes(fig, "Position", [panelLeft 0.340 panelRight - panelLeft 0.175]);
drawYearLines(axB, stats);
addPanelLabel(axB, "b");

%% Panel c and d
bottomY = 0.062;
bottomH = 0.168;
gapCD = 0.046;
panelDWidth = 0.270;
panelCWidth = panelRight - panelLeft - gapCD - panelDWidth;
panelDLeft = panelLeft + panelCWidth + gapCD;

axC = axes(fig, "Position", [panelLeft bottomY panelCWidth bottomH]);
drawPointCountStacked(axC, rawRecordBins, rawRecordTotals, "Number of records");
addPanelLabel(axC, "c");

axD = axes(fig, "Position", [panelDLeft bottomY panelDWidth bottomH]);
drawOverlapUpSetOnly(axD, overlapCounts);
addPanelLabel(axD, "d");

alignVerticalAxisLabels(fig, [axB, axC], 0.010);

drawnow;
hideAxesToolbars(fig);
savefig(fig, outBase + ".fig");
exportFigureCompat(fig, outBase + ".png", "png", 600);
exportFigureCompat(fig, outBase + ".pdf", "pdf", 600);
fprintf("Saved figure:\n  %s.fig\n  %s.png\n  %s.pdf\n", outBase, outBase, outBase);

end

%% Local functions
function stats = buildFigure1StatsFromStandardizedEventISO(event_ISO)
    n = numel(event_ISO);
    seq = (1:n)';
    country_iso3 = strings(n,1);
    country_name = strings(n,1);
    yearVals = nan(n,1);
    lat = nan(n,1);
    lon = nan(n,1);
    magnitude = nan(n,1);
    n_region = zeros(n,1);
    n_E = zeros(n,1);
    n_W = zeros(n,1);
    n_G = zeros(n,1);
    n_EWG = zeros(n,1);
    n_EW = zeros(n,1);
    n_EG = zeros(n,1);
    n_WG = zeros(n,1);
    n_E_only = zeros(n,1);
    n_W_only = zeros(n,1);
    n_G_only = zeros(n,1);
    total_system_region = zeros(n,1);
    total_any_region = zeros(n,1);
    eventDate = NaT(n,1);
    decimalYearVals = nan(n,1);
    event_key = strings(n,1);
    event_title = strings(n,1);

    for i = 1:n
        Ei = event_ISO(i);
        country_name(i) = getFieldString(Ei, "Country");
        country_iso3(i) = normalizeIso3(getFieldString(Ei, "ISO"));
        magnitude(i) = getFieldNumber(Ei, "magnitude");
        event_title(i) = getFieldString(Ei, "usgsEventTitle");

        dt = parseDateCell(getFieldValue(Ei, "OccTime"));
        eventDate(i) = dt;
        if ~isnat(dt)
            yearVals(i) = year(dt);
            decimalYearVals(i) = decimalYear(dt);
        end

        [lat(i), lon(i)] = parseEpicenterCoordinate(getFieldString(Ei, "Epicenter"));
        event_key(i) = makeStandardizedEventKey(Ei, dt, magnitude(i), lat(i), lon(i));

        if ~isfield(Ei, "Region") || isempty(Ei.Region)
            continue;
        end

        n_region(i) = numel(Ei.Region);
        for j = 1:numel(Ei.Region)
            R = Ei.Region(j);
            hasE = hasSystemObservation(R, "Power");
            hasW = hasSystemObservation(R, "Water");
            hasG = hasSystemObservation(R, "Gas");

            n_E(i) = n_E(i) + double(hasE);
            n_W(i) = n_W(i) + double(hasW);
            n_G(i) = n_G(i) + double(hasG);
            total_any_region(i) = total_any_region(i) + double(hasE || hasW || hasG);

            if hasE && hasW && hasG
                n_EWG(i) = n_EWG(i) + 1;
            elseif hasE && hasW
                n_EW(i) = n_EW(i) + 1;
            elseif hasE && hasG
                n_EG(i) = n_EG(i) + 1;
            elseif hasW && hasG
                n_WG(i) = n_WG(i) + 1;
            elseif hasE
                n_E_only(i) = n_E_only(i) + 1;
            elseif hasW
                n_W_only(i) = n_W_only(i) + 1;
            elseif hasG
                n_G_only(i) = n_G_only(i) + 1;
            end
        end
        total_system_region(i) = n_E(i) + n_W(i) + n_G(i);
    end

    stats = table(seq, country_iso3, country_name, yearVals, lat, lon, magnitude, ...
        n_region, n_E, n_W, n_G, n_EWG, n_EW, n_EG, n_WG, ...
        n_E_only, n_W_only, n_G_only, total_system_region, total_any_region, ...
        eventDate, decimalYearVals, event_key, event_title, ...
        'VariableNames', {'seq','country_iso3','country_name','year','lat','lon','magnitude', ...
        'n_region','n_E','n_W','n_G','n_EWG','n_EW','n_EG','n_WG', ...
        'n_E_only','n_W_only','n_G_only','total_system_region','total_any_region', ...
        'eventDate','decimalYear','event_key','event_title'});
end
function printFigure1DatasetTotals(stats, currentRecords)
    keys = string(stats.event_key);
    keys = keys(strlength(keys) > 0);
    nUniqueEvents = numel(unique(keys));
    nEarthquakeCountry = height(stats);

    regionUnits = stats.total_any_region;
    regionUnits = regionUnits(isfinite(regionUnits));
    nEarthquakeRegionUnits = sum(regionUnits);

    if nargin >= 2 && ~isempty(currentRecords)
        nEarthquakeRegionSystem = numel(currentRecords);
    else
        records = stats.total_system_region;
        records = records(isfinite(records));
        nEarthquakeRegionSystem = sum(records);
    end

    fprintf("\nStandardized Figure 1 dataset totals:\n");
    fprintf("  %d distinct earthquake events\n", nUniqueEvents);
    fprintf("  %d earthquake-ISO records\n", nEarthquakeCountry);
    fprintf("  %d earthquake-region units with >=1 system\n", nEarthquakeRegionUnits);
    fprintf("  %d earthquake-region-system records\n", nEarthquakeRegionSystem);
    fprintf("\n");
end
function M = buildCountryCountMap(stats)
    valid = strlength(stats.country_iso3) > 0 & strlength(stats.event_key) > 0;
    canonIso = strings(height(stats), 1);
    for i = 1:height(stats)
        canonIso(i) = canonicalCountryIso3(stats.country_iso3(i));
    end
    valid = valid & strlength(canonIso) > 0;
    codes = unique(canonIso(valid));
    M = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for i = 1:numel(codes)
        hit = valid & canonIso == codes(i);
        keys = unique(stats.event_key(hit));
        keys = keys(strlength(keys) > 0);
        if isempty(keys)
            continue;
        end
        code = char(codes(i));
        if isKey(M, code)
            M(code) = M(code) + numel(keys);
        else
            M(code) = numel(keys);
        end
    end
end

function [C, totalBySys] = countObservationBinsFromStandardizedEventISO(event_ISO, sysOrder)
    C = zeros(numel(sysOrder), 3);
    totalBySys = zeros(numel(sysOrder), 1);
    sysLong = ["Power","Water","Gas"];
    for i = 1:numel(event_ISO)
        if ~isfield(event_ISO(i), "Region") || isempty(event_ISO(i).Region)
            continue;
        end
        for j = 1:numel(event_ISO(i).Region)
            R = event_ISO(i).Region(j);
            for s = 1:numel(sysOrder)
                nObs = systemObservationCount(R, sysLong(s));
                if nObs < 1
                    continue;
                end
                totalBySys(s) = totalBySys(s) + 1;
                if nObs <= 3
                    C(s,1) = C(s,1) + 1;
                elseif nObs <= 5
                    C(s,2) = C(s,2) + 1;
                else
                    C(s,3) = C(s,3) + 1;
                end
            end
        end
    end
end
function key = makeStandardizedEventKey(Ei, dt, mag, lat, lon)
    eqId = getFieldNumber(Ei, "event_ID");
    if isfinite(eqId)
        key = "eq" + string(sprintf("%06d", round(eqId)));
        return;
    end
    key = lower(strtrim(getFieldString(Ei, "usgsEventTitle")));
    if key == ""
        key = lower(strtrim(getFieldString(Ei, "EMDATDisNo")));
    end
    if key == ""
        if isnat(dt)
            d = "unknown-date";
        else
            d = string(datestr(dt, "yyyy-mm-dd HH:MM:SS"));
        end
        key = lower(sprintf('%s|%.2f|%.3f|%.3f', d, mag, lat, lon));
    end
end

function tf = hasSystemObservation(R, systemName)
    tf = systemObservationCount(R, systemName) >= 1;
end

function nObs = systemObservationCount(R, systemName)
    nObs = 0;
    systemName = char(systemName);
    if ~isfield(R, systemName) || isempty(R.(systemName))
        return;
    end
    S = R.(systemName);
    if isstruct(S) && isfield(S, "nObservations")
        nObs = double(S.nObservations);
        if ~isscalar(nObs) || ~isfinite(nObs)
            nObs = 0;
        end
        return;
    end

    if isstruct(S)
        curveFields = {"outagecountCurve", "functionalityCurve", "restorationratioCurve"};
        for k = 1:numel(curveFields)
            if isfield(S, curveFields{k})
                nObs = nObs + countCurvePointsForCoverage(S.(curveFields{k}));
            end
        end
        summaryFields = {"maxOutageCount", "maxOutageFraction"};
        for k = 1:numel(summaryFields)
            if isfield(S, summaryFields{k})
                v = S.(summaryFields{k});
                if isnumeric(v) && isscalar(v) && isfinite(v)
                    nObs = nObs + 1;
                end
            end
        end
    end
end

function n = countCurvePointsForCoverage(curve)
    n = 0;
    if isempty(curve)
        return;
    end
    if isstruct(curve)
        for k = 1:numel(curve)
            v = getFieldNumber(curve(k), "value");
            d = getFieldNumber(curve(k), "day");
            if isfinite(v) && (isfinite(d) || getFieldNumber(curve(k), "day_censor_code") ~= 0)
                n = n + 1;
            end
        end
    elseif isnumeric(curve) && size(curve, 2) >= 2
        n = sum(isfinite(curve(:,1)) & isfinite(curve(:,2)));
    end
end
function O = countOverlapOnly(stats)
    O.combos = ["E+W only","E+G only","W+G only","E+W+G"];
    O.counts = [sum(stats.n_EW), sum(stats.n_EG), sum(stats.n_WG), sum(stats.n_EWG)];
    O.totalUnits = sum(stats.total_any_region, "omitnan");
end

function value = getFieldValue(S, fieldName)
    if isfield(S, fieldName)
        value = S.(fieldName);
    else
        value = [];
    end
end

function s = getFieldString(S, fieldName)
    s = "";
    if ~isfield(S, fieldName), return; end
    v = S.(fieldName);
    if iscell(v)
        if isempty(v), return; end
        v = v{1};
    end
    if isempty(v), return; end
    s = string(v);
    if isempty(s) || ismissing(s), s = ""; else, s = strtrim(s(1)); end
end

function x = getFieldNumber(S, fieldName)
    x = NaN;
    if ~isfield(S, fieldName), return; end
    x = cellNumber(S.(fieldName));
end

function x = cellNumber(v)
    if iscell(v)
        if isempty(v), x = NaN; return; end
        v = v{1};
    end
    if isnumeric(v) && isscalar(v)
        x = double(v);
    elseif ischar(v) || isstring(v)
        x = str2double(string(v));
    else
        x = NaN;
    end
end

function dt = parseDateCell(v)
    dt = NaT;
    if iscell(v)
        if isempty(v), return; end
        v = v{1};
    end
    if isdatetime(v)
        dt = v;
    elseif isnumeric(v) && isscalar(v) && isfinite(v)
        try
            dt = datetime(v, "ConvertFrom", "excel");
        catch
            dt = NaT;
        end
    elseif ischar(v) || isstring(v)
        s = string(v);
        fmts = ["yyyy-MM-dd HH:mm:ss", "yyyy-MM-dd HH:mm", "yyyy-MM-dd", ...
            "MM/dd/yyyy HH:mm:ss", "MM/dd/yyyy"];
        for k = 1:numel(fmts)
            try
                dt = datetime(s, "InputFormat", fmts(k));
                if ~isnat(dt)
                    return;
                end
            catch
            end
        end
        try
            dt = datetime(s);
        catch
            dt = NaT;
        end
    end
end

function y = decimalYear(dt)
    y0 = dateshift(dt, "start", "year");
    y1 = dateshift(dt, "start", "year", "next");
    y = year(dt) + days(dt - y0) ./ days(y1 - y0);
end

function [lat, lon] = parseEpicenterCoordinate(s)
    lat = NaN;
    lon = NaN;

    s = string(s);

    if strlength(strtrim(s)) == 0 || ismissing(s)
        return;
    end


    expression = [ ...
        '([+-]?[0-9]*\.?[0-9]+)[^0-9A-Za-z+\-\.]*([NS])[,;\s]+' ...
        '([+-]?[0-9]*\.?[0-9]+)[^0-9A-Za-z+\-\.]*([EW])' ...
        ];

    tok = regexpi( ...
        char(s), ...
        expression, ...
        'tokens', ...
        'once');

    if isempty(tok)
        nums = regexp( ...
            char(s), ...
            '[-+]?[0-9]*\.?[0-9]+', ...
            'match');

        if numel(nums) >= 2
            lat = str2double(nums{1});
            lon = str2double(nums{2});
        end
        return;
    end

    lat = str2double(tok{1});
    lon = str2double(tok{3});

    if strcmpi(tok{2}, 'S')
        lat = -abs(lat);
    end

    if strcmpi(tok{4}, 'W')
        lon = -abs(lon);
    end
end

function countryShp = findNaturalEarthCountryShapefile(scriptDir)

    countryShp = "";

    fileName = "ne_10m_admin_0_countries.shp";

    searchDirs = [
        fullfile(scriptDir, "ne_10m_admin_0_countries")
        fullfile(scriptDir, "natural_earth_vector")];

    for d = 1:numel(searchDirs)
        if ~isfolder(searchDirs(d))
            continue;
        end

        directCandidate = fullfile(searchDirs(d), fileName);
        if isCompleteShapefile(directCandidate)
            countryShp = directCandidate;
            return;
        end

        matches = dir(fullfile(searchDirs(d), "**", fileName));
        for i = 1:numel(matches)
            candidate = string(fullfile(matches(i).folder, matches(i).name));

            if isCompleteShapefile(candidate)
                countryShp = candidate;
                return;
            end
        end
    end
end

function tf = isCompleteShapefile(shpPath)

    tf = false;

    if strlength(string(shpPath)) == 0 || ~isfile(shpPath)
        return;
    end

    [folder, name, ~] = fileparts(shpPath);

    tf = ...
        isfile(fullfile(folder, name + ".dbf")) && ...
        isfile(fullfile(folder, name + ".shx"));
end

function drawWorldCountries(ax, countries, countryCountMap)
    noEventColor = [0.90 0.90 0.88];
    edgeColor = [0.78 0.78 0.78];
    vals = cell2mat(values(countryCountMap));
    maxCount = max(vals);
    clim(ax, [-0.5 maxCount]);
    cmap = datasetMapColormap(maxCount);
    nC = size(cmap, 1);

    for k = 1:numel(countries)
        code = normalizeIso3(string(countries(k).ADM0_A3));
        if code == "" || code == "-99"
            code = normalizeIso3(string(countries(k).ISO_A3));
        end
        code = canonicalCountryIso3(code);

        if code ~= "" && isKey(countryCountMap, char(code))
            c = countryCountMap(char(code));
            idx = 1 + round((nC-1) * c / max(maxCount, eps));
            idx = max(1, min(nC, idx));
            fc = cmap(idx, :);
        else

            fc = noEventColor;
        end

        geoshow(ax, countries(k), "DisplayType", "polygon", ...
            "FaceColor", fc, "EdgeColor", edgeColor, "LineWidth", 0.25);
    end
end

function drawJapanCountry(ax, countries, countryCountMap)
    edgeColor = [0.45 0.45 0.45];
    vals = cell2mat(values(countryCountMap));
    maxCount = max(vals);
    clim(ax, [-0.5 maxCount]);
    cmap = datasetMapColormap(maxCount);
    nC = size(cmap, 1);

    for k = 1:numel(countries)
        code = normalizeIso3(string(countries(k).ADM0_A3));
        if code == "" || code == "-99"
            code = normalizeIso3(string(countries(k).ISO_A3));
        end
        code = canonicalCountryIso3(code);
        if code ~= "JPN"
            continue;
        end

        c = countryCountMap(char(code));
        idx = 1 + round((nC-1) * c / max(maxCount, eps));
        idx = max(1, min(nC, idx));
        geoshow(ax, countries(k), "DisplayType", "polygon", ...
            "FaceColor", cmap(idx, :), "EdgeColor", edgeColor, "LineWidth", 0.35);
    end
end

function drawJapanBox(ax)
    x = [126 147 147 126 126];
    y = [30 30 46 46 30];
    plot(ax, x, y, "k-", "LineWidth", 1.2);
end

function drawEpicenters(ax, stats)
    keep = isfinite(stats.lat) & isfinite(stats.lon);
    [~, uniqueIdx] = unique(string(stats.event_key(keep)), "stable");
    keepRows = find(keep);
    keep = false(height(stats), 1);
    keep(keepRows(uniqueIdx)) = true;
    markerSizes = magnitudeMarkerSizes(stats.magnitude(keep));
    scatter(ax, stats.lon(keep), stats.lat(keep), markerSizes, ...
        "Marker", "o", "MarkerFaceColor", [0.16 0.16 0.16], ...
        "MarkerFaceAlpha", 0.42, ...
        "MarkerEdgeColor", [0.92 0.92 0.92], ...
        "MarkerEdgeAlpha", 0.62, "LineWidth", 0.40);
end

function drawJapanInset(fig, countries, stats, countryCountMap, shiftX, shiftY)
    ax = axes(fig, "Position", [0.405 + shiftX 0.510 + shiftY 0.225 0.160]);
    ax.Color = "none";
    hold(ax, "on");
    drawJapanCountry(ax, countries, countryCountMap);
    colormap(ax, datasetMapColormap(max(cell2mat(values(countryCountMap)))));
    sub = stats(stats.lon >= 128 & stats.lon <= 147 & stats.lat >= 30 & stats.lat <= 46, :);
    if ~isempty(sub)
        drawEpicenters(ax, sub);
        legend(ax, "off");
    end
    axis(ax, "equal");
    xlim(ax, [128 147]);
    ylim(ax, [30 46]);
    ax.XTick = [];
    ax.YTick = [];
    ax.Box = "on";
    ax.LineWidth = 0.9;
    text(ax, 146.2, 30.8, "Japan", "FontSize", 16, "FontWeight", "normal", "Color", [0.12 0.12 0.12], ...
        "HorizontalAlignment", "right", "VerticalAlignment", "bottom");
end

function x = mapVisibleLeft(axMap, fig)
    pos = axMap.Position;
    xRange = diff(axMap.XLim);
    yRange = diff(axMap.YLim);
    figPos = fig.Position;
    targetPhysW = pos(4) * figPos(4) * xRange / yRange;
    targetNormW = targetPhysW / figPos(3);
    plotW = min(pos(3), targetNormW);
    x = pos(1) + max(0, (pos(3) - plotW) / 2);
end

function drawYearLines(ax, stats)
    [eventStats, Y] = aggregateStatsByEarthquakeEvent(stats);
    valid = isfinite(eventStats.decimalYear);
    [xRaw, ord] = sort(eventStats.decimalYear(valid));
    axisStart = min(xRaw);
    axisBreaks = selectYearAxisBreaks(xRaw);
    x = makeUniqueX(compressYearAxis(xRaw, axisStart, axisBreaks), 0.018);
    Y = Y(valid, :);
    Y = Y(ord, :);
    yrs = eventStats.year(valid);
    yrs = yrs(ord);

    cols = systemBaseColors();
    hold(ax, "on");
    for i = 1:numel(x)
        y0 = 0;
        for s = 1:3
            if Y(i, s) > 0
                line(ax, [x(i) x(i)], [y0 y0 + Y(i, s)], ...
                    "Color", cols(s, :), "LineWidth", 2.05);
                y0 = y0 + Y(i, s);
            end
        end
    end

    hE = plot(ax, NaN, NaN, "-", "Color", cols(1, :), "LineWidth", 5);
    hW = plot(ax, NaN, NaN, "-", "Color", cols(2, :), "LineWidth", 5);
    hG = plot(ax, NaN, NaN, "-", "Color", cols(3, :), "LineWidth", 5);
    ylabel(ax, "Number of records", "FontSize", 16);
    xlabel(ax, "Year", "FontSize", 16);
    lgd = legend(ax, [hE hW hG], ["Electric power", "Water supply", "Natural gas"], ...
        "Location", "none", ...
        "Orientation", "vertical", ...
        "NumColumns", 1, ...
        "Box", "off");
    lgd.FontSize = 16;
    lgd.FontWeight = "normal";

    lgd.Units = "normalized";
    pos = ax.Position;
    lgd.Position = [pos(1) + pos(3)*0.018, pos(2) + pos(4)*0.54, ...
        pos(3)*0.22, pos(4)*0.31];

    tickYears = unique([floor(min(yrs)) 1980 1990 2000 2010 2020 ceil(max(yrs))]);
    ax.XTick = compressYearAxis(tickYears, axisStart, axisBreaks);
    ax.XTickLabel = string(tickYears);
    ax.TickLength = [0 0];
    xlim(ax, [floor(min(x)) - 0.35, ceil(max(x)) + 0.35]);
    ax.Box = "off";
    ax.FontSize = 16;
    ylim(ax, [0, max(10, ceil(max(sum(Y, 2)) / 10) * 10)]);
    drawYearAxisBreakMarks(ax, axisBreaks, axisStart);
end

function [eventStats, Y] = aggregateStatsByEarthquakeEvent(stats)
    keys = string(stats.event_key);
    validKey = strlength(keys) > 0;
    keysUnique = unique(keys(validKey), "stable");
    eventStats = stats(false, :);
    Y = zeros(numel(keysUnique), 3);
    for i = 1:numel(keysUnique)
        hit = find(keys == keysUnique(i));
        if isempty(hit)
            continue;
        end
        first = hit(1);
        eventStats = [eventStats; stats(first,:)];
        Y(i,:) = [sum(stats.n_E(hit)), sum(stats.n_W(hit)), sum(stats.n_G(hit))];
        eventStats.n_E(end) = Y(i,1);
        eventStats.n_W(end) = Y(i,2);
        eventStats.n_G(end) = Y(i,3);
        eventStats.total_system_region(end) = sum(Y(i,:));
    end
end
function mag = distinctEventMagnitudes(stats)
    mag = double(stats.magnitude);
    keys = string(stats.event_key);
    keep = strlength(keys) > 0 & isfinite(mag);
    keys = keys(keep);
    mag = mag(keep);
    if isempty(keys)
        return;
    end
    [~, ia] = unique(keys, "stable");
    mag = mag(ia);
end

function drawCountryRecordColorbarLabel(fig, cb)
    ax = axes(fig, "Position", [0 0 1 1], "Color", "none", ...
        "XLim", [0 1], "YLim", [0 1], "Visible", "off", ...
        "HitTest", "off", "PickableParts", "none", ...
        "Tag", "CountryRecordColorbarLabel");

    cbPos = cb.Position;
    labelX = cbPos(1) + cbPos(3) + 0.032;
    labelY = cbPos(2) + cbPos(4)/2;

    text(ax, labelX, labelY, "Number of earthquakes", ...
        "Units", "normalized", ...
        "HorizontalAlignment", "center", ...
        "VerticalAlignment", "middle", ...
        "Rotation", 90, ...
        "FontName", "Arial", ...
        "FontSize", 16, ...
        "FontWeight", "normal", ...
        "Color", [0.12 0.12 0.12], ...
        "Interpreter", "none");
end

function drawMagnitudeSizeLegend(ax, stats)
    representative = [4.6 5.5 6.5 7.5 8.5 9.1];
    baseLabels = ["<5", "5-6", "6-7", "7-8", "8-9", "≥9"];

    mag = distinctEventMagnitudes(stats);
    mag = mag(isfinite(mag));
    [edges, ~, ~] = magnitudeBinScheme();
    counts = histcounts(mag, edges);

    countStrings = string(counts);

    legendX = -151;
    labelRightX = -137.5;

    leftParenX = -134.4;

    dy = 6.8;

    bottomMargin = 5;

    yBottom = ax.YLim(1) + bottomMargin;
    y = yBottom + (5:-1:0) * dy;

    mwGap = 8;
    mwY = y(1) + mwGap;

    scatter(ax, ...
        repmat(legendX,1,6), ...
        y, ...
        magnitudeMarkerSizes(representative), ...
        "o", ...
        "filled", ...
        "MarkerFaceColor", [0.16 0.16 0.16], ...
        "MarkerFaceAlpha", 0.52, ...
        "MarkerEdgeColor", [0.92 0.92 0.92], ...
        "MarkerEdgeAlpha", 0.78, ...
        "LineWidth", 0.40);

    for i = 1:6

        text(ax, ...
            labelRightX, ...
            y(i), ...
            baseLabels(i), ...
            "HorizontalAlignment", "right", ...
            "VerticalAlignment", "middle", ...
            "FontSize", 16, ...
            "FontWeight", "normal", ...
            "Color", [0.10 0.10 0.10]);

        text(ax, ...
            leftParenX, ...
            y(i), ...
            "(" + countStrings(i) + ")", ...
            "HorizontalAlignment", "left", ...
            "VerticalAlignment", "middle", ...
            "FontSize", 16, ...
            "FontWeight", "normal", ...
            "Color", [0.10 0.10 0.10]);
    end

    text(ax, ...
        -143.5, ...
        mwY, ...
        "Mw", ...
        "HorizontalAlignment", "center", ...
        "VerticalAlignment", "middle", ...
        "FontSize", 16, ...
        "FontWeight", "bold", ...
        "Color", [0.10 0.10 0.10]);
end

function cols = systemBaseColors()
    cols = [ ...
        0.75 0.22 0.25;   % E: electric power
        0.20 0.50 0.78;   % W: water supply
        0.86 0.55 0.16];  % G: natural gas
end

function cols = observationCountBinColors()
    cols = [ ...
        232 198 106
        143 188  90
         79 169 181] / 255;
end

function drawPointCountStacked(ax, C, totalBySys, xlabelText)
    sysNames = ["E","W","G"];
    binLabels = ["1-3", "4-5", "≥6"];
    segmentColors = observationCountBinColors();

    cla(ax, "reset");
    hold(ax, "on");
    y = 1 + (0:size(C,1)-1) * 0.50;

    b = barh(ax, y, C, "stacked", "EdgeColor", "none", "BarWidth", 0.44);
    for seg = 1:3
        b(seg).FaceColor = "flat";
        b(seg).CData = repmat(segmentColors(seg,:), numel(sysNames), 1);
    end

    maxTotal = max(sum(C, 2));
    if maxTotal <= 0
        maxTotal = 1;
    end
    for ii = 1:numel(sysNames)
        text(ax, sum(C(ii,:)) + maxTotal*0.012, y(ii), sprintf("%d", round(totalBySys(ii))), ...
            "HorizontalAlignment", "left", "VerticalAlignment", "middle", ...
            "FontSize", 16, "FontWeight", "normal", "Color", [0.10 0.10 0.10]);
    end

    set(ax, "YTick", y, "YTickLabel", sysNames, "YDir", "reverse", ...
        "TickLength", [0 0], "FontSize", 16);
    xlabel(ax, xlabelText, "FontSize", 16);
    ylabel(ax, "System", "FontSize", 16);
    xlim(ax, [0 maxTotal*1.13]);
    ylim(ax, [min(y)-0.28 max(y)+0.28]);
    ax.Box = "off";

    drawPointLegend(ax, segmentColors, binLabels);
end

function drawPointLegend(ax, colors, labels)
    fig = ancestor(ax, "figure");
    pos = ax.Position;
    legAx = axes(fig, "Position", [pos(1)+pos(3)*0.61 pos(2)+pos(4)*0.08 pos(3)*0.33 pos(4)*0.16]);
    axis(legAx, "off");
    xlim(legAx, [0 1]);
    ylim(legAx, [0 1]);
    hold(legAx, "on");

    x0 = 0.05;
    dx = 0.31;
    sw = 0.105;
    sh = 0.34;
    y0 = 0.34;
    for k = 1:numel(labels)
        thisX = x0 + (k-1)*dx;
        rectangle(legAx, ...
            "Position", [thisX y0 sw sh], ...
            "FaceColor", colors(k,:), ...
            "EdgeColor", "none");
        text(legAx, thisX + sw + 0.035, y0 + sh/2, labels(k), ...
            "FontName", "Arial", ...
            "FontSize", 16, ...
            "FontWeight", "normal", ...
            "HorizontalAlignment", "left", ...
            "VerticalAlignment", "middle", ...
            "Color", [0.10 0.10 0.10], ...
            "Interpreter", "none");
    end
end

function drawOverlapUpSetOnly(ax, O)
    combos = O.combos;
    counts = O.counts;
    nComb = numel(counts);
    cla(ax, "reset");
    hold(ax, "on");
bar(ax, 1:nComb, counts, 0.48, ...
    "FaceColor", [63 127 134] / 255, ...
    "EdgeColor", "none");
    for i = 1:nComb
        text(ax, i, counts(i) + max(counts) * 0.070, string(counts(i)), ...
            "HorizontalAlignment", "center", "FontSize", 16, "FontWeight", "normal");
    end

    nMulti = sum(counts);
    if isfield(O, "totalUnits") && isfinite(O.totalUnits) && O.totalUnits > 0
        summaryText = sprintf("%s of %s units (%.1f%%)", ...
            formatIntegerWithCommas(nMulti), ...
            formatIntegerWithCommas(O.totalUnits), ...
            100*nMulti/O.totalUnits);
        text(ax, 0.98, 0.88, summaryText, ...
            "Units", "normalized", ...
            "HorizontalAlignment", "right", "VerticalAlignment", "middle", ...
            "FontSize", 16, "FontWeight", "normal", ...
            "Color", [0.10 0.10 0.10]);
    end

    ylh = ylabel(ax, "Number of units", "FontSize", 16);
    ylh.Units = "normalized";
    yPos = ylh.Position;
    ylh.Position = [yPos(1) yPos(2)+0.06 yPos(3)];
    ax.XLim = [0.5, nComb + 0.5];
    ax.XTick = 1:nComb;
    ax.XTickLabel = combos;
    ax.TickLength = [0 0];
    ax.Box = "off";
    ax.FontSize = 16;
    xlabel(ax, "System combination", "FontSize", 16);
    ylim(ax, [0 max(counts) * 1.30]);
end

function alignVerticalAxisLabels(fig, axList, xNorm)
    axOverlay = axes(fig, "Position", [0 0 1 1], "Color", "none", ...
        "XLim", [0 1], "YLim", [0 1], "Visible", "off", ...
        "HitTest", "off", "PickableParts", "none", ...
        "Tag", "AlignedYAxisLabels");
    for i = 1:numel(axList)
        yl = axList(i).YLabel;
        labelStr = string(yl.String);
        labelFs = yl.FontSize;
        labelFw = yl.FontWeight;
        labelFn = yl.FontName;
        yl.String = "";
        pos = axList(i).Position;
        text(axOverlay, xNorm, pos(2) + pos(4)/2, labelStr, ...
            "Units", "normalized", ...
            "Rotation", 90, ...
            "HorizontalAlignment", "center", ...
            "VerticalAlignment", "middle", ...
            "FontSize", max(labelFs, 16), ...
            "FontWeight", labelFw, ...
            "FontName", labelFn, ...
            "Color", [0 0 0], ...
            "Interpreter", "none");
    end
end

function breaks = selectYearAxisBreaks(years)
    years = sort(unique(double(years(isfinite(years)))));
    gaps = diff(years);
    minGapYears = 5;
    maxBreaks = 2;
    displayWidth = 1.1;
    candidates = find(gaps >= minGapYears);
    if isempty(candidates)
        breaks = struct("start", {}, "end", {}, "displayWidth", {});
        return;
    end
    [~, order] = sort(gaps(candidates), "descend");
    keep = sort(candidates(order(1:min(maxBreaks, numel(order)))));
    breaks = repmat(struct("start", NaN, "end", NaN, ...
        "displayWidth", displayWidth), numel(keep), 1);
    for i = 1:numel(keep)
        breaks(i).start = years(keep(i));
        breaks(i).end = years(keep(i) + 1);
    end
end

function x = compressYearAxis(years, dataStart, breaks)
    years = double(years);

    if nargin < 2 || isempty(dataStart)
        dataStart = min(years);
    end
    if nargin < 3
        breaks = struct("start", {}, "end", {}, "displayWidth", {});
    end

    x = years;
    removed = 0;
    for i = 1:numel(breaks)
        b0 = breaks(i).start;
        b1 = breaks(i).end;
        w = breaks(i).displayWidth;
        shift = max(0, (b1 - b0) - w);
        inside = years > b0 & years < b1;
        after = years >= b1;
        x(inside) = b0 - removed + ...
            (years(inside) - b0) ./ max(eps, b1 - b0) .* w;
        x(after) = years(after) - removed - shift;
        removed = removed + shift;
    end
end

function drawYearAxisBreakMarks(ax, breaks, axisStart)
    if isempty(breaks)
        return;
    end
    y0 = ax.YLim(1);
    xr = range(ax.XLim);
    yr = range(ax.YLim);
    dx = 0.0028 * xr;
    pairGap = 0.0025 * xr;
    dy = 0.018 * yr;
    for i = 1:numel(breaks)
        edgeX = compressYearAxis([breaks(i).start breaks(i).end], ...
            axisStart, breaks);
        x0 = mean(edgeX);
        eraseHalfWidth = pairGap + 1.35 * dx;
        line(ax, [x0 - eraseHalfWidth, x0 + eraseHalfWidth], [y0 y0], "Color", "w", ...
            "LineWidth", 0.9, "Clipping", "off", "HandleVisibility", "off");
        line(ax, [x0 - pairGap - dx, x0 - pairGap + dx], ...
            [y0 - dy, y0 + dy], "Color", [0 0 0], ...
            "LineWidth", 0.9, "Clipping", "off", ...
            "HandleVisibility", "off");
        line(ax, [x0 + pairGap - dx, x0 + pairGap + dx], ...
            [y0 - dy, y0 + dy], "Color", [0 0 0], ...
            "LineWidth", 0.9, "Clipping", "off", ...
            "HandleVisibility", "off");
    end
end

function xOut = makeUniqueX(x, step)
    xOut = x;
    [u, ~, g] = unique(x);
    for k = 1:numel(u)
        idx = find(g == k);
        if numel(idx) > 1
            offsets = ((1:numel(idx)) - (numel(idx)+1)/2) * step;
            xOut(idx) = x(idx) + offsets(:);
        end
    end
end

function sizes = magnitudeMarkerSizes(magnitude)
    magnitude = double(magnitude(:));
    bins = discretize(magnitude, [-inf 5 6 7 8 9 inf]);
    sizeLevels = [16 28 44 63 86 114];
    sizes = repmat(sizeLevels(1), size(magnitude));
    valid = isfinite(bins);
    sizes(valid) = sizeLevels(bins(valid));
end

function textValue = formatIntegerWithCommas(value)
    digits = char(string(round(value)));
    parts = strings(0,1);
    while strlength(digits) > 3
        parts = [string(digits(end-2:end)); parts];
        digits = digits(1:end-3);
    end
    parts = [string(digits); parts];
    textValue = char(strjoin(parts, ","));
end

function cmap = datasetMapColormap(maxCount)
    % Panel-a count palette from the previous final Figure 1 script.

    if nargin < 1 || ~isfinite(maxCount) || maxCount <= 0
        maxCount = 1;
    end

    nMap = 256;
    anchors = [ ...
        0.98 0.92 0.76
        0.92 0.63 0.36
        0.80 0.32 0.43
        0.50 0.18 0.55
        0.18 0.17 0.49];
    x = linspace(0, 1, size(anchors, 1));
    xi = linspace(0, 1, nMap);
    cmap = [ ...
        interp1(x, anchors(:,1), xi, "linear")', ...
        interp1(x, anchors(:,2), xi, "linear")', ...
        interp1(x, anchors(:,3), xi, "linear")'];
end

function [edges, labels, cols] = magnitudeBinScheme()

    edges = [-inf 5 6 7 8 9 inf];
    labels = ["<5", "5-6", "6-7", "7-8", "8-9", "≥9"];

    cols = [ ...
        233 239 167
        196 222 151
        156 206 135
        94 178 101
        40 160 88
        7 140 77] / 255;
end

function addPanelLabel(ax, str)
    fig = ancestor(ax, "figure");
    pos = ax.Position;
    label = string(str);

    if any(label == ["a","b","c"])
        x = 0.001;
    elseif label == "d"
        x = pos(1) - 0.052;
    else
        x = pos(1) - 0.040;
    end

    if label == "a"
        y = pos(2) + pos(4) - 0.006;
    elseif label == "b"
        y = pos(2) + pos(4) + 0.014;
    elseif any(label == ["c","d"])
        y = pos(2) + pos(4) + 0.024;
    else
        y = pos(2) + pos(4) + 0.012;
    end

    annotation(fig, "textbox", [x y 0.03 0.03], ...
        "String", char(label), ...
        "EdgeColor", "none", ...
        "HorizontalAlignment", "left", ...
        "VerticalAlignment", "bottom", ...
        "FontName", "Arial", ...
        "FontSize", 16, ...
        "FontWeight", "bold", ...
        "Interpreter", "none", ...
        "FitBoxToText", "off");
end

function code = normalizeIso3(txt)
    code = upper(strtrim(string(txt)));
    code = regexprep(code, "[^A-Z0-9]", "");
    if strlength(code) > 3
        code = extractBefore(code, 4);
    end
end

function code = canonicalCountryIso3(code)
    code = normalizeIso3(code);
    if code == ""
        return;
    end
    switch code
        case ["ASM","GUM","PRI","MNP","VIR","UMI"]
            code = "USA";
        case ["MTQ","GLP","GUF","REU","MYT","BLM","MAF","PYF","NCL","WLF","SPM","ATF"]
            code = "FRA";
        case ["TWN","HKG","MAC"]
            code = "CHN";
        case ["COK","NIU","TKL"]
            code = "NZL";
    end
end

function hideAxesToolbars(fig)
    axList = findall(fig, "Type", "axes");
    for i = 1:numel(axList)
        if isprop(axList(i), "Toolbar") && ~isempty(axList(i).Toolbar)
            axList(i).Toolbar.Visible = "off";
        end
    end
end

function exportFigureCompat(fig, outPath, fmt, dpi)
    if nargin < 4 || ~isfinite(dpi)
        dpi = 600;
    end
    if ~isgraphics(fig, "figure")
        error("Figure handle is invalid before export: %s", outPath);
    end

    % Keep interactive and batch exports on the same fixed pixel canvas.
    set(fig, "Units", "pixels", "Resize", "off");
    pos = fig.Position;
    pos(3:4) = [1600 1320];
    fig.Position = pos;
    drawnow;

    try
        switch lower(string(fmt))
            case "png"
                exportgraphics(fig, outPath, "Resolution", dpi);
            case "pdf"
                exportgraphics(fig, outPath, "ContentType", "image", "Resolution", dpi);
            otherwise
                error("Unsupported format for exportFigureCompat: %s", fmt);
        end
        return;
    catch ME
        warning("exportgraphics failed for %s (%s). Falling back to print.", outPath, ME.message);
    end

    if ~isgraphics(fig, "figure")
        error("Figure handle became invalid before print fallback: %s", outPath);
    end

    set(fig, "PaperPositionMode", "auto");
    switch lower(string(fmt))
        case "png"
            print(fig, outPath, "-dpng", sprintf("-r%d", round(dpi)));
        case "pdf"
            print(fig, outPath, "-dpdf", sprintf("-r%d", round(dpi)));
        otherwise
            error("Unsupported format for print fallback: %s", fmt);
    end
end
