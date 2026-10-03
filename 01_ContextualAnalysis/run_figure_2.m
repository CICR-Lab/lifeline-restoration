function run_figure_2(varargin)
% Render Figure 2 from packaged or generated figure-data CSV files.
moduleRoot = string(fileparts(mfilename("fullpath")));
addpath(fullfile(moduleRoot, "Code"));
 [dataSource, outputDir, dataRoot] = parseFigureOptions(moduleRoot, varargin{:});
 dataDir = resolveFigureDataRoot(moduleRoot, dataSource, dataRoot);
 render_contextual_figures("", "", dataDir, outputDir);
end

function [dataSource, outputDir, dataDir] = parseFigureOptions(moduleRoot, varargin)
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
    error("ContextualFigure:InvalidDataSource", ...
        "DataSource must be 'packaged' or 'generated'.");
end
outputDir = string(p.Results.OutputDir);
if strlength(outputDir) == 0
    outputDir = fullfile(moduleRoot, "outputs", "figures", "main");
end
dataDir = string(p.Results.DataRoot);
end

function dataDir = resolveFigureDataRoot(moduleRoot, dataSource, dataRoot)
if strlength(dataRoot) > 0
    dataDir = dataRoot;
elseif dataSource == "packaged"
    dataDir = fullfile(moduleRoot, "figure_data", "main");
else
    dataDir = fullfile(moduleRoot, "outputs", "analysis", "figure_data");
end

required = ["figure2a_D0_records.csv", "figure2b_T90_records.csv", ...
    "figure2c_fitted_r2.csv", "figure2d_adjusted_T90_factors.csv"];
if any(~isfile(fullfile(dataDir, required)))
    error("ContextualFigure:MissingData", ...
        "Required Figure 2 data are missing from %s for DataSource='%s'.", ...
        dataDir, dataSource);
end
end

function tf = isTextScalar(value)
tf = (ischar(value) && (isrow(value) || isempty(value))) || ...
    (isstring(value) && isscalar(value));
end
