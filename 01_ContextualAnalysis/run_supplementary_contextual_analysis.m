function run_supplementary_contextual_analysis(varargin)
% Render Supplementary Figures 2-5 from packaged or generated CSV files.
moduleRoot = string(fileparts(mfilename("fullpath")));
addpath(fullfile(moduleRoot, "Code"));
 [dataSource, outputDir, dataRoot] = parseFigureOptions(moduleRoot, varargin{:});
 dataDir = resolveFigureDataRoot(moduleRoot, dataSource, dataRoot);
 render_contextual_supplement("", "", dataDir, outputDir);
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
p.addParameter("OutputDir", fullfile(moduleRoot, "outputs", "figures", "supplementary"), ...
    @(x) isTextScalar(x));
p.addParameter("DataRoot", "", @(x) isTextScalar(x));
p.parse(varargin{:});

dataSource = lower(char(string(p.Results.DataSource)));
if ~ismember(dataSource, {'packaged', 'generated'})
    error("ContextualSupplementary:InvalidDataSource", ...
        "DataSource must be 'packaged' or 'generated'.");
end
outputDir = string(p.Results.OutputDir);
if strlength(outputDir) == 0
    outputDir = fullfile(moduleRoot, "outputs", "figures", "supplementary");
end
dataDir = string(p.Results.DataRoot);
end

function dataDir = resolveFigureDataRoot(moduleRoot, dataSource, dataRoot)
if strlength(dataRoot) > 0
    dataDir = dataRoot;
elseif dataSource == "packaged"
    dataDir = fullfile(moduleRoot, "figure_data", "supplementary");
else
    dataDir = fullfile(moduleRoot, "outputs", "supplementary", "figure_data");
end

required = ["s1a_D0_fitted_R2_sensitivity.csv", ...
    "s1b_T90_fitted_R2_sensitivity.csv", ...
    "s1c_paired_fitted_delta_R2_sensitivity.csv", ...
    "s2a_primary_LOEO_performance.csv", ...
    "s2b_paired_LOEO_delta_sensitivity.csv", ...
    "s3a_system_factor_sensitivity.csv", ...
    "s3b_D0_factor_sensitivity.csv", ...
    "s4b_D0_calibration.csv", ...
    "s4_full_T90_residuals.csv", "s4_full_T90_qq.csv", ...
    "s4_aligned_plus_D0_residuals.csv", ...
    "s4_aligned_plus_D0_qq.csv"];
if any(~isfile(fullfile(dataDir, required)))
    error("ContextualSupplementary:MissingData", ...
        "Required supplementary figure data are missing from %s for DataSource='%s'.", ...
        dataDir, dataSource);
end
end

function tf = isTextScalar(value)
tf = (ischar(value) && (isrow(value) || isempty(value))) || ...
    (isstring(value) && isscalar(value));
end
