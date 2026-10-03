function run_supplementary_emergent_simplicity(varargin)
% Render Supplementary Figures S6-S12 from packaged or generated inputs.

moduleRoot = string(fileparts(mfilename("fullpath")));
[dataSource, outputDir, dataRoot] = parseFigureOptions(moduleRoot, varargin{:});
if strlength(dataRoot) == 0
    if dataSource == "packaged"
        figureDataRoot = fullfile(moduleRoot, "figure_data", "supplementary");
    else
        figureDataRoot = fullfile(moduleRoot, "outputs", "figure_data", "supplementary");
    end
else
    figureDataRoot = dataRoot;
end

required = ["SourceData_FigS6.xlsx", "SourceData_FigS7.xlsx", ...
    "SourceData_FigS8_S11.xlsx", "SourceData_FigS12.xlsx", ...
    "unique3_records.mat", "restoration_shape_ge5_fit_result.mat", ...
    "badfit_weibull_gompertz_loglogistic_ge6_candidates.xlsx"];
if any(~isfile(fullfile(figureDataRoot, required)))
    error("EmergentSupplementary:MissingData", ...
        "Required supplementary figure data are missing from %s for DataSource='%s'.", ...
        figureDataRoot, dataSource);
end

if ~isfolder(figureDataRoot)
    error("Supplementary figure-data directory not found: %s", figureDataRoot);
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

addpath(fullfile(moduleRoot, "Code", "supplementary"));

plot_FigS6_well_resolved_trajectory_coverage(outputDir, figureDataRoot);
plot_FigS7_curve_resolution_sensitivity(outputDir, figureDataRoot);
plot_FigS8_S11_scurve_examples_exceptions(outputDir, figureDataRoot);
plot_FigS12(outputDir, figureDataRoot);

fprintf("\nSupplementary Figures S6-S12 saved under:\n  %s\n", outputDir);

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
p.addParameter("OutputDir", fullfile(moduleRoot, "outputs", "figures", "supplementary"), ...
    @(x) isTextScalar(x));
p.addParameter("DataRoot", "", @(x) isTextScalar(x));
p.parse(varargin{:});

dataSource = lower(char(string(p.Results.DataSource)));
if ~ismember(dataSource, {'packaged', 'generated'})
    error("EmergentSupplementary:InvalidDataSource", ...
        "DataSource must be 'packaged' or 'generated'.");
end
outputDir = string(p.Results.OutputDir);
if strlength(outputDir) == 0
    outputDir = fullfile(moduleRoot, "outputs", "figures", "supplementary");
end
dataRoot = string(p.Results.DataRoot);
end

function tf = isTextScalar(value)
tf = (ischar(value) && (isrow(value) || isempty(value))) || ...
    (isstring(value) && isscalar(value));
end
