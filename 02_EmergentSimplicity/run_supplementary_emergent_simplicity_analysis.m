function run_supplementary_emergent_simplicity_analysis()

root = string(fileparts(mfilename("fullpath")));
outRoot = fullfile(root, "outputs", "supplementary");

ensureDir(outRoot);
ensureDir(fullfile(outRoot, "figures"));
ensureDir(fullfile(outRoot, "source_data"));

sourceDataDir = fullfile(outRoot, "source_data");
codeRoot = fullfile(root, "Code");

deleteIfFile(fullfile(sourceDataDir, "TableS8_S11_exception_classes.csv"));
deleteIfFile(fullfile(sourceDataDir, "TableS8_S11_exception_classes.xlsx"));
migrateLegacySourceData(sourceDataDir);
deleteLegacySourceData(sourceDataDir);

suppCode = fullfile(codeRoot, "supplementary");

addpath(suppCode);
make_TableS3_successive_milestone_ratio_resolution_sensitivity(outRoot);
plot_FigS6_well_resolved_trajectory_coverage(outRoot);
plot_FigS7_curve_resolution_sensitivity(outRoot);
plot_FigS8_S11_scurve_examples_exceptions(outRoot);

make_TableS4_loss_compression_specification_sensitivity(outRoot);
make_TableS5_late_loss_fraction_resolution_sensitivity(outRoot);

make_TableS6_cross_lifeline_timing_summary(outRoot);
make_TableS7_Cmax_summary(outRoot);

plot_FigS12(outRoot);

checkSupplementaryOutputs(outRoot);
exportSupplementaryFigureData(root, outRoot);

fprintf("\nSupplementary outputs saved under:\n  %s\n", outRoot);

end

function exportSupplementaryFigureData(moduleRoot, outRoot)
% Copy only direct figure inputs into the regenerable figure-data area.
figureDataDir = fullfile(moduleRoot, "outputs", "figure_data", "supplementary");
ensureDir(figureDataDir);

sourceFiles = [
    fullfile(outRoot, "source_data", "SourceData_FigS6.xlsx")
    fullfile(outRoot, "source_data", "SourceData_FigS7.xlsx")
    fullfile(outRoot, "source_data", "SourceData_FigS8_S11.xlsx")
    fullfile(outRoot, "source_data", "SourceData_FigS12.xlsx")
    fullfile(moduleRoot, "outputs", "intermediate", ...
        "01_restoration_shape_phase_ratio", "unique3_records.mat")
    fullfile(moduleRoot, "outputs", "intermediate", ...
        "01_restoration_shape_phase_ratio", "restoration_shape_ge5_fit_result.mat")
    fullfile(moduleRoot, "outputs", "intermediate", ...
        "01_restoration_shape_phase_ratio", ...
        "badfit_weibull_gompertz_loglogistic_ge6_candidates.xlsx")];

for i = 1:numel(sourceFiles)
    sourceFile = string(sourceFiles(i));
    if ~isfile(sourceFile)
        error("EmergentSupplementary:MissingFigureData", ...
            "Expected generated figure input was not found: %s", sourceFile);
    end
    [~, baseName, extension] = fileparts(char(sourceFile));
    copyfile(sourceFile, fullfile(figureDataDir, [baseName extension]), "f");
end
end

function checkSupplementaryOutputs(outRoot)

sourceDataDir = fullfile(outRoot, "source_data");
figureDir = fullfile(outRoot, "figures");

requiredFiles = [
    fullfile(sourceDataDir, "SourceData_TableS3.xlsx")
    fullfile(sourceDataDir, "SourceData_FigS6.xlsx")
    fullfile(sourceDataDir, "SourceData_FigS7.xlsx")
    fullfile(sourceDataDir, "SourceData_TableS4.xlsx")
    fullfile(sourceDataDir, "SourceData_TableS5.xlsx")
    fullfile(sourceDataDir, "SourceData_TableS6.xlsx")
    fullfile(sourceDataDir, "SourceData_TableS7.xlsx")
    fullfile(sourceDataDir, "SourceData_FigS8_S11.xlsx")
    fullfile(sourceDataDir, "SourceData_FigS12.xlsx")
    fullfile(figureDir, "FigS6_well_resolved_trajectory_coverage.png")
    fullfile(figureDir, "FigS6_well_resolved_trajectory_coverage.pdf")
    fullfile(figureDir, "FigS7_curve_resolution_sensitivity.png")
    fullfile(figureDir, "FigS7_curve_resolution_sensitivity.pdf")
    fullfile(figureDir, "FigS8_scurve_examples.png")
    fullfile(figureDir, "FigS8_scurve_examples.pdf")
    fullfile(figureDir, "FigS9_exception_long_tail.png")
    fullfile(figureDir, "FigS9_exception_long_tail.pdf")
    fullfile(figureDir, "FigS10_exception_near_threshold.png")
    fullfile(figureDir, "FigS10_exception_near_threshold.pdf")
    fullfile(figureDir, "FigS11_exception_unknown.png")
    fullfile(figureDir, "FigS11_exception_unknown.pdf")
    fullfile(figureDir, "FigS12_contextual_stability.png")
    fullfile(figureDir, "FigS12_contextual_stability.pdf")];

missing = requiredFiles(~isfile(requiredFiles));
if ~isempty(missing)
    error("Missing supplementary output file(s):\n%s", strjoin(missing, newline));
end

unwantedPaths = [
    fullfile(outRoot, "tables")
    fullfile(outRoot, "cohort_scope_summary.csv")
    fullfile(sourceDataDir, "TableS8_S11_exception_classes.csv")
    fullfile(sourceDataDir, "TableS8_S11_exception_classes.xlsx")
    fullfile(sourceDataDir, "contextual_stability_data.mat")
    fullfile(sourceDataDir, "contextual_stability_results.mat")
    fullfile(sourceDataDir, "Contextual_stability_results.xlsx")
    fullfile(sourceDataDir, "Contextual_contrasts_plot_source.csv")
    fullfile(sourceDataDir, "Contextual_R2_plot_source.csv")
    fullfile(sourceDataDir, "contextual_stability_data.xlsx")
    fullfile(figureDir, "ExtendedData_contextual_stability.png")
    fullfile(figureDir, "ExtendedData_contextual_stability.pdf")];

bad = unwantedPaths(arrayfun(@(p) isfolder(p) || isfile(p), unwantedPaths));
if ~isempty(bad)
    error("Obsolete supplementary output path(s) still exist:\n%s", strjoin(bad, newline));
end

end

function ensureDir(pathText)
if ~isfolder(pathText)
    mkdir(pathText);
end
end

function deleteIfFile(pathText)
if isfile(pathText)
    delete(pathText);
end
end

function migrateLegacySourceData(sourceDataDir)
pairs = [
    "SourceData_TableS3_successive_milestone_ratio_resolution_sensitivity.xlsx", "SourceData_TableS3.xlsx"
    "SourceData_TableS4_loss_compression_specification_sensitivity.xlsx", "SourceData_TableS4.xlsx"
    "SourceData_TableS5_late_loss_fraction_resolution_sensitivity.xlsx", "SourceData_TableS5.xlsx"
    "SourceData_TableS6_cross_lifeline_timing_summary.xlsx", "SourceData_TableS6.xlsx"
    "SourceData_TableS7_Cmax_summary.xlsx", "SourceData_TableS7.xlsx"
    "SourceData_FigS6_well_resolved_trajectory_coverage.xlsx", "SourceData_FigS6.xlsx"
    "SourceData_FigS7_curve_resolution_sensitivity.xlsx", "SourceData_FigS7.xlsx"
    "SourceData_FigS8_S11_scurve_examples_exceptions.xlsx", "SourceData_FigS8_S11.xlsx"
    "SourceData_FigS12_contextual_stability.xlsx", "SourceData_FigS12.xlsx"];

for i = 1:size(pairs, 1)
    oldPath = fullfile(sourceDataDir, pairs(i, 1));
    newPath = fullfile(sourceDataDir, pairs(i, 2));
    if isfile(oldPath) && ~isfile(newPath)
        movefile(oldPath, newPath);
    end
end
end

function deleteLegacySourceData(sourceDataDir)
oldNames = [
    "SourceData_TableS3_successive_milestone_ratio_resolution_sensitivity.xlsx"
    "SourceData_TableS4_loss_compression_specification_sensitivity.xlsx"
    "SourceData_TableS5_late_loss_fraction_resolution_sensitivity.xlsx"
    "SourceData_TableS6_cross_lifeline_timing_summary.xlsx"
    "SourceData_TableS7_Cmax_summary.xlsx"
    "SourceData_FigS6_well_resolved_trajectory_coverage.xlsx"
    "SourceData_FigS7_curve_resolution_sensitivity.csv"
    "SourceData_FigS7_curve_resolution_sensitivity.xlsx"
    "SourceData_FigS7.csv"
    "SourceData_FigS8_S11_scurve_examples_exceptions.xlsx"
    "SourceData_FigS12_contextual_stability.xlsx"];

for i = 1:numel(oldNames)
    deleteIfFile(fullfile(sourceDataDir, oldNames(i)));
end
end
