function run_emergent_simplicity_analysis()

emergentRoot = string(fileparts(mfilename("fullpath")));
analysisRoot = string(fileparts(emergentRoot));
inputMat = fullfile(analysisRoot, "00_Dataset", "standardized_dataset.mat");
outputDir = fullfile(emergentRoot, "outputs");

if ~isfolder(outputDir)
    mkdir(outputDir);
end

if ~isfile(inputMat)
    error("EmergentSimplicity:MissingSource", ...
        "Canonical standardized dataset is missing: %s\n" + ...
        "Request access at https://doi.org/10.5281/zenodo.23113551. " + ...
        "After author approval, download standardized_dataset.mat and " + ...
        "place it in 00_Dataset/ at the path above; see README.md.", inputMat);
end

runMainAnalysisScripts(emergentRoot);

figureDataDir = fullfile(outputDir, "figure_data", "main");
if ~isfolder(figureDataDir)
    mkdir(figureDataDir);
end
sourceDataOut = fullfile(figureDataDir, "SourceData_Figure3.xlsx");

writeFigure3SourceDataWorkbook(emergentRoot, sourceDataOut);

fprintf("Emergent simplicity analysis output saved:\n%s\n", sourceDataOut);

end

function runMainAnalysisScripts(emergentRoot)

codeRoot = fullfile(emergentRoot, "Code");

shapeCode = fullfile(codeRoot, ...
    "01_restoration_shape_phase_ratio");

lossCode = fullfile(codeRoot, ...
    "02_loss_compression_late_loss");

crossCode = fullfile(codeRoot, ...
    "03_cross_lifeline_timing_coupling");

pathsToAdd = char(join([shapeCode, lossCode, crossCode], pathsep));
addpath(pathsToAdd);
cleanupPath = onCleanup(@() rmpath(pathsToAdd));

fprintf("\nRunning run_00_prepare_restoration_shape_phase_ratio_data\n");
run_00_prepare_restoration_shape_phase_ratio_data();

fprintf("\nRunning run_01_fit_restoration_shape_models\n");
run_01_fit_restoration_shape_models();

fprintf("\nRunning run_00_prepare_loss_tail_dataset\n");
run_00_prepare_loss_tail_dataset();

fprintf("\nRunning run_01_make_loss_compression_late_loss_source_data\n");
run_01_make_loss_compression_late_loss_source_data();

fprintf("\nRunning run_01_make_cross_lifeline_timing_coupling_source_data\n");
run_01_make_cross_lifeline_timing_coupling_source_data();
end

function writeFigure3SourceDataWorkbook(emergentRoot, outFile)

intermediateRoot = fullfile(emergentRoot, "outputs", "intermediate");

shapeSource = fullfile(intermediateRoot, ...
    "01_restoration_shape_phase_ratio", ...
    "SourceData_Fig3a_restoration_shape.xlsx");

phaseSource = fullfile(intermediateRoot, ...
    "01_restoration_shape_phase_ratio", ...
    "SourceData_Fig3b_phase_ratios.xlsx");

lossSource = fullfile(intermediateRoot, ...
    "02_loss_compression_late_loss", ...
    "SourceData_Fig3c_loss_compression.xlsx");

lateLossSource = fullfile(intermediateRoot, ...
    "02_loss_compression_late_loss", ...
    "SourceData_Fig3d_tail_loss90.xlsx");

timingSource = fullfile(intermediateRoot, ...
    "03_cross_lifeline_timing_coupling", ...
    "SourceData_Fig3e_cross_lifeline_timing_ratios.xlsx");

couplingSource = fullfile(intermediateRoot, ...
    "03_cross_lifeline_timing_coupling", ...
    "SourceData_Fig3f_cross_lifeline_coupling_Cmax90.xlsx");

if isfile(outFile)
    delete(outFile);
end

sources = {
    shapeSource, "Fig3a_source_data", "Shape_model_summary"
    shapeSource, "Core_model_winner_rate", "Core_model_winner_rate"
    phaseSource, "Fig3a_normalized_trajectories", "Fig3a_normalized_trajectories"
    phaseSource, "record_level", "Fig3b_record_level"
    phaseSource, "summary", "Fig3b_summary"
    lossSource, "record_level", "Fig3c_record_level"
    lossSource, "fit_summary", "Fig3c_fit_summary"
    lossSource, "ratio_summary", "Fig3c_ratio_summary"
    lateLossSource, "record_level", "Fig3d_record_level"
    lateLossSource, "summary", "Fig3d_summary"
    timingSource, "record_level", "Fig3e_record_level"
    timingSource, "summary", "Fig3e_summary"
    couplingSource, "record_level", "Fig3f_record_level"
    couplingSource, "summary", "Fig3f_summary"};

for i = 1:size(sources, 1)
    sourceFile = string(sources{i, 1});
    sourceSheet = string(sources{i, 2});
    outSheet = string(sources{i, 3});

    assertFileExists(sourceFile);
    T = readtable(sourceFile, "Sheet", sourceSheet, ...
        "TextType", "string", "VariableNamingRule", "preserve");
    T = standardizeSystemNamesForOutput(T);
    writetable(T, outFile, "Sheet", outSheet);
end
end

function T = standardizeSystemNamesForOutput(T)
vars = string(T.Properties.VariableNames);
if any(vars == "group") && ~any(vars == "system")
    T.Properties.VariableNames{vars == "group"} = "system";
    vars = string(T.Properties.VariableNames);
end
for j = 1:numel(vars)
    name = lower(vars(j));
    if name == "system" || name == "lifeline_system"
        T.(char(vars(j))) = standardizeSystemText(T.(char(vars(j))));
    end
end
end

function x = standardizeSystemText(x)
if iscell(x)
    x = string(x);
end
if isstring(x) || iscategorical(x)
    s = string(x);
    s = strtrim(string(s));
    s(s == "All" | s == "Pooled" | s == "All systems") = "All systems";
    s(s == "Power" | s == "Electric power") = "Electric power";
    s(s == "Water" | s == "Water supply") = "Water supply";
    s(s == "Gas" | s == "Natural gas") = "Natural gas";

    known = ismissing(s) | s == "" | ismember( ...
        s, ...
        ["All systems", "Electric power", "Water supply", "Natural gas"]);
    if any(~known)
        error( ...
            "Unexpected system name(s) in source data: %s", ...
            strjoin(unique(s(~known)), ", "));
    end
    x = s;
end
end

function txt = formatPairSampleSize(T, pairVar, pairOrder)
    pairValues = string(T.(char(pairVar)));
    parts = strings(1, numel(pairOrder));

    for i = 1:numel(pairOrder)
        parts(i) = sprintf("n_%s=%d", ...
            pairSampleLabel(pairOrder(i)), nnz(pairValues == pairOrder(i)));
    end

    txt = strjoin(parts, "; ");
end

function label = pairSampleLabel(pairValue)
    switch string(pairValue)
        case "E->W"
            label = "E,W";
        case "E->G"
            label = "E,G";
        otherwise
            label = erase(string(pairValue), ["/", "->"]);
    end
end


function T = readtableChecked(filename, sheetName)
    assertFileExists(filename);
    if nargin < 2
        T = readtable(filename, "TextType", "string", ...
            "VariableNamingRule", "preserve");
    else
        T = readtable(filename, "Sheet", sheetName, ...
            "TextType", "string", "VariableNamingRule", "preserve");
    end
    T = normalizeTableStrings(T);
end

function T = normalizeTableStrings(T)
    for i = 1:width(T)
        if iscellstr(T.(i)) || isstring(T.(i)) || ischar(T.(i))
            T.(i) = strtrim(string(T.(i)));
        end
    end
end

function assertFileExists(pathValue)
    if ~isfile(pathValue)
        error("Required file not found:\n%s", pathValue);
    end
end

