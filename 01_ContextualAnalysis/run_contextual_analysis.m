function run_contextual_analysis(projectRoot, runDir, nReplicates)
% Run the contextual analysis from the canonical dataset in Module 00.
if nargin < 1 || isempty(projectRoot) || strlength(string(projectRoot)) == 0
    projectRoot = string(fileparts(fileparts(fileparts(mfilename("fullpath")))));
end
codeDir = string(fileparts(mfilename("fullpath")));
codeRoot = string(fileparts(codeDir));
addpath(fullfile(codeDir, "Code"));
if nargin < 2 || isempty(runDir) || strlength(string(runDir)) == 0
    runDir = fullfile(codeDir, "outputs", "analysis");
end
if nargin < 3 || isempty(nReplicates)
    nReplicates = 2000;
end
projectRoot = string(projectRoot); runDir = string(runDir);
configPath = fullfile(codeDir, "contextual_config.json");
config = jsondecode(fileread(configPath));
sourceMat = fullfile(codeRoot, "00_Dataset", "standardized_dataset.mat");
extractedDir = fullfile(codeDir, "outputs", "intermediate", "extracted");
extractedMat = fullfile(extractedDir, "D0_Tp_broad_context.mat");
analysisInputDir = fullfile(codeDir, "outputs", "intermediate", "analysis_input");
canonicalPath = fullfile(analysisInputDir, "input", "canonical_analysis_table.mat");
if ~isfile(sourceMat)
    error("ContextualAnalysis:MissingSource", ...
        "Canonical standardized dataset is missing: %s\n" + ...
        "Request access at https://doi.org/10.5281/zenodo.23113551. " + ...
        "After author approval, download standardized_dataset.mat and " + ...
        "place it in 00_Dataset/ at the path above; see README.md.", sourceMat);
end
if isfolder(runDir)
    error("ContextualAnalysis:ExistingRun", "Analysis output already exists: %s", runDir);
end
config.project_root = projectRoot;
config.valid_systems = ["E"; "W"; "G"];
config.event_grouping = "event_ID";
config.event_iso_id_role = "event_ISO_ID";
if ~isfile(extractedMat)
    extract_broad_context(projectRoot, sourceMat, extractedDir);
end
ctx_build_inputs(config, extractedMat, analysisInputDir);
run_contextual_main_fit(projectRoot, canonicalPath, runDir, configPath, nReplicates, "");
supplementaryDir = fullfile(codeDir, "outputs", "supplementary");
run_contextual_supplementary(projectRoot, canonicalPath, configPath, runDir, supplementaryDir, nReplicates);
fprintf("Contextual analysis complete: %s\n", runDir);
end
