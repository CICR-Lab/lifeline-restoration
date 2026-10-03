function ctx_build_inputs(config, extractedMat, runDir)
% Contextual analysis: canonical adapter for extracted  contextual data.
% Preserves the validated contextual-analysis design: event_ID grouping, positive-finite-T90-as-exact,
% D0_source audit-only. Does not modify existing pipeline functions.

root = string(config.project_root);
moduleRoot = ctx_module_root();
if nargin < 2 || isempty(extractedMat), extractedMat = fullfile(moduleRoot, "outputs", "intermediate", "extracted", "D0_Tp_broad_context.mat"); end
if nargin < 3 || isempty(runDir), runDir = fullfile(moduleRoot, "outputs", "intermediate", "analysis_input"); end
extractedMat = string(extractedMat); runDir = string(runDir);
if ~isfile(extractedMat), error("ContextualAnalysis:MissingSource", "Extracted  MAT file is missing: %s", extractedMat); end
for folder = ["input", "tables", "manifest"], if ~isfolder(fullfile(runDir, folder)), mkdir(fullfile(runDir, folder)); end, end
S = load(extractedMat, "D0andTp");
N = S.D0andTp;
required = ["record_ID", "event_ISO_ID", "event_ID", "Country", "ISO", "IncomeGroup", ...
    "OccTime", "magnitude", "regionIndex", "regionName", "System", "Depth", "Epicenter", ...
    "EMDATDisNo", "usgsEventTitle", "avg_PGA", "population", "popDensity_km2", ...
    "popWeightedPGA_g", "popMeanPGA_g", "pgaMedian_g", "pgaMax_g", "D0", "D0_source", ...
    "T10", "T20", "T30", "T40", "T50", "T60", "T70", "T80", "T90", "T95", "T100"];
missing = setdiff(required, string(N.Properties.VariableNames));
if ~isempty(missing), error("ContextualAnalysis:MissingColumns", "Missing  columns: %s", strjoin(missing, ", ")); end
% ---- Build canonical T with pipeline-compatible field layout ----
T = table();
T.record_id = strip(string(N.record_ID));
T.source_row = transpose(2:height(N)+1);
T.eid = strip(string(N.event_ID));
T.event_iso_id = double(N.event_ISO_ID);
T.rid = strip(string(N.regionName));
sys = upper(strip(string(N.System)));
sys(sys == "POWER") = "E";
sys(sys == "WATER") = "W";
sys(sys == "GAS") = "G";
T.sys = sys;
T.country = strip(string(N.Country));
T.ISO = strip(string(N.ISO));
if isdatetime(N.OccTime)
    occTime = N.OccTime;
else
    if isdatetime(N.OccTime)
    occTime = N.OccTime;
else
    occTime = datetime(strip(string(N.OccTime)), "InputFormat", "yyyy-MM-dd HH:mm:ss");
end
end
T.year = year(occTime);
T.magnitude = double(N.magnitude);

numericCtx = ["avg_PGA","population","popDensity_km2","popWeightedPGA_g", ...
    "popMeanPGA_g","pgaMedian_g","pgaMax_g"];
for name = numericCtx
    T.(name) = double(N.(name));
end

T.development_level = strip(string(N.IncomeGroup));
T.outage0 = double(N.D0);
T.outage0_source = strip(string(N.D0_source));

for threshold = [80 90 95]
    timeName = "T" + threshold;
    censorName = "C" + threshold;
    boundName = timeName + "_bound_type";
    val = double(N.(timeName));
    positiveFinite = isfinite(val) & val > 0;
    censor = nan(height(T),1);
    censor(positiveFinite) = 0;
    bound = repmat("missing", height(T), 1);
    bound(positiveFinite) = "exact";
    T.(timeName) = val;
    T.(censorName) = censor;
    T.(boundName) = bound;
end

% Retain extra T10-T100 for supplementary/audit only (not used in P1-P3 main pipeline).
for threshold = [10 20 30 40 50 60 70 100]
    name = "T" + threshold;
    T.(name) = double(N.(name));
end

% Development-level audit (self-contained; IncomeGroup embedded).
developmentAudit = table(T.eid, T.country, T.ISO, T.development_level, ...
    (strlength(T.development_level)>0), ...
    'VariableNames', {'eid','country','iso','embedded_development_level','develop_ok'});
writetable(developmentAudit, fullfile(runDir, "input", "development_level_audit.csv"));

T.crosswalk_matched = strlength(T.development_level) > 0;

% ---- Record-set flags (identical rules to ctx_build_inputs.m) ----
validSystems = string(config.valid_systems(:));
T.valid_system = ismember(T.sys, validSystems);
T.valid_d0 = isfinite(T.outage0) & T.outage0 >= 0 & T.outage0 <= 1;
T.valid_context = isfinite(T.year) & isfinite(T.avg_PGA) & T.avg_PGA >= 0 & ...
    isfinite(T.population) & T.population > 0 & ...
    strlength(T.development_level) > 0 & strlength(T.eid) > 0;

T.log_pga = nan(height(T),1);
T.log_pop_density = nan(height(T),1);
T.log_d0 = nan(height(T),1);
T.log1p_d0 = nan(height(T),1);
T.log_t90 = nan(height(T),1);
T.log_pga(T.avg_PGA > 0) = log(T.avg_PGA(T.avg_PGA > 0));
T.log_pop_density(T.population > 0) = log(T.population(T.population > 0));
T.log_d0(T.outage0 > 0) = log(T.outage0(T.outage0 > 0));
T.log1p_d0(T.valid_d0) = log1p(T.outage0(T.valid_d0));
T.log_t90(T.T90 > 0) = log(T.T90(T.T90 > 0));

T.d0_descriptive_eligible = T.valid_system & T.valid_d0;
T.d0_model_eligible = T.d0_descriptive_eligible & T.valid_context;
T.t80_descriptive_eligible = exactOutcome(T.T80, T.C80, T.T80_bound_type);
T.t90_descriptive_eligible = T.valid_system & exactOutcome(T.T90, T.C90, T.T90_bound_type);
T.t95_descriptive_eligible = exactOutcome(T.T95, T.C95, T.T95_bound_type);
T.t90_context_eligible = T.t90_descriptive_eligible & T.valid_context;
T.t90_primary_eligible = T.t90_descriptive_eligible & T.valid_context & T.valid_d0 & T.outage0 > 0;
T.t90_zero_inclusive_eligible = T.t90_descriptive_eligible & T.valid_context & T.valid_d0;
T.fold_eid = T.eid;

T.d0_exclusion_reason = d0Reason(T);
T.t90_exclusion_reason = t90Reason(T);
T.primary_exclusion_reason = primaryReason(T);

if numel(unique(T.record_id)) ~= height(T)
    error("ContextualAnalysis:DuplicateRecordId", "record_id is not unique in .");
end

scientificKey = T.eid + "|" + T.rid + "|" + T.sys;
duplicateScientificKey = duplicatedScientificKey(scientificKey);
keyAudit = T(duplicateScientificKey, ["record_id","source_row","eid","rid","sys","country"]);
writetable(keyAudit, fullfile(runDir, "input", "scientific_key_duplicates.csv"));

writetable(T, fullfile(runDir, "input", "canonical_analysis_table.csv"));
save(fullfile(runDir, "input", "canonical_analysis_table.mat"), "T", "-v7.3");
writetable(T(:, ["record_id","source_row","eid","event_iso_id","rid","sys","country", ...
    "outage0_source","crosswalk_matched","d0_descriptive_eligible", ...
    "d0_model_eligible","t90_descriptive_eligible","t90_primary_eligible", ...
    "t90_zero_inclusive_eligible","d0_exclusion_reason","t90_exclusion_reason", ...
    "primary_exclusion_reason"]), fullfile(runDir, "input", "record_set_manifest.csv"));

summary = recordSetSummary(T);
writetable(summary, fullfile(runDir, "input", "record_set_summary.csv"));

primaryComposition = T(T.t90_primary_eligible,:);
[Gc,countryLabel] = findgroups(primaryComposition.country);
countryRecords = splitapply(@numel, primaryComposition.record_id, Gc);
countryEvents = splitapply(@(x) numel(unique(x)), primaryComposition.eid, Gc);
countryComposition = table(countryLabel, countryRecords, countryEvents, ...
    'VariableNames', {'country','record_count','event_count'});
countryComposition = sortrows(countryComposition, 'record_count', 'descend');
writetable(countryComposition, fullfile(runDir, "tables", "table_s1_country_composition.csv"));

% Unmatched event crosswalk is not applicable (event_ID embedded); keep an empty audit file.
unmatched = T(~T.crosswalk_matched, ["record_id","eid","country"]);
writetable(unmatched, fullfile(runDir, "input", "unmatched_event_crosswalk.csv"));

% ---- Input manifest ----
source = struct();
source.mat = char(extractedMat);
source.mat_sha256 = fileHash(extractedMat);
source.source_format = "d0_tp_broad_context";
source.source_rows = height(T);
source.schema_fields = cellstr(required);
source.context_source_mode = "self_contained";
source.predictor_variant = "population";
source.population_field = "population";
source.event_grouping = char(config.event_grouping);
source.event_iso_id_role = char(config.event_iso_id_role);
writeJson(fullfile(runDir, "manifest", "input_manifest.json"), source);
fprintf("Contextual analysis complete: %d rows; D0 model n=%d; T90 primary n=%d.\n", ...
    height(T), sum(T.d0_model_eligible), sum(T.t90_primary_eligible));
end

% ---- Helper functions (copied from ctx_build_inputs.m to keep adapter independent) ----

function value = exactOutcome(time, censor, boundType)
value = isfinite(time) & time > 0 & isfinite(censor) & censor == 0 & ...
    lower(strip(string(boundType))) == "exact";
end

function key = normalizeKey(value)
key = lower(regexprep(strip(string(value)), "\s+", " "));
end

function result = duplicatedScientificKey(value)
[~,~,g] = unique(value);
counts = accumarray(g, 1);
result = counts(g) > 1;
end

function reason = d0Reason(T)
reason = repmat("eligible", height(T), 1);
reason(~T.valid_system) = "invalid_system";
reason(T.valid_system & ~T.valid_d0) = "missing_or_out_of_range_d0";
reason(T.d0_descriptive_eligible & ~T.valid_context) = "incomplete_context";
end

function reason = t90Reason(T)
reason = repmat("eligible", height(T), 1);
reason(~T.valid_system) = "invalid_system";
reason(T.valid_system & ~(isfinite(T.T90) & T.T90 > 0)) = "missing_or_nonpositive_t90";
reason(T.valid_system & isfinite(T.T90) & T.T90 > 0 & ...
    ~(isfinite(T.C90) & T.C90 == 0)) = "censored_t90";
reason(T.valid_system & isfinite(T.T90) & T.T90 > 0 & isfinite(T.C90) & T.C90 == 0 & ...
    lower(strip(T.T90_bound_type)) ~= "exact") = "nonexact_t90";
end

function reason = primaryReason(T)
reason = T.t90_exclusion_reason;
reason(T.t90_descriptive_eligible & ~T.valid_context) = "incomplete_context";
reason(T.t90_descriptive_eligible & T.valid_context & ~T.valid_d0) = ...
    "missing_or_out_of_range_d0";
reason(T.t90_descriptive_eligible & T.valid_context & T.valid_d0 & ~(T.outage0 > 0)) = ...
    "nonpositive_d0";
reason(T.t90_primary_eligible) = "eligible";
end

function output = recordSetSummary(T)
names = ["source","d0_descriptive","d0_model","t90_descriptive", ...
    "t90_primary","t90_zero_inclusive"];
flags = {true(height(T),1), T.d0_descriptive_eligible, T.d0_model_eligible, ...
    T.t90_descriptive_eligible, T.t90_primary_eligible, T.t90_zero_inclusive_eligible};
rows = cell(0,6);
for i = 1:numel(names)
    subset = T(flags{i},:);
    rows(end+1,:) = {names(i), height(subset), numel(unique(subset.eid)), ...
        sum(subset.sys=="E"), sum(subset.sys=="W"), sum(subset.sys=="G")}; %#ok<AGROW>
end
output = cell2table(rows, 'VariableNames', ...
    {'record_set_id','record_count','event_count','electricity_count','water_count','gas_count'});
end

function hash = fileHash(path)
fid = fopen(path, 'r');
if fid < 0, error("ContextualAnalysis:HashRead", "Unable to read %s", path); end
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, '*uint8')';
digest = java.security.MessageDigest.getInstance("SHA-256");
digest.update(bytes);
hash = lower(reshape(dec2hex(typecast(digest.digest(), 'uint8'))', 1, []));
end
function writeJson(path, value)
fid = fopen(path, "w");
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, jsonencode(value, "PrettyPrint", true), "char");
end
