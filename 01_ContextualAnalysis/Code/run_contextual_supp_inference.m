function [systemResult, d0FactorResult] = run_contextual_supp_inference(T, cfg)
% Coefficient inference shared by the supplementary workflow and validation.
% Fixed models use CR2/Satterthwaite; random intercepts use model-based
% Satterthwaite inference. R2 bootstrap calculations stay in the caller.
[pairSpecs, systemSpecs] = specifications(cfg);
systemResult = fixedSystemFactors(T, systemSpecs, cfg);
d0FactorResult = fixedD0Factors(T, pairSpecs, cfg);
[mixedSystem, mixedD0] = mixedFactors(T, cfg);
systemResult = [systemResult; mixedSystem];
d0FactorResult = [d0FactorResult; mixedD0];
end

function [pairSpecs, systemSpecs] = specifications(cfg)
pairSpecs = repmat(baseSpec(), 0, 1);
pairSpecs(end+1) = spec("primary", "log2(D0) (primary)", "D0 representation", "T90", "avg_PGA", "raw", "population", "log2_d0", nan);
pairSpecs(end+1) = spec("raw_d0", "Raw D0", "D0 representation", "T90", "avg_PGA", "raw", "population", "raw_d0", nan);
pairSpecs(end+1) = spec("log2_d0_complete", "log2(D0) + complete-disruption indicator", "D0 representation", "T90", "avg_PGA", "raw", "population", "log2_d0_complete", nan);
pairSpecs(end+1) = spec("spline_log2_d0", "Spline in log2(D0)", "D0 representation", "T90", "avg_PGA", "raw", "population", "spline_log2_d0", nan);
pairSpecs(end+1) = spec("log2p1_d0_zero_inclusive", "log2(1 + D0)", "zero-inclusive D0", "T90", "avg_PGA", "raw", "population", "log2p1_d0", nan);
for out = ["T80", "T95"]
    pairSpecs(end+1) = spec(out, out, "restoration milestone", out, "avg_PGA", "raw", "population", "log2_d0", nan); %#ok<AGROW>
end
for field = string(cfg.alternative_pga(:))'
    pairSpecs(end+1) = spec(field, pgaLabel(field), "PGA summary", "T90", field, "raw", "population", "log2_d0", nan); %#ok<AGROW>
end
pairSpecs(end+1) = spec("ln1p_pga_over_0p1", "ln[1 + PGA / (0.1 g)]", "PGA functional form", "T90", "avg_PGA", "ln1p", "population", "log2_d0", nan);
pairSpecs(end+1) = spec("population_density", "Population density", "population measure", "T90", "avg_PGA", "raw", "density", "log2_d0", nan);
pairSpecs(end+1) = spec("exclude_bottom_1pct", "Exclude bottom 1% of positive D0", "small-positive D0", "T90", "avg_PGA", "raw", "population", "log2_d0", 0.01);
pairSpecs(end+1) = spec("exclude_bottom_5pct", "Exclude bottom 5% of positive D0", "small-positive D0", "T90", "avg_PGA", "raw", "population", "log2_d0", 0.05);

systemSpecs = repmat(baseSpec(), 0, 1);
systemSpecs(end+1) = spec("primary", "Primary model", "system-and-context model", "T90", "avg_PGA", "raw", "population", "none", nan);
for field = string(cfg.alternative_pga(:))'
    systemSpecs(end+1) = spec(field, pgaLabel(field), "PGA summary", "T90", field, "raw", "population", "none", nan); %#ok<AGROW>
end
systemSpecs(end+1) = spec("ln1p_pga_over_0p1", "ln[1 + PGA / (0.1 g)]", "PGA functional form", "T90", "avg_PGA", "ln1p", "population", "none", nan);
systemSpecs(end+1) = spec("population_density", "Population density", "population measure", "T90", "avg_PGA", "raw", "density", "none", nan);
systemSpecs(end+1) = spec("T80", "T80", "restoration milestone", "T80", "avg_PGA", "raw", "population", "none", nan);
systemSpecs(end+1) = spec("T95", "T95", "restoration milestone", "T95", "avg_PGA", "raw", "population", "none", nan);
end

function s = baseSpec()
s = struct("id", "", "label", "", "family", "", "outcome", "T90", "pga_field", "avg_PGA", ...
    "pga_transform", "raw", "context_kind", "population", "d0_kind", "none", "trim_quantile", nan, "model_kind", "linear");
end

function s = spec(id, label, family, outcome, pga, pgaTransform, context, d0, trim)
s = baseSpec(); s.id = string(id); s.label = string(label); s.family = string(family);
s.outcome = string(outcome); s.pga_field = string(pga); s.pga_transform = string(pgaTransform);
s.context_kind = string(context); s.d0_kind = string(d0); s.trim_quantile = trim;
end

function label = pgaLabel(field)
switch string(field)
    case "popWeightedPGA_g", label = "Population-weighted PGA";
    case "popMeanPGA_g", label = "Mean PGA across populated cells";
    case "pgaMedian_g", label = "Median PGA across populated cells";
    case "pgaMax_g", label = "Maximum PGA across populated cells";
    otherwise, label = replace(string(field), "_", " ");
end
end

function R = fixedSystemFactors(T, specs, cfg)
rows = cell(numel(specs), 18);
for j = 1:numel(specs)
    D = selectTime(T, specs(j), false); levels = incomeLevels(D, cfg);
    [X, names, ~] = timeDesign(D, specs(j), false, levels, []); y = log(D.(specs(j).outcome));
    A = [ones(height(D), 1), X]; C = zeros(size(A,2), 2);
    C(1 + find(names == "water_vs_electric", 1), 1) = 1;
    C(1 + find(names == "gas_vs_electric", 1), 2) = 1;
    inf = ctx_cr2_satterthwaite(A, y, D.eid, C, ["water"; "gas"]);
    rows(j,:) = {specs(j).id, specs(j).label, specs(j).family, height(D), numel(unique(D.eid)), ...
        exp(inf.estimate(1)), exp(inf.estimate(2)), exp(inf.ci95_low(1)), exp(inf.ci95_high(1)), ...
        exp(inf.ci95_low(2)), exp(inf.ci95_high(2)), inf.p_value(1), inf.p_value(2), ...
        inf.satterthwaite_df(1), inf.satterthwaite_df(2), nan, nan, "CR2 Satterthwaite"};
end
R = cell2table(rows, 'VariableNames', {'specification','label','family','record_count','earthquake_count', ...
    'water_factor','gas_factor','water_ci95_low','water_ci95_high','gas_ci95_low','gas_ci95_high', ...
    'water_p_value','gas_p_value','water_df','gas_df','water_bootstrap_successes','gas_bootstrap_successes','ci_method'});
end

function R = fixedD0Factors(T, specs, cfg)
rows = cell(numel(specs), 13);
for j = 1:numel(specs)
    D = selectTime(T, specs(j), true); levels = incomeLevels(D, cfg);
    [X, ~, knots] = timeDesign(D, specs(j), true, levels, []); y = log(D.(specs(j).outcome));
    lo = D; hi = D; lo.outage0(:) = 0.10; hi.outage0(:) = 0.20;
    Xlo = timeDesign(lo, specs(j), true, levels, knots); Xhi = timeDesign(hi, specs(j), true, levels, knots);
    contrast = [0; mean(Xhi - Xlo, 1)'];
    inf = ctx_cr2_satterthwaite([ones(height(D),1), X], y, D.eid, contrast, "D0_0.10_to_0.20");
    rows(j,:) = {specs(j).id, specs(j).label, specs(j).family, height(D), numel(unique(D.eid)), ...
        exp(inf.estimate), exp(inf.ci95_low), exp(inf.ci95_high), inf.p_value, ...
        inf.satterthwaite_df, nan, "CR2 Satterthwaite", "D0 0.10 to 0.20"};
end
R = cell2table(rows, 'VariableNames', {'specification','label','family','record_count','earthquake_count', ...
    'estimate','ci95_low','ci95_high','p_value','satterthwaite_df','bootstrap_successes','ci_method','reported_contrast'});
end

function [systemRow, d0Row] = mixedFactors(T, cfg)
s = spec("random_intercept", "Earthquake random intercept", "earthquake random-intercept model", "T90", "avg_PGA", "raw", "population", "none", nan);
D = selectTime(T, s, false); model = fitMixed(D, s, false, cfg);
[w, wl, wh, wp, wdf] = mixedFactor(model, "sys_cat_W");
[g, gl, gh, gp, gdf] = mixedFactor(model, "sys_cat_G");
systemRow = cell2table({"random_intercept","Earthquake random intercept","earthquake random-intercept model", ...
    height(D),numel(unique(D.eid)),w,g,wl,wh,gl,gh,wp,gp,wdf,gdf,nan,nan,"Wald Satterthwaite"}, ...
    'VariableNames', {'specification','label','family','record_count','earthquake_count','water_factor','gas_factor', ...
    'water_ci95_low','water_ci95_high','gas_ci95_low','gas_ci95_high','water_p_value','gas_p_value','water_df','gas_df', ...
    'water_bootstrap_successes','gas_bootstrap_successes','ci_method'});

s.d0_kind = "log2_d0"; D = selectTime(T, s, true); model = fitMixed(D, s, true, cfg);
[v, lo, hi, p, df] = mixedFactor(model, "d0_x");
d0Row = cell2table({"random_intercept","Earthquake random intercept","earthquake random-intercept model", ...
    height(D),numel(unique(D.eid)),v,lo,hi,p,df,nan,"Wald Satterthwaite","D0 0.10 to 0.20"}, ...
    'VariableNames', {'specification','label','family','record_count','earthquake_count','estimate','ci95_low','ci95_high', ...
    'p_value','satterthwaite_df','bootstrap_successes','ci_method','reported_contrast'});
end

function model = fitMixed(D, s, includeD0, cfg)
levels = incomeLevels(D, cfg); D.sys_cat = categorical(D.sys, ["E";"W";"G"]);
D.income_cat = categorical(string(D.development_level), levels); D.event_group = categorical(D.eid);
D.log_outcome = log(D.(s.outcome));
if s.pga_transform == "ln1p", D.pga_x = log1p(D.(s.pga_field)/0.1); else, D.pga_x = D.(s.pga_field)/0.1; end
D.year_x = (D.year-2000)/10;
if s.context_kind == "density", D.context_x = log2(D.popDensity_km2); else, D.context_x = log2(D.population); end
formula = "log_outcome ~ 1 + sys_cat + income_cat + pga_x + year_x + context_x";
if includeD0, D.d0_x = log2(D.outage0); formula = formula + " + d0_x"; end
model = fitlme(D, char(formula + " + (1|event_group)"), "FitMethod", "ML");
end

function [factor, low, high, p, df] = mixedFactor(model, name)
q = model.Coefficients; idx = string(q.Name) == name;
ci = coefCI(model, "Alpha", 0.05, "DFMethod", "satterthwaite");
H = zeros(1, height(q)); H(idx) = 1;
[p, ~, ~, df] = coefTest(model, H, 0, "DFMethod", "satterthwaite");
factor = exp(q.Estimate(idx)); low = exp(ci(idx,1)); high = exp(ci(idx,2));
end

function D = selectTime(T, s, aligned)
out = s.outcome;
if aligned
    if s.d0_kind == "log2p1_d0", D = T(T.t90_zero_inclusive_eligible, :);
    elseif out == "T90", D = T(T.t90_primary_eligible, :);
    else, D = T(exactOutcome(T,out) & T.valid_context & T.valid_d0 & T.outage0 > 0, :); end
else
    if out == "T90", D = T(T.t90_context_eligible, :);
    else, D = T(exactOutcome(T,out) & T.valid_context, :); end
end
D = filterContext(D, s);
if aligned && isfinite(s.trim_quantile)
    threshold = quantile(D.outage0, s.trim_quantile); D = D(D.outage0 > threshold, :);
end
end

function D = filterContext(D, s)
ok = isfinite(D.(s.pga_field)) & D.(s.pga_field) >= 0;
if s.context_kind == "density", ok = ok & isfinite(D.popDensity_km2) & D.popDensity_km2 > 0;
else, ok = ok & isfinite(D.population) & D.population > 0; end
D = D(ok, :);
end

function tf = exactOutcome(T, out)
c = "C" + extractAfter(out, "T"); b = out + "_bound_type";
tf = T.valid_system & isfinite(T.(out)) & T.(out) > 0 & T.(c) == 0 & lower(string(T.(b))) == "exact";
end

function X = baseDesign(D, s, levels)
if s.pga_transform == "ln1p", pga = log1p(D.(s.pga_field)/0.1); else, pga = D.(s.pga_field)/0.1; end
if s.context_kind == "density", context = log2(D.popDensity_km2); else, context = log2(D.population); end
X = [pga,(D.year-2000)/10,context,double(D.sys=="W"),double(D.sys=="G")];
for j = 2:numel(levels), X(:,end+1) = double(string(D.development_level)==levels(j)); end
end

function [X, names, knots] = timeDesign(D, s, includeD0, levels, knots)
X = baseDesign(D, s, levels); names = ["pga";"year";"context";"water_vs_electric";"gas_vs_electric"];
for j = 2:numel(levels), names(end+1,1) = "income_" + matlab.lang.makeValidName(levels(j)); end
if ~includeD0, knots = []; return; end
switch s.d0_kind
    case "log2_d0", X(:,end+1) = log2(D.outage0); names(end+1,1) = "d0_main";
    case "raw_d0", X(:,end+1) = D.outage0; names(end+1,1) = "d0_main";
    case "log2_d0_complete"
        X(:,end+1) = log2(D.outage0); X(:,end+1) = double(D.outage0==1);
        names(end+1:end+2,1) = ["d0_main";"d0_complete"];
    case "log2p1_d0", X(:,end+1) = log2(1+D.outage0); names(end+1,1) = "d0_main";
    case "spline_log2_d0"
        z = log2(D.outage0); if isempty(knots), knots = quantile(z,[.05 .35 .65 .95]); end
        basis = rcsBasis(z,knots); X = [X,basis];
        for k = 1:size(basis,2), names(end+1,1) = "d0_spline_" + k; end
    otherwise, error("ContextualSupplementaryInference:D0Kind", "Unknown D0 representation %s", s.d0_kind);
end
end

function basis = rcsBasis(x, knots)
knots = unique(knots);
if numel(knots) < 4, error("ContextualSupplementaryInference:Spline", "Insufficient spline knots."); end
last = knots(end); penult = knots(end-1); basis = x;
for j = 1:numel(knots)-2
    h = max(x-knots(j),0).^3 - ((last-knots(j))/(last-penult))*max(x-penult,0).^3 + ...
        ((penult-knots(j))/(last-penult))*max(x-last,0).^3;
    basis(:,end+1) = h/(last-knots(1))^2;
end
end

function levels = incomeLevels(D, cfg)
levels = sort(unique(string(D.development_level))); reference = string(cfg.reference_income_group);
if ~ismember(reference, levels)
    error("ContextualSupplementaryInference:IncomeReference", "Reference income group is absent.");
end
levels = [reference; levels(levels ~= reference)];
end
