function [lo,hi] = ctx_reuse_loeo_baseline(D,y,contextPrediction,plusPrediction,cfg,baseline,kind)
% Reuse a primary interval only after checking the complete prediction identity.
if ~isequaln(cfg,baseline.config)
    error("ContextualLOEO:BaselineMismatch","Primary and repeated baseline configurations differ.");
end
if kind=="d0"
    source=baseline.d0;
    expected=source.prediction;
    actual=contextPrediction;
    lo=source.summary.r2_ci95_low;
    hi=source.summary.r2_ci95_high;
elseif kind=="paired"
    source=baseline;
    expected=[source.context_prediction,source.plus_prediction];
    actual=[contextPrediction,plusPrediction];
    lo=source.summary.delta_ci95_low;
    hi=source.summary.delta_ci95_high;
else
    error("ContextualLOEO:BaselineKind","Unknown baseline kind: %s",kind);
end
if ~isequal(D.record_id,source.record_id) || ~isequal(D.eid,source.eid) || ...
        ~isequal(y,source.observed) || ~isequal(size(actual),size(expected))
    error("ContextualLOEO:BaselineMismatch","Records, earthquake IDs or response scale differ from the primary analysis.");
end
if any(~isfinite(actual),"all") || any(~isfinite(expected),"all") || ...
        any(abs(actual-expected)>1e-10*(1+abs(expected)),"all")
    error("ContextualLOEO:BaselineMismatch","Out-of-fold predictions differ from the primary analysis.");
end
if ~isscalar(lo) || ~isscalar(hi) || ~isfinite(lo) || ~isfinite(hi) || lo>hi
    error("ContextualLOEO:BaselineInterval","The primary confidence interval is invalid.");
end
end
