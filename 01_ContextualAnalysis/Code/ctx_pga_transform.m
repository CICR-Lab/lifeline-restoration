function x = ctx_pga_transform(pga, transform)
% PGA feature definitions used throughout the  contextual analysis.
pga = double(pga);
if any(~isfinite(pga) | pga < 0)
    error("ContextualAnalysis:InvalidPGA", "PGA values must be finite and nonnegative.");
end

switch string(transform)
    case "raw"
        x = pga;
    case "ln1p_pga_over_0p1"
        x = log(1 + pga ./ 0.1);
    otherwise
        error("ContextualAnalysis:UnknownPGATransform", ...
            "Unknown PGA transformation: %s", string(transform));
end
end
