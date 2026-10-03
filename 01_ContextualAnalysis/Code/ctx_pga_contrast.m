function [delta, label] = ctx_pga_contrast(transform)
% Reporting contrasts are defined on the corresponding model feature scale.
switch string(transform)
    case "raw"
        delta = 0.1;
        label = "0.1 g service-region mean PGA";
    case "ln1p_pga_over_0p1"
        delta = log(2);
        label = "PGA from 0 to 0.1 g";
    otherwise
        error("ContextualAnalysis:UnknownPGATransform", ...
            "Unknown PGA transformation: %s", string(transform));
end
end
