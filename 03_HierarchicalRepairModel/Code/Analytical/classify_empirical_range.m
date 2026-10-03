function C = classify_empirical_range(values, range)
C = nan(size(values));
if ~isfield(range,'available') || ~range.available
    return;
end
valid = isfinite(values);
insideCentral = valid & values >= range.centralLow & values <= range.centralHigh;
if ~isfinite(range.outerLow) || ~isfinite(range.outerHigh)
    return;
end
insideOuter = valid & values >= range.outerLow & values <= range.outerHigh;
C(valid & ~insideOuter) = 3;
C(insideOuter & ~insideCentral) = 2;
C(insideCentral) = 1;
end
