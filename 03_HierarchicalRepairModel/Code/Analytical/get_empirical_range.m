function range = get_empirical_range(T, metricName)
range = struct('outerLow',NaN,'centralLow',NaN,'median',NaN, ...
    'centralHigh',NaN,'outerHigh',NaN,'available',false,'metricName','');
if isempty(T), return; end
row = strcmpi(T.Metric,metricName);
if sum(row) ~= 1, return; end
range.metricName=char(metricName);
v = [T.OuterLow(row),T.CentralLow(row),T.Median(row),T.CentralHigh(row),T.OuterHigh(row)];
if any(~isfinite(v)) || ~(v(1)<=v(2) && v(2)<=v(3) && v(3)<=v(4) && v(4)<=v(5))
    return;
end
range.outerLow=v(1); range.centralLow=v(2); range.median=v(3);
range.centralHigh=v(4); range.outerHigh=v(5); range.available=true;
end
