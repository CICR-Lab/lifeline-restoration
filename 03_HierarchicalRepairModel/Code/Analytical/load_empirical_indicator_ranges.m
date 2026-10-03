function [T, isComplete] = load_empirical_indicator_ranges(filePath, requiredMetrics)

if nargin < 2 || isempty(requiredMetrics)
    requiredMetrics = {'tau80_tau50','tau90_tau80','tau95_tau90','kappa','eta90'};
end
if ~exist(filePath,'file')
    T = table(); isComplete = false; return;
end
phase = readtable(filePath,'Sheet','Fig3b_summary','TextType','string', ...
    'PreserveVariableNames',true);
loss = readtable(filePath,'Sheet','Fig3c_ratio_summary','TextType','string', ...
    'PreserveVariableNames',true);
late = readtable(filePath,'Sheet','Fig3d_summary','TextType','string', ...
    'PreserveVariableNames',true);

sourceMetrics = {'tau80_over_tau50','tau90_over_tau80','tau95_over_tau90'};
metricNames = {'tau80_tau50','tau90_tau80','tau95_tau90','kappa','eta90'};
values = nan(5,5);
for k = 1:3
    row = find_one(phase,1,sourceMetrics{k},2,'All systems',filePath);
    values(k,:) = phase{row,[7 5 4 6 8]};
end
row = find_one(loss,1,'All systems',0,'',filePath);
values(4,:) = loss{row,[6 4 3 5 7]};
row = find_one(late,1,'All systems',0,'',filePath);
values(5,:) = late{row,[6 4 3 5 7]};

T = table(string(metricNames(:)),values(:,1),values(:,2),values(:,3), ...
    values(:,4),values(:,5),'VariableNames', ...
    {'Metric','OuterLow','CentralLow','Median','CentralHigh','OuterHigh'});
isComplete = true;
for k = 1:numel(requiredMetrics)
    row = strcmpi(T.Metric,requiredMetrics{k});
    if sum(row) ~= 1
        isComplete = false; continue;
    end
    values = [T.OuterLow(row),T.CentralLow(row),T.Median(row),T.CentralHigh(row),T.OuterHigh(row)];
    if any(~isfinite(values)) || ~(values(1)<=values(2) && values(2)<=values(3) && values(3)<=values(4) && values(4)<=values(5))
        isComplete = false;
    end
end
end

function row = find_one(T,keyColumn,key,groupColumn,group,filePath)
mask = strcmpi(string(T{:,keyColumn}),key);
if groupColumn>0
    mask = mask & strcmpi(string(T{:,groupColumn}),group);
end
row = find(mask);
if numel(row)~=1
    error('Expected one %s / %s row in %s.',key,group,filePath);
end
end
