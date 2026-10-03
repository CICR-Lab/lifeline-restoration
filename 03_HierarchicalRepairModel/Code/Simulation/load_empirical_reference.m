function empirical=load_empirical_reference(filePath)

if ~isfile(filePath)
    error('Empirical reference file not found: %s',filePath);
end
phase=readtable(filePath,'Sheet','Fig3b_summary','TextType','string', ...
    'PreserveVariableNames',true);
loss=readtable(filePath,'Sheet','Fig3c_ratio_summary','TextType','string', ...
    'PreserveVariableNames',true);
late=readtable(filePath,'Sheet','Fig3d_summary','TextType','string', ...
    'PreserveVariableNames',true);
if width(phase)<8 || width(loss)<7 || width(late)<7
    error('Empirical reference workbook has an unexpected column layout.');
end

sourceMetrics={'tau80_over_tau50','tau90_over_tau80','tau95_over_tau90'};
values=nan(5,3);
for k=1:3
    row=find_one(phase,1,sourceMetrics{k},2,'All systems',filePath);
    values(k,:)=phase{row,[7 4 8]};
end
row=find_one(loss,1,'All systems',0,'',filePath);
values(4,:)=loss{row,[6 3 7]};
row=find_one(late,1,'All systems',0,'',filePath);
values(5,:)=late{row,[6 3 7]};

if any(~isfinite(values),'all') || ...
        any(values(:,1)>values(:,2) | values(:,2)>values(:,3))
    error('Empirical reference values must be finite and ordered Q05 <= median <= Q95.');
end
empirical.metric={'tau80_tau50','tau90_tau80','tau95_tau90','kappa','eta90'};
empirical.lower=values(:,1)';
empirical.median=values(:,2)';
empirical.upper=values(:,3)';
empirical.intervalLabel='pooled empirical 5th-95th percentile interval';
[~,name,extension]=fileparts(filePath);
empirical.sourceFile=[name extension];
end

function row=find_one(t,keyColumn,key,groupColumn,group,filePath)
mask=strcmpi(string(t{:,keyColumn}),key);
if groupColumn>0
    mask=mask & strcmpi(string(t{:,groupColumn}),group);
end
row=find(mask);
if numel(row)~=1
    error('Expected one %s / %s row in %s.',key,group,filePath);
end
end
