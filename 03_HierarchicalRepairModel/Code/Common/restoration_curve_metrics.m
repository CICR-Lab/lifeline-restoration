function metrics=restoration_curve_metrics(t,R,method,tolerance,milestones)
if nargin<3 || isempty(method), method='linear'; end
if nargin<4 || isempty(tolerance), tolerance=1e-10; end
if nargin<5 || isempty(milestones), milestones=[0.50 0.80 0.90 0.95]; end
t=t(:); R=R(:);
if numel(t)~=numel(R) || numel(t)<2 || abs(t(1))>tolerance
    error('RestorationMetric:InvalidOrigin', ...
        'The restoration curve must start at time zero.');
end
[t,~,group]=unique(t,'stable');
R=accumarray(group,R,[],@max);
R=min(max(cummax(R),0),1);
if numel(milestones)~=4 || abs(milestones(3)-0.90)>tolerance
    error('RestorationMetric:MilestoneContract', ...
        'Metrics require four milestones with 0.90 in the third position.');
end
tau=milestone_time(t,R,milestones,method,tolerance);
if any(~isfinite(tau)) || any(tau<=0) || ...
        any(diff(tau)<-tolerance*max(1,max(abs(tau))))
    error('RestorationMetric:InvalidMilestone', ...
        'Milestones must be finite, positive and nondecreasing.');
end
switch lower(method)
    case 'linear'
        lossArea=trapz(t,1-R);
        after=t>tau(3);
        lateLossArea=trapz([tau(3);t(after)],[0.10;1-R(after)]);
    case 'step'
        lossArea=sum((1-R(1:end-1)).*diff(t));
        idx90=find(R>=0.90-tolerance,1,'first');
        lateLossArea=sum((1-R(idx90:end-1)).*diff(t(idx90:end)));
    otherwise
        error('RestorationMetric:UnknownMethod', ...
            'method must be ''linear'' or ''step''.');
end
if ~isfinite(lossArea) || lossArea<=0
    error('RestorationMetric:InvalidLoss', ...
        'Cumulative restoration loss must be positive and finite.');
end
lateLossArea=min(max(lateLossArea,0),lossArea);
metrics=struct('tau50',tau(1),'tau80',tau(2),'tau90',tau(3), ...
    'tau95',tau(4),'tau80_tau50',tau(2)/tau(1), ...
    'tau90_tau80',tau(3)/tau(2),'tau95_tau90',tau(4)/tau(3), ...
    'lossArea',lossArea,'kappa',lossArea/tau(1), ...
    'lateLossArea',lateLossArea,'eta90',lateLossArea/lossArea);
end
