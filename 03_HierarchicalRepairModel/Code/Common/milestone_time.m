function tau=milestone_time(t,R,p,method,tolerance)
if nargin<4 || isempty(method), method='linear'; end
if nargin<5 || isempty(tolerance), tolerance=1e-10; end
t=t(:); R=R(:); outputSize=size(p); p=p(:);
if numel(t)~=numel(R) || numel(t)<2 || any(~isfinite(t)) || ...
        any(~isfinite(R)) || any(diff(t)<0) || ...
        any(diff(R)<-tolerance) || any(p<0 | p>1)
    error('RestorationMetric:InvalidCurve', ...
        'The curve must be finite, time-ordered and monotone.');
end
[t,~,group]=unique(t,'stable');
R=accumarray(group,R,[],@max);
R=min(max(cummax(R),0),1);
tau=nan(size(p));
for k=1:numel(p)
    idx=find(R>=p(k)-tolerance,1,'first');
    if isempty(idx), continue; end
    if strcmpi(method,'step') || idx==1 || ...
            R(idx)-R(idx-1)<=tolerance
        tau(k)=t(idx);
    elseif strcmpi(method,'linear')
        fraction=(p(k)-R(idx-1))/(R(idx)-R(idx-1));
        tau(k)=t(idx-1)+min(max(fraction,0),1)*(t(idx)-t(idx-1));
    else
        error('RestorationMetric:UnknownMethod', ...
            'method must be ''linear'' or ''step''.');
    end
end
tau=reshape(tau,outputSize);
end
