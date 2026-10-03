function [t,R] = curve_from_connection_times(tConn,w)
tConn = tConn(:); w = w(:);
if numel(tConn)~=numel(w), error('tConn and w must have equal length.'); end
if any(~isfinite(tConn)), error('All connection times must be finite.'); end
[ts,ord] = sort(tConn);
ws = w(ord);
[tu,~,ic] = unique(ts);
wu = accumarray(ic,ws,[],@sum);
R = cumsum(wu);
t = tu;
t = [0;t]; R = [0;R];
R = min(max(cummax(R),0),1);
R(end)=1;
end
