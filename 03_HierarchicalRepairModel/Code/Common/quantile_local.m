function q=quantile_local(x,p)
x=sort(x(isfinite(x)));
if isempty(x), q=nan(size(p)); return; end
q=zeros(size(p));
for k=1:numel(p)
    position=1+(numel(x)-1)*p(k);
    lo=floor(position); hi=ceil(position);
    if lo==hi, q(k)=x(lo);
    else, q(k)=x(lo)+(position-lo)*(x(hi)-x(lo)); end
end
end
