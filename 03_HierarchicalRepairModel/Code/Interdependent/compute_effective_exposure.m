function exposure = compute_effective_exposure(net,w,directGate)

N=net.N;
directGate=logical(directGate(:));
w=w(:);
if numel(directGate)~=N || numel(w)~=N
    error('directGate and w must have length net.N.');
end
connected=false(N,1);
for l=0:net.L
    idx=net.layerIndices{l+1};
    for q=1:numel(idx)
        i=idx(q);
        if directGate(i)
            connected(i)=false;
        else
            connected(i)=prereq_satisfied(net,i,connected);
        end
    end
end
exposed=~connected;
exposure.directGate=directGate;
exposure.exposed=exposed;
exposure.gBdir=sum(w(directGate));
exposure.gB=sum(w(exposed));
exposure.directGateCount=nnz(directGate);
if exposure.gBdir>0
    exposure.propagationMultiplier=exposure.gB/exposure.gBdir;
else
    exposure.propagationMultiplier=NaN;
end
end
