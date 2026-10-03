function exposure = construct_gate_set_for_target(net,w,targetGB,location,tolerance,nTrials)

N=net.N;
w=w(:);
pool=gate_candidate_pool(net,location);
if isempty(pool), error('Gate candidate pool is empty for location %s.',location); end

bestDiff=inf;
bestExposure=[];
for tr=1:max(1,nTrials)
    ord=pool(randperm(numel(pool)));
    lo=0; hi=numel(ord);
    while lo<hi
        mid=floor((lo+hi)/2);
        if mid==0
            g=0;
        else
            gate=false(N,1); gate(ord(1:mid))=true;
            tmp=compute_effective_exposure(net,w,gate); g=tmp.gB;
        end
        if g>=targetGB
            hi=mid;
        else
            lo=mid+1;
        end
    end
    candidates=unique(max(0,min(numel(ord),[lo-1 lo lo+1])));
    for k=candidates
        gate=false(N,1);
        if k>0, gate(ord(1:k))=true; end
        tmp=compute_effective_exposure(net,w,gate);
        diff=abs(tmp.gB-targetGB);
        if diff<bestDiff
            bestDiff=diff;
            bestExposure=tmp;
        end
    end
    if bestDiff<=tolerance, break; end
end

if isempty(bestExposure)
    error('Failed to construct a gate set.');
end
bestExposure.targetGB=targetGB;
bestExposure.targetTolerance=tolerance;
bestExposure.targetMet=(bestDiff<=tolerance);
bestExposure.gateLocation=string(location);
exposure=bestExposure;
end

function pool=gate_candidate_pool(net,location)
switch lower(location)
    case 'shallow'
        pool=find(net.layer==0);
    case 'intermediate'
        pool=find(net.layer==1);
    case 'deep'
        pool=find(net.layer>=max(2,net.L-1));
    case 'mixed'
        pool=(1:net.N)';
    otherwise
        error('Unknown gate location: %s',location);
end
pool=pool(:);
end
