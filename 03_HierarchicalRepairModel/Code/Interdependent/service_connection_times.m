function tConn = service_connection_times(net, repairFinishTime, directGate, supportTime)

N = net.N;
repairFinishTime = repairFinishTime(:);
directGate = logical(directGate(:));
supportTime = supportTime(:);
if numel(repairFinishTime)~=N || numel(directGate)~=N || numel(supportTime)~=N
    error('Input vectors must have length net.N.');
end

tConn = nan(N,1);
for l = 0:net.L
    idx = net.layerIndices{l+1};
    for q = 1:numel(idx)
        i = idx(q);
        p = net.parents{i};
        if isempty(p)
            parentTime = 0;
        elseif net.parentLogic(i) == "OR"
            parentTime = min(tConn(p));
        else
            parentTime = max(tConn(p));
        end
        gateTime = 0;
        if directGate(i)
            gateTime = supportTime(i);
        end
        tConn(i) = max([repairFinishTime(i), parentTime, gateTime]);
    end
end
end
