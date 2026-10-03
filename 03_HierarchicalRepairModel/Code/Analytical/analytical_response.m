function out = analytical_response(limitName,t,depths,model)

if nargin<4 || isempty(model), model=struct; end
if ~isstruct(model) || ~isscalar(model)
    error('model must be a scalar structure.');
end
if ~isnumeric(t) || ~isreal(t)
    error('t must be a real numeric array.');
end
if ~isnumeric(depths) || isempty(depths) || any(~isfinite(depths(:))) || ...
        any(depths(:)<0) || any(depths(:)~=floor(depths(:)))
    error('depths must contain non-negative integers.');
end

depths=depths(:)';
tColumn=max(t(:),0);
meanDuration=get_option(model,'meanDuration',1);
if ~isscalar(meanDuration) || ~isfinite(meanDuration) || meanDuration<=0
    error('model.meanDuration must be a positive finite scalar.');
end
releaseSpec=get_option(model,'releaseSpec',struct('type','identity'));

nT=numel(tColumn);
nL=numel(depths);
connected=zeros(nT,nL);

switch lower(char(string(limitName)))
    case 'frontier'
        theta=tColumn./meanDuration;
        maxDepth=max(depths);
        allConnected=zeros(nT,maxDepth+1);
        term=exp(-theta);
        survival=term;
        allConnected(:,1)=1-survival;
        for l=1:maxDepth
            term=term.*theta/l;
            survival=survival+term;
            allConnected(:,l+1)=1-survival;
        end
        connected=allConnected(:,depths+1);
        task=connected;

    case 'parallel'
        distName=get_option(model,'distName','exponential');
        durationCV=get_option(model,'durationCV',1);
        F=duration_cdf(tColumn,distName,meanDuration,durationCV);
        task=repmat(F,1,nL);
        for k=1:nL
            connected(:,k)=F.^(depths(k)+1);
        end

    otherwise
        error('Unknown analytical repair limit: %s',limitName);
end

restored=zeros(nT,nL);
for k=1:nL
    restored(:,k)=service_release(connected(:,k),releaseSpec,depths(k)+1);
end

out=struct;
out.limitName=lower(char(string(limitName)));
out.time=tColumn;
out.depths=depths;
out.taskFraction=min(max(task,0),1);
out.connectedFraction=min(max(connected,0),1);
out.restoredFraction=min(max(restored,0),1);
out.systemRestoredFraction=[];

if isfield(model,'pi') && ~isempty(model.pi)
    piVec=model.pi(:);
    if numel(piVec)~=nL || any(~isfinite(piVec)) || any(piVec<0) || ...
            sum(piVec)<=0
        error('model.pi must be non-negative and match the supplied depths.');
    end
    piVec=piVec/sum(piVec);
    out.systemRestoredFraction=out.restoredFraction*piVec;
end
end

function value=get_option(s,name,defaultValue)
if isfield(s,name) && ~isempty(s.(name))
    value=s.(name);
else
    value=defaultValue;
end
end
