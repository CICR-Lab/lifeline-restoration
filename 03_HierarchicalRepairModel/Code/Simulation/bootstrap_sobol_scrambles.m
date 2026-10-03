function boot = bootstrap_sobol_scrambles(data,nBootstrap,bootstrapSeed)

arguments
    data struct
    nBootstrap (1,1) double {mustBeInteger,mustBePositive} = 2000
    bootstrapSeed (1,1) double {mustBeInteger,mustBeNonnegative} = 260823
end
if ~isfield(data,'scrambleID') || numel(data.scrambleID)~=size(data.A,1)
    error('AllFactorGSA:MissingScrambleID', ...
        'data.scrambleID must identify every matched Sobol pair.');
end
scrambleID = data.scrambleID(:);
blocks = unique(scrambleID,'stable');
if numel(blocks)<2
    error('AllFactorGSA:TooFewScrambles', ...
        'At least two independent scrambles are required for resampling.');
end

stream = RandStream('mt19937ar','Seed',bootstrapSeed);
template = compute_sobol_indices(data);
nTarget = template.nTargets;
nOutcome = template.nOutcomes;
first = nan(nBootstrap,nTarget,nOutcome);
total = nan(nBootstrap,nTarget,nOutcome);
parameterShare = nan(nBootstrap,nOutcome);
stochasticShare = nan(nBootstrap,nOutcome);

stats=block_sufficient_statistics(data,scrambleID,blocks);
for b = 1:nBootstrap
    selectedIndex=randi(stream,numel(blocks),[numel(blocks),1]);
    multiplicity=accumarray(selectedIndex,1,[numel(blocks),1]);
    e=estimate_from_block_statistics(stats,multiplicity,template);
    first(b,:,:) = e.firstOrderRaw;
    total(b,:,:) = e.totalOrderRaw;
    parameterShare(b,:) = e.parameterShare;
    stochasticShare(b,:) = e.stochasticShare;
end

boot.nBootstrap = nBootstrap;
boot.bootstrapSeed = bootstrapSeed;
boot.bootstrapUnit = 'complete_scramble_block';
boot.firstOrderP025 = percentile_first_dimension(first,0.025);
boot.firstOrderP975 = percentile_first_dimension(first,0.975);
boot.totalOrderP025 = percentile_first_dimension(total,0.025);
boot.totalOrderP975 = percentile_first_dimension(total,0.975);
boot.parameterShareP025 = percentile_first_dimension(parameterShare,0.025);
boot.parameterShareP975 = percentile_first_dimension(parameterShare,0.975);
boot.stochasticShareP025 = percentile_first_dimension(stochasticShare,0.025);
boot.stochasticShareP975 = percentile_first_dimension(stochasticShare,0.975);
end

function stats=block_sufficient_statistics(data,scrambleID,blocks)
A=double(data.A); B=double(data.B); H=double(data.hybrid);
if ismatrix(A), A=reshape(A,size(A,1),size(A,2),1); end
if ismatrix(B), B=reshape(B,size(B,1),size(B,2),1); end
if ndims(H)==3, H=reshape(H,size(H,1),size(H,2),size(H,3),1); end
if any(~isfinite(A),'all') || any(~isfinite(B),'all') || ...
        any(~isfinite(H),'all')
    error('AllFactorGSA:NonfiniteBootstrapInput', ...
        'Fast scramble bootstrap requires finite A/B/hybrid values.');
end
nBlock=numel(blocks); nTarget=size(H,2); nOutcome=size(A,3);
stats.count=zeros(nBlock,1);
stats.x1=zeros(nBlock,nOutcome); stats.x2=zeros(nBlock,nOutcome);
stats.x1x2=zeros(nBlock,nOutcome);
stats.b1=zeros(nBlock,nOutcome); stats.b2=zeros(nBlock,nOutcome);
stats.h1=zeros(nBlock,nTarget,nOutcome); stats.h2=stats.h1;
stats.b1h2=zeros(nBlock,nTarget,nOutcome); stats.b2h1=stats.b1h2;
stats.totalProduct=zeros(nBlock,nTarget,nOutcome);
stats.stochastic=zeros(nBlock,nOutcome);
for s=1:nBlock
    idx=scrambleID==blocks(s); stats.count(s)=nnz(idx);
    for o=1:nOutcome
        a1=A(idx,1,o); a2=A(idx,2,o);
        b1=B(idx,1,o); b2=B(idx,2,o);
        stats.x1(s,o)=sum(a1)+sum(b1);
        stats.x2(s,o)=sum(a2)+sum(b2);
        stats.x1x2(s,o)=sum(a1.*a2)+sum(b1.*b2);
        stats.b1(s,o)=sum(b1); stats.b2(s,o)=sum(b2);
        for k=1:nTarget
            h1=H(idx,k,1,o); h2=H(idx,k,2,o);
            stats.h1(s,k,o)=sum(h1); stats.h2(s,k,o)=sum(h2);
            stats.b1h2(s,k,o)=sum(b1.*h2);
            stats.b2h1(s,k,o)=sum(b2.*h1);
            stats.totalProduct(s,k,o)=sum((a1-h1).*(a2-h2));
        end
    end
    stats.stochastic(s,:)=sum(data.rowStochasticVariance(idx,:),1);
end
end

function e=estimate_from_block_statistics(s,w,template)
w=double(w(:)); n=sum(w.*s.count);
if n<=0, error('AllFactorGSA:EmptyBootstrap','Empty scramble resample.'); end
weighted=@(x) sum(x.*reshape(w,[numel(w),ones(1,ndims(x)-1)]),1);
m1=reshape(weighted(s.x1)./(2*n),1,[]);
m2=reshape(weighted(s.x2)./(2*n),1,[]);
v=reshape(weighted(s.x1x2)./(2*n),1,[])-m1.*m2;
if any(~isfinite(v) | v<=0)
    error('AllFactorGSA:BootstrapNonpositiveParameterVariance', ...
        ['A complete-scramble bootstrap replicate has nonpositive ' ...
         'cross-batch parameter variance.']);
end
nTarget=size(s.h1,2); nOutcome=numel(v);
S1=nan(nTarget,nOutcome); ST=S1;
for o=1:nOutcome
    for k=1:nTarget
        eb1h2=weighted(s.b1h2(:,k,o))./n;
        eb2h1=weighted(s.b2h1(:,k,o))./n;
        eb1=weighted(s.b1(:,o))./n; eb2=weighted(s.b2(:,o))./n;
        eh1=weighted(s.h1(:,k,o))./n; eh2=weighted(s.h2(:,k,o))./n;
        first=0.5*((eb1h2-m1(o)*eh2-m2(o)*eb1+m1(o)*m2(o))+ ...
            (eb2h1-m2(o)*eh1-m1(o)*eb2+m1(o)*m2(o)));
        total=0.5*weighted(s.totalProduct(:,k,o))./n;
        S1(k,o)=first./v(o); ST(k,o)=total./v(o);
    end
end
sv=reshape(weighted(s.stochastic)./n,1,[]);
totalVariance=max(v,0)+sv;
e=template;
e.parameterVariance=v;
e.firstOrderRaw=S1; e.totalOrderRaw=ST;
e.parameterShare=max(v,0)./totalVariance;
e.stochasticShare=sv./totalVariance;
end

function q = percentile_first_dimension(x,p)
sz = size(x);
flat = reshape(x,sz(1),[]);
qFlat = nan(1,size(flat,2));
for j=1:size(flat,2)
    values=sort(flat(isfinite(flat(:,j)),j));
    if isempty(values), continue; end
    position=1+(numel(values)-1)*p;
    lo=floor(position); hi=ceil(position);
    if lo==hi
        qFlat(j)=values(lo);
    else
        qFlat(j)=values(lo)+(position-lo)*(values(hi)-values(lo));
    end
end
q=reshape(qFlat,[1,sz(2:end)]);
q=squeeze(q);
end
