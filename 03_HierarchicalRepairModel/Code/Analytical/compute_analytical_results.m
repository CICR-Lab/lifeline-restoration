function output=compute_analytical_results(root,mode)
cfg=config_analytical(root,mode);
if ~isfolder(cfg.resultDir), mkdir(cfg.resultDir); end
theta=cfg.theta.grid;
identity=struct('type','identity');
parallelModel=struct('distName',cfg.parallel.durationFamily, ...
    'meanDuration',cfg.parallel.meanDuration, ...
    'durationCV',cfg.parallel.durationCV,'releaseSpec',identity);

%% Figure 4b: depth-specific normalized restoration trajectories
b.depths=0:3;
b.theta=theta;
b.frontier=analytical_response('frontier',theta,b.depths, ...
    struct('meanDuration',1,'releaseSpec',identity)).restoredFraction;
b.parallel=analytical_response('parallel',theta,b.depths, ...
    parallelModel).restoredFraction;
b.uFrontier=zeros(size(b.frontier)); b.uParallel=zeros(size(b.parallel));
for k=1:numel(b.depths)
    b.uFrontier(:,k)=theta./milestone_time(theta,b.frontier(:,k),.5,'linear');
    b.uParallel(:,k)=theta./milestone_time(theta,b.parallel(:,k),.5,'linear');
end

%% Figure 4c: restoration milestone ratios by dependency depth
c.depths=0:7;
f=analytical_response('frontier',theta,c.depths, ...
    struct('meanDuration',1,'releaseSpec',identity)).restoredFraction;
p=analytical_response('parallel',theta,c.depths,parallelModel).restoredFraction;
c.frontier=zeros(numel(c.depths),3); c.parallel=c.frontier;
for k=1:numel(c.depths)
    mf=restoration_curve_metrics(theta,f(:,k),'linear');
    mp=restoration_curve_metrics(theta,p(:,k),'linear');
    c.frontier(k,:)=[mf.tau80_tau50 mf.tau90_tau80 mf.tau95_tau90];
    c.parallel(k,:)=[mp.tau80_tau50 mp.tau90_tau80 mp.tau95_tau90];
end
[c.empirical,complete]=load_empirical_indicator_ranges(cfg.data.empiricalRanges, ...
    {'tau80_tau50','tau90_tau80','tau95_tau90'});
if ~complete, error('Incomplete pooled empirical milestone ranges.'); end

%% Figure 4d: restoration trajectories by layer service-share profile
d.profiles=service_profiles(3);
d.uGrid=linspace(0,10,801)';
L=numel(d.profiles(1).pi)-1;
lf=analytical_response('frontier',theta,0:L, ...
    struct('meanDuration',1,'releaseSpec',identity)).restoredFraction;
lp=analytical_response('parallel',theta,0:L,parallelModel).restoredFraction;
n=numel(d.profiles);
d.frontier=zeros(numel(d.uGrid),n); d.parallel=d.frontier;
d.metrics=repmat(struct('name','','pi',[],'frontier',[], ...
    'parallel',[]),n,1);
for k=1:n
    pi=d.profiles(k).pi(:);
    rf=lf*pi; rp=lp*pi;
    t50f=milestone_time(theta,rf,.5,'linear');
    t50p=milestone_time(theta,rp,.5,'linear');
    d.frontier(:,k)=interp1(theta/t50f,rf,d.uGrid,'pchip');
    d.parallel(:,k)=interp1(theta/t50p,rp,d.uGrid,'pchip');
    d.metrics(k)=struct('name',d.profiles(k).name, ...
        'pi',d.profiles(k).pi,'frontier',restoration_curve_metrics(theta,rf,'linear'), ...
        'parallel',restoration_curve_metrics(theta,rp,'linear'));
end

%% Supplementary Figures S13-S16
[s13,s14,s15,s16]=compute_analytical_supplementary(cfg);
sourceFile=fullfile(cfg.resultDir,'Analytical_results.mat');
save(sourceFile,'b','c','d','s13','s14','s15','s16','cfg','-v7.3');
output=struct('sourceFile',sourceFile,'cfg',cfg);
end
