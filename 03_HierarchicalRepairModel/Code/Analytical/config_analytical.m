function cfg=config_analytical(root,mode)
if nargin<2 || isempty(mode), mode='full'; end
cfg.rootDir=root; cfg.mode=lower(char(mode)); cfg.seed=20260721;
cfg.resultDir=fullfile(root,'outputs',cfg.mode,'Analytical');
cfg.data.empiricalRanges=fullfile(fileparts(root), ...
    '02_EmergentSimplicity','figure_data','main','SourceData_Figure3.xlsx');
cfg.output.main=fullfile(root,'outputs','figures','main');
cfg.output.supp=fullfile(root,'outputs','figures','supplementary');
cfg.figure.fontName='Arial'; cfg.figure.fontSize=9; cfg.figure.lineWidth=1.6;
cfg.figure.mainPanelSizeCm=[7.4 7.2]; cfg.figure.widthCm=20;
cfg.figure.analytical01HeightCm=16.5; cfg.figure.analytical02HeightCm=30;
cfg.figure.analytical03HeightCm=14.5; cfg.figure.dpi=300;
cfg.colors.blue=[.1059 .4706 .7137]; cfg.colors.orange=[1 .498 .0549];
cfg.colors.green=[.1725 .6275 .1725]; cfg.colors.red=[.8392 .1529 .1569];
cfg.colors.purple=[.5804 .4039 .7412]; cfg.colors.gray=[.45 .45 .45];
cfg.colors.classMap=[.31 .68 .39;.98 .76 .23;.85 .23 .23];
cfg.theta.max=35;
if strcmp(cfg.mode,'full'), cfg.theta.n=14001; else, cfg.theta.n=5001; end
cfg.theta.grid=linspace(0,cfg.theta.max,cfg.theta.n)';
cfg.parallel=struct('durationFamily','exponential','meanDuration',1,'durationCV',1);
cfg.solver.relTol=1e-8; cfg.solver.absTol=1e-10;
cfg.solver.nEval=cfg.theta.n; cfg.solver.tMax=cfg.theta.max;
cfg.sensitivity.depthValues=1:7;
if strcmp(cfg.mode,'full')
    cfg.sensitivity.depthNormGrid=linspace(0,1,101);
else
    cfg.sensitivity.depthNormGrid=linspace(0,1,51);
end
cfg.sensitivity.indicatorNames={'tau80_tau50','tau90_tau80', ...
    'tau95_tau90','kappa','eta90'};
cfg.sensitivity.indicatorLabels={'\tau_{80}/\tau_{50}', ...
    '\tau_{90}/\tau_{80}','\tau_{95}/\tau_{90}','\kappa','\eta_{90}'};
cfg.sensitivity.repairLimits={'frontier','parallel'};
cfg.sensitivity.repairLimitLabels={'Frontier repair','Parallel repair'};
cfg.sensitivity.parallel=cfg.parallel;
cfg.sensitivity.rateNames={'Uniform','Shallow fast','Deep fast', ...
    'Middle bottleneck','Alternating'};
cfg.sensitivity.rateProfiles=[1 1 1 1;1.4 1.2 .8 .6;.6 .8 1.2 1.4; ...
    1.4 .6 .6 1.4;1.3 .7 1.3 .7];
cfg.sensitivity.accessNames={'Single prerequisite','2 mandatory', ...
    '3 mandatory','2 alternatives','3 alternatives'};
cfg.sensitivity.accessTypes={'single','mandatory','mandatory','alternative','alternative'};
cfg.sensitivity.accessM=[1 2 3 2 3];
cfg.sensitivity.releaseDurationNames={'Power q=0.75','Identity, exponential', ...
    'Power q=1.25','Threshold h=0.70','Threshold h=0.85', ...
    'Gamma CV=1.50','Lognormal CV=1.50'};
cfg.sensitivity.releaseDurationSpecs={release('power',.75,.70,20), ...
    release('identity',1,.70,20),release('power',1.25,.70,20), ...
    release('threshold',1,.70,20),release('threshold',1,.85,20), ...
    release('identity',1,.70,20),release('identity',1,.70,20)};
cfg.sensitivity.releaseDurationFamilies={'exponential','exponential', ...
    'exponential','exponential','exponential','gamma','lognormal'};
cfg.sensitivity.releaseDurationCV=[1 1 1 1 1 1.5 1.5];
cfg.sensitivity.releaseDurationTimeMax=80;
cfg.fit.uMax=4; cfg.fit.nStarts=25; cfg.fit.maxIter=2500;
end

function value=release(type,q,h,s)
value=struct('type',type,'q',q,'h',h,'s',s);
end
