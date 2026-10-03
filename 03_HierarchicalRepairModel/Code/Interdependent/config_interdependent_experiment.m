function cfg=config_interdependent_experiment(rootDir,mode,modelDir)
cfg.rootDir=rootDir;
cfg.mode=lower(string(mode));
cfg.outputDir=fullfile(rootDir,'outputs',char(cfg.mode),'Interdependent');
cfg.cacheDir=cfg.outputDir;
cfg.baseSeed=731904;
cfg.algorithmVersion='interdependent-v6-shared-qf-scheduler';
cfg.forceRerun=false;
cfg.saveRawMetrics=false;
cfg.figureResolution=450;
cfg.analyticalTheta=linspace(0,40,8001);

cfg.coupled.NE=1000; cfg.coupled.NB=1000;
cfg.coupled.crewFractionE=0.1; cfg.coupled.crewFractionB=0.1;
cfg.coupled.CVdE=0.5; cfg.coupled.CVdB=0.5;
cfg.coupled.CVwE=0.5; cfg.coupled.CVwB=0.5;
cfg.coupled.rhoDWE=0; cfg.coupled.rhoDWB=0;
cfg.coupled.qFE=1; cfg.coupled.qFB=1;
cfg.coupled.dependencyModeE='single_parent';
cfg.coupled.dependencyModeB='single_parent';
cfg.coupled.qE=0.90;
cfg.coupled.gateLocation='mixed';
cfg.coupled.supportRepresentation='aggregate_threshold';
cfg.coupled.gammaGrid=[0.25 0.5 1 2 4 8];
cfg.coupled.gBGrid=[0.1 0.3 0.5 0.7];
cfg.coupled.gammaPlotMax=9;
cfg.coupled.gBTolerance=0.02;
cfg.coupled.mainContourLevels=[0.05 0.10 0.20 0.40];
cfg.coupled.rho90ContourLevels=[1 2 4 8];

switch cfg.mode
    case "quick"
        cfg.nRepCoupledCore=8;
        cfg.nRepCoupledSensitivity=6;
        cfg.coupled.gateSearchTrials=3;
    case "full"
        cfg.nRepCoupledCore=500;
        cfg.nRepCoupledSensitivity=200;
        cfg.coupled.gateSearchTrials=8;
    otherwise
        error('Unknown mode "%s". Use quick or full.',mode);
end

cfg.generator.modelDir=modelDir;
cfg.generator.maximumDepthE=3; cfg.generator.maximumDepthB=3;
cfg.generator.taskLayerProfileE='uniform'; cfg.generator.taskLayerProfileB='uniform';
cfg.generator.taskCountConcentrationE=10; cfg.generator.taskCountConcentrationB=10;
cfg.generator.parentConcentrationE=10; cfg.generator.parentConcentrationB=10;
cfg.generator.serviceLayerProfileE='intermediate';
cfg.generator.serviceLayerProfileB='intermediate';
cfg.generator.durationDepthRatioE=1; cfg.generator.durationDepthRatioB=1;
cfg.coupled.piE=layer_profile_weights( ...
    cfg.generator.serviceLayerProfileE,cfg.generator.maximumDepthE)';
cfg.coupled.piB=layer_profile_weights( ...
    cfg.generator.serviceLayerProfileB,cfg.generator.maximumDepthB)';
validate_service_profile(cfg.coupled.piE,cfg.generator.maximumDepthE,'E');
validate_service_profile(cfg.coupled.piB,cfg.generator.maximumDepthB,'B');

generatorCfg=config_all_factor_experiment(char(cfg.mode),cfg.outputDir);
generatorCfg.nSeeds=max(cfg.nRepCoupledCore,cfg.nRepCoupledSensitivity);
generatorCfg.nRandomBlocks=2;
cfg.generator.config=generatorCfg;
end

function validate_service_profile(pi,L,systemName)
if numel(pi)~=L+1 || any(~isfinite(pi)) || any(pi<=0) || ...
        abs(sum(pi)-1)>1e-12
    error('Invalid layer service shares for system %s.',systemName);
end
end
