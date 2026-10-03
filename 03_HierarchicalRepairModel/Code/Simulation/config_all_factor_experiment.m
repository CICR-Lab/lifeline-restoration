function cfg=config_all_factor_experiment(mode,outputDirectory)

if nargin<1||isempty(mode),mode='quick';end
if nargin<2||isempty(outputDirectory),outputDirectory=fullfile(pwd,'outputs',mode);end
mode=lower(char(mode));
if ~ismember(mode,{'quick','full'}),error('Mode must be quick or full.');end
cfg.modelVersion='finite-network-v3';
cfg.mode=mode;
cfg.masterSeed=260818;
cfg.bootstrapSeed=260824;
cfg.matlabRelease=version('-release');
if verLessThan('matlab','9.7'),error('MATLAB R2019b or later is required.');end
cfg.randomGenerator='mrg32k3a';
cfg.sobolSkip=1024;
[factorInfo,groupInfo]=all_factor_metadata();
cfg.nFactors=height(factorInfo);
cfg.factorNames=factorInfo.factorName';
cfg.factorLabels=factorInfo.factorLabel';
cfg.groupNames=groupInfo.groupName';
cfg.groupFactorIDs=cell(1,height(groupInfo));
for g=1:height(groupInfo)
    cfg.groupFactorIDs{g}=factorInfo.factorID(factorInfo.groupID==g)';
end

cfg.nTasksBounds=[100,2000];
cfg.maximumDepth=[2,4,6,8];
cfg.taskLayerProfiles={'shallow','uniform','deep'};
cfg.taskCountConcentrationBounds=[0.5,10];
cfg.dependencyLogic={'single','and2','or2'};
cfg.parentConcentrationBounds=[0.5,10];
cfg.serviceLayerProfiles={ ...
    'shallow','uniform','intermediate','deep','deep_tail'};
cfg.durationDepthRatioBounds=[0.5,2.0];
cfg.cvDurationBounds=[0,2];
cfg.cvWeightBounds=[0,2];
cfg.rhoDurationWeightBounds=[-0.6,0.6];
cfg.crewFractionBounds=[0.01,1.0];
cfg.qFBounds=[0,1];
cfg.capacityAnchors=[0.01,0.03,0.10,0.30,1.00];
cfg.qFAnchors=[0,0.5,1];

switch cfg.mode
    case 'quick'
        cfg.nScrambles=1;
        cfg.nBasePerScramble=64;
        cfg.nSeeds=2;
        cfg.nResponseScrambles=1;
        cfg.nBackgroundsPerResponseScramble=4;
        cfg.useParallel=false;
        cfg.maxWorkers=1;
        cfg.nBootstrap=0;
    case 'full'
        cfg.nScrambles=8;
        cfg.nBasePerScramble=512;
        cfg.nSeeds=50;
        cfg.nResponseScrambles=8;
        cfg.nBackgroundsPerResponseScramble=32;
        cfg.useParallel=true;
        cfg.maxWorkers=8;
        cfg.nBootstrap=2000;
end
cfg.nResponseBackgrounds=cfg.nResponseScrambles* ...
    cfg.nBackgroundsPerResponseScramble;
cfg.seedBatchSize=cfg.nSeeds/2;
cfg.nBase=cfg.nScrambles*cfg.nBasePerScramble;
cfg.nSobolMatrices=2+cfg.nFactors+numel(cfg.groupNames);
cfg.nResponseCells=numel(cfg.taskLayerProfiles)* ...
    numel(cfg.serviceLayerProfiles);
if numel(cfg.capacityAnchors)*numel(cfg.qFAnchors)~=cfg.nResponseCells
    error('RepairTask:AllFactorResponseCellContract', ...
        'The two matched response modules must contain the same 15 cells.');
end
cfg.expectedSobolCases=cfg.nBase*cfg.nSobolMatrices;
cfg.expectedTaskServiceCases=cfg.nResponseBackgrounds*cfg.nResponseCells;
cfg.expectedCapacityQCases=cfg.nResponseBackgrounds*cfg.nResponseCells;
cfg.expectedCases=cfg.expectedSobolCases+ ...
    cfg.expectedTaskServiceCases+cfg.expectedCapacityQCases;
cfg.expectedTrajectories=cfg.expectedCases*cfg.nSeeds;
cfg.basePointsPerCheckpoint=64;
cfg.nRandomBlocks=cfg.nBase+2*cfg.nResponseBackgrounds;

cfg.milestones=[0.50,0.80,0.90,0.95];
cfg.simultaneousTolerance=1e-11;
cfg.numericTolerance=1e-10;
cfg.minimumTaskDuration=1e-8;
cfg.minTasksPerLayer=2;
cfg.durationFamily='gamma';
cfg.gammaSampler='inverse_cdf_gammaincinv';
cfg.serviceReleaseRule='task_weight_on_reconnection';
cfg.durationLayerNormalization='exact_layer_mean_affine_floor';

cfg.outputRoot=char(outputDirectory);
end
