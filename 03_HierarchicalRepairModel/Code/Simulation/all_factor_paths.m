function paths = all_factor_paths(cfg)

paths.outputRoot=cfg.outputRoot;
paths.checkpointDir=fullfile(paths.outputRoot,'checkpoints');
paths.failureDir=fullfile(paths.outputRoot,'failures');
if ispc && numel(paths.checkpointDir)>180
    warning('RepairTask:LongAllFactorOutputPath', ...
        ['The checkpoint path is long for Windows. Pass a short ' ...
         'OutputDirectory, for example C:\\rtgsa.']);
end
end
