function ensure_output_dirs(cfg)
if ~exist(cfg.outputDir, 'dir')
    mkdir(cfg.outputDir);
end
if ~exist(cfg.cacheDir, 'dir')
    mkdir(cfg.cacheDir);
end
end
