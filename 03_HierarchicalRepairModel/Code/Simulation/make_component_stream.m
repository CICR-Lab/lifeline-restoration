function stream = make_component_stream(cfg, blockID, componentID)

substream = 1 + double(blockID) * 16 + double(componentID);
try
    stream = RandStream(cfg.randomGenerator, 'Seed', cfg.masterSeed);
catch ME
    error('RepairTask:RandomGeneratorUnavailable', ...
        'Required generator %s is unavailable: %s', ...
        cfg.randomGenerator, ME.message);
end
stream.Substream = substream;
end
