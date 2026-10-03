function out=run_coupled_scenario_cached(scenario,cfg)

cacheTag=sprintf('coupled_%s_%s_n%d',scenario.id,char(cfg.mode),scenario.nRep);
cacheFile=fullfile(cfg.cacheDir,[cacheTag '.mat']);
signature=scenario_signature(scenario,cfg.algorithmVersion);
if isfile(cacheFile) && ~cfg.forceRerun
    saved=load(cacheFile,'out','cacheMeta');
    if isfield(saved,'out') && isfield(saved,'cacheMeta') && ...
            strcmp(saved.cacheMeta.scenarioSignature,signature)
        out=saved.out;
        return
    end
end
out=run_coupled_scenario_ensemble(scenario,cfg);
cacheMeta=struct('algorithmVersion',cfg.algorithmVersion, ...
    'scenarioSignature',signature,'createdAt',datestr(now,30));
save(cacheFile,'out','cacheMeta','-v7.3');
if cfg.saveRawMetrics
    writetable(out.metrics,fullfile(cfg.outputDir,[cacheTag '_raw_metrics.csv']));
end
end
