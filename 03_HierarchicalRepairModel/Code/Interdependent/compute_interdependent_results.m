function output=compute_interdependent_results(root,mode)
modelDir=fullfile(root,'Code','Simulation');
cfg=config_interdependent_experiment(root,mode,modelDir);
ensure_output_dirs(cfg);
catalog=build_coupled_experiment_catalog(cfg);
results=cell(numel(catalog),1);
for i=1:numel(catalog)
    fprintf('  Coupled [%3d/%3d] %s\n',i,numel(catalog),catalog(i).id);
    results{i}=run_coupled_scenario_cached(catalog(i),cfg);
end
analytical=compute_analytical_limits(cfg);
coupling=compute_coupling_analytical(analytical,cfg);
sourceFile=fullfile(cfg.outputDir,'Interdependent_results.mat');
for i=1:numel(results)
    if strcmp(catalog(i).id,'core_gam1_gb0p5')
        results{i}=struct('metrics',results{i}.metrics,'example',results{i}.example);
    else
        results{i}=struct('metrics',results{i}.metrics);
    end
end
save(sourceFile,'catalog','results','coupling','cfg','-v7');
output=struct('sourceFile',sourceFile,'cfg',cfg);
end
