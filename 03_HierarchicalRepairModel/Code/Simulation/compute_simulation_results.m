function output=compute_simulation_results(root,mode)
outputRoot=fullfile(root,'outputs',mode,'Simulation');
if ~isfolder(outputRoot), mkdir(outputRoot); end
cfg=config_all_factor_experiment(mode,outputRoot);
design=build_all_factor_design(cfg);
cfg.experimentID=compute_sha256(jsonencode(table2struct(design.cases)));
cfg.paths=all_factor_paths(cfg);
plan=run_all_factor_checkpoints(cfg,design);
empirical=load_empirical_reference(fullfile(fileparts(root), ...
    '02_EmergentSimplicity','figure_data','main','SourceData_Figure3.xlsx'));
summary=aggregate_all_factor_results(cfg,design,empirical);
pub=prepare_publication_data(summary);
sourceDir=fullfile(outputRoot,'source_data');
if ~isfolder(sourceDir), mkdir(sourceDir); end
writetable(pub.factorSobol,fullfile(sourceDir, ...
    'Fig4e_global_sensitivity.csv'));
writetable(pub.taskService,fullfile(sourceDir, ...
    'Fig4f_task_service.csv'));
writetable(pub.capacityQ,fullfile(sourceDir, ...
    'Fig4g_capacity_priority.csv'));
writetable(pub.consistency(1:6,:),fullfile(sourceDir, ...
    'FigS17a_consistency.csv'));
writetable(pub.consistency(7,:),fullfile(sourceDir, ...
    'FigS17a_allfive.csv'));
writetable(pub.modelDataAudit,fullfile(sourceDir, ...
    'FigS17b_model_data_audit.csv'));
writetable(pub.variance,fullfile(sourceDir, ...
    'FigS17c_variance_structure.csv'));
writetable(pub.groupSobol,fullfile(sourceDir, ...
    'FigS17d_group_sobol.csv'));
tableDir=fullfile(root,'outputs',mode,'SupplementaryTables');
if ~isfolder(tableDir), mkdir(tableDir); end
tableFile=fullfile(tableDir,'Supplementary_Table_S8.xlsx');
if isfile(tableFile), delete(tableFile); end
writetable(pub.tableS8,tableFile,'Sheet','Inputs');
output=struct('cfg',cfg,'plan',plan,'sourceDir',sourceDir,'tableFile',tableFile);
end
