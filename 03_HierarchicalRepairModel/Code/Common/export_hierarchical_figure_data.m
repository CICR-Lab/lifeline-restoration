function export_hierarchical_figure_data(root, mode)
% Copy direct plotting inputs into the regenerable generated figure-data area.
sourceRoot = fullfile(root, 'outputs', mode);
targetRoot = fullfile(root, 'outputs', 'figure_data', mode);
required = {
    fullfile('Analytical', 'Analytical_results.mat')
    fullfile('Interdependent', 'Interdependent_results.mat')
    fullfile('Simulation', 'source_data', 'Fig4e_global_sensitivity.csv')
    fullfile('Simulation', 'source_data', 'Fig4f_task_service.csv')
    fullfile('Simulation', 'source_data', 'Fig4g_capacity_priority.csv')
    fullfile('Simulation', 'source_data', 'FigS17a_allfive.csv')
    fullfile('Simulation', 'source_data', 'FigS17a_consistency.csv')
    fullfile('Simulation', 'source_data', 'FigS17b_model_data_audit.csv')
    fullfile('Simulation', 'source_data', 'FigS17c_variance_structure.csv')
    fullfile('Simulation', 'source_data', 'FigS17d_group_sobol.csv')
    fullfile('SupplementaryTables', 'Supplementary_Table_S8.xlsx')};
for i = 1:numel(required)
    sourceFile = fullfile(sourceRoot, required{i});
    if ~isfile(sourceFile)
        error('HierarchicalRepair:MissingFigureInput', ...
            'Expected figure input was not generated: %s', sourceFile);
    end
    targetFile = fullfile(targetRoot, required{i});
    targetDir = fileparts(targetFile);
    if ~isfolder(targetDir), mkdir(targetDir); end
    copyfile(sourceFile, targetFile, 'f');
end
end
