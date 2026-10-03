function file=plot_interdependent_figure4h(root,mode,dataRoot,outputRoot)
if nargin<3 || isempty(dataRoot), dataRoot=resolve_figure_data_root(root,mode); end
if nargin<4 || isempty(outputRoot), outputRoot=resolve_figure_output_root(root,mode); end
sourceFile=fullfile(dataRoot,'Interdependent', ...
    'Interdependent_results.mat');
if ~isfile(sourceFile)
    error('Missing Figure 4h source data. Run run_hierarchical_repair_analysis first.');
end
s=load(sourceFile,'catalog','results','coupling','cfg');
s.cfg.outputDir=fullfile(outputRoot,'main');
if ~isfolder(s.cfg.outputDir), mkdir(s.cfg.outputDir); end
plot_figure4h(s.catalog,s.results,s.coupling,s.cfg);
file=fullfile(s.cfg.outputDir,'Figure4h_coupled_network_contours.fig');
end
