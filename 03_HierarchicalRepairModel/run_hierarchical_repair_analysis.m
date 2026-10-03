function output=run_hierarchical_repair_analysis(mode)
if nargin<1 || isempty(mode), mode='full'; end
mode=validate_mode(mode);
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'Code','Analytical'));
addpath(fullfile(root,'Code','Simulation'));
addpath(fullfile(root,'Code','Interdependent'));
addpath(fullfile(root,'Code','Common'));

output.analytical=compute_analytical_results(root,mode);
output.simulation=compute_simulation_results(root,mode);
output.interdependent=compute_interdependent_results(root,mode);
export_hierarchical_figure_data(root,mode);
fprintf('\nMain and supplementary calculations completed in:\n  %s\n', ...
    fullfile(root,'outputs',mode));
end

function mode=validate_mode(mode)
mode=lower(char(mode));
if ~ismember(mode,{'quick','full'})
    error('Mode must be ''quick'' or ''full''.');
end
end
