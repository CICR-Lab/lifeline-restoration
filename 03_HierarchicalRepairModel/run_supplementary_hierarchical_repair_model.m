function figures=run_supplementary_hierarchical_repair_model(varargin)
% Render Supplementary Figures S13-S18 from packaged or generated inputs.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'Code','Analytical'));
addpath(fullfile(root,'Code','Simulation'));
addpath(fullfile(root,'Code','Interdependent'));
addpath(fullfile(root,'Code','Common'));
[mode,dataSource,dataRoot,outputRoot] = parseFigureOptions(root,varargin{:});
dataRoot=resolve_figure_data_root(root,mode,dataSource,dataRoot);
figures.figureS13=plot_supplementary_figure_s13(root,mode,dataRoot,outputRoot);
figures.figureS14=plot_supplementary_figure_s14(root,mode,dataRoot,outputRoot);
figures.figureS15=plot_supplementary_figure_s15(root,mode,dataRoot,outputRoot);
figures.figureS16=plot_supplementary_figure_s16(root,mode,dataRoot,outputRoot);
figures.figureS17=plot_supplementary_figure_s17(root,mode,dataRoot,outputRoot);
figures.figureS18=plot_supplementary_figure_s18(root,mode,dataRoot,outputRoot);
fprintf('\nSupplementary Figures S13-S18 completed in:\n  %s\n', ...
    fullfile(outputRoot,'supplementary'));
end

function [mode,dataSource,dataRoot,outputRoot] = parseFigureOptions(root,varargin)
mode='full';
if ~isempty(varargin) && isTextScalar(varargin{1}) && ...
        ismember(lower(char(varargin{1})),{'quick','full'})
    mode=lower(char(varargin{1}));
    varargin=varargin(2:end);
end
p=inputParser;
p.addParameter('Mode',mode,@(x)isTextScalar(x));
p.addParameter('DataSource','packaged',@(x)isTextScalar(x));
p.addParameter('DataRoot','',@(x)isTextScalar(x));
p.addParameter('OutputDir','',@(x)isTextScalar(x));
p.parse(varargin{:});
mode=lower(char(p.Results.Mode));
dataSource=lower(char(p.Results.DataSource));
if ~ismember(mode,{'quick','full'}), error('Mode must be ''quick'' or ''full''.'); end
if ~ismember(dataSource,{'packaged','generated'})
    error('DataSource must be ''packaged'' or ''generated''.');
end
dataRoot=char(p.Results.DataRoot);
outputRoot=char(p.Results.OutputDir);
if isempty(outputRoot), outputRoot=resolve_figure_output_root(root,mode); end
end

function tf=isTextScalar(value)
tf=(ischar(value) && (isrow(value) || isempty(value))) || ...
    (isstring(value) && isscalar(value));
end
