function outputRoot = resolve_figure_output_root(root, mode, outputRootOverride)
% Resolve the rendered-figure directory without mixing quick and full runs.
if nargin >= 3 && ~isempty(outputRootOverride)
    outputRoot = char(outputRootOverride);
    return;
end
mode = lower(char(mode));
switch mode
    case 'full'
        outputRoot = fullfile(root, 'outputs', 'figures');
    case 'quick'
        outputRoot = fullfile(root, 'outputs', 'validation', 'quick_figures');
    otherwise
        error('Mode must be ''quick'' or ''full''.');
end
end
