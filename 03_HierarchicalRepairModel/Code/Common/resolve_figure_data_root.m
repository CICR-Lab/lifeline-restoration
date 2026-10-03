function dataRoot = resolve_figure_data_root(root, mode, dataSource, dataRootOverride)
% Resolve one explicit figure-data source without automatic fallback.
mode = lower(char(mode));
if nargin < 3 || isempty(dataSource)
    dataSource = 'packaged';
end
dataSource = lower(char(dataSource));

if nargin >= 4 && ~isempty(dataRootOverride)
    dataRoot = char(dataRootOverride);
elseif strcmp(dataSource, 'packaged')
    dataRoot = fullfile(root, 'figure_data', mode);
elseif strcmp(dataSource, 'generated')
    dataRoot = fullfile(root, 'outputs', 'figure_data', mode);
else
    error('HierarchicalRepair:InvalidDataSource', ...
        'DataSource must be ''packaged'' or ''generated''.' );
end

if ~isfolder(dataRoot)
    error('HierarchicalRepair:MissingFigureData', ...
        'Figure-data directory is missing for DataSource=%s, mode=%s: %s', ...
        dataSource, mode, dataRoot);
end
end
