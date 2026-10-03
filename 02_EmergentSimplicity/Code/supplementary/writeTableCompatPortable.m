function writeTableCompatPortable(T, filename, sheetName)

if nargin < 3 || strlength(string(sheetName)) == 0
    sheetName = "Sheet1";
end

filename = string(filename);
[folderPath, ~, extension] = fileparts(filename);

if strlength(folderPath) > 0 && exist(folderPath, "dir") ~= 7
    mkdir(folderPath);
end

C = tableToCellWithHeader(T);
C = sanitizeCell(C);

if strcmpi(extension, ".csv")
    writeCellCsv(C, filename);
else
    writeCellExcel(C, filename, sheetName);
end

end

function C = tableToCellWithHeader(T)

variables = T.Properties.VariableNames;

C = cell(height(T) + 1, width(T));
C(1, :) = variables;

for c = 1:width(T)

    column = T.(variables{c});

    for r = 1:height(T)

        if iscell(column)
            value = column{r};
        else
            value = column(r, :);
        end

        C{r + 1, c} = excelSafeValue(value);
    end
end

end

function writeCellExcel(C, filename, sheetName)

filename = char(filename);
sheetName = char(sheetName);

try

    writecell(C, filename, "Sheet", sheetName, "Range", "A1");

catch firstError

    try

        xlswrite(filename, C, sheetName, "A1");

    catch secondError

        error( ...
            "Excel output could not be written:\n%s\n" + ...
            "writecell: %s\n" + ...
            "xlswrite: %s", ...
            filename, ...
            firstError.message, ...
            secondError.message);
    end
end

end

function writeCellCsv(C, filename)

fid = fopen(char(filename), "w", "n", "UTF-8");

if fid < 0
    error("Cannot write CSV:\n%s", filename);
end

cleanupObject = onCleanup(@() fclose(fid));

for r = 1:size(C, 1)

    values = strings(1, size(C, 2));

    for c = 1:size(C, 2)
        values(c) = csvEscape(C{r, c});
    end

    fprintf(fid, "%s\n", strjoin(values, ","));
end

end

function C = sanitizeCell(C)

for r = 1:size(C, 1)
    for c = 1:size(C, 2)
        C{r, c} = excelSafeValue(C{r, c});
    end
end

end

function value = excelSafeValue(value)

if isempty(value)

    value = [];

elseif iscell(value)

    value = excelSafeValue(value{1});

elseif isstring(value)

    value(ismissing(value)) = "";
    value = char(strjoin(value(:).', " "));

elseif ischar(value)


elseif isnumeric(value)

    if ~isscalar(value)
        value = char(strjoin(string(value(:).'), " "));
    elseif ~isfinite(value)
        value = [];
    end

elseif islogical(value)

    value = double(value);

else

    try

        s = string(value);
        s(ismissing(s)) = "";
        value = char(strjoin(s(:).', " "));
    catch
        value = [];
    end
end

end

function s = csvEscape(value)

value = excelSafeValue(value);

if isempty(value)

    s = "";

elseif isnumeric(value)

    s = string(sprintf("%.15g", value));

else

    s = string(value);
end

s = replace(s, """", """""");

if any(contains(s, [",", newline, """"]))
    s = """" + s + """";
end

end
