% Update path to your .mat files folder
folder = '/Users/farhan/Downloads/Farhan_PPP';
files = dir(fullfile(folder, '*.mat'));
filename = 'output_params.xlsx';

for k = 1:length(files)
    matfile = fullfile(folder, files(k).name);
    data = load(matfile);
    params = data.params;
    
    % Sheet name from filename (truncate to 31 chars)
    sheet = strrep(files(k).name(1:end-4), ' ', '_');
    sheet = sheet(1:min(31, length(sheet)));
    
    % Handle struct or other types
    if isstruct(params)
        % Convert scalar/array of structs to table (fields as columns)
        param_table = struct2table(params);
    else
        % For numeric/cell, etc., convert to table
        param_table = array2table(params);
    end
    
    % Write table to specific sheet
    writetable(param_table, filename, 'Sheet', sheet, 'WriteMode', 'append');
end
