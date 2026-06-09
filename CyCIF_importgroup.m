% Define the directory where your files are located
inputDir = '/Volumes/HITS/lsp-analysis/cycif-production/73-CRC1_reconstruction_Tcell_KQ/e13_Tcell_activation_abval/csvs';
% Loop through each file

% Get a list of all files in the directory with the specified format (if
% from server)
files = dir(fullfile(inputDir, '*', 'quantification', '*--unmicst_cellRing.csv'));

% Initialize a cell array to hold the data tables
dataStruct1 = struct();

% Loop through each file (for LSP*_ file names)
for i = 1:length(files)
    % Get the full filename
    fullFileName = fullfile(files(i).folder, files(i).name);

    % Extract the part of the filename before the first underscore/period
    underscoreIndex = strfind(files(i).name, '_');
    if isempty(underscoreIndex)
       disp(['Skipping file (no underscore found): ', files(i).name]);
       continue; % Skip this file if there's no underscore
    end
    newFileName = files(i).name(1:underscoreIndex(1)-1);
    disp(newFileName);

    % Call your custom import function (0.325 orion, 0.65 cycif)
    dataTable = CycIF_importMcMicro(fullFileName, 0.65);

    % Remove zeros and convert the table values to uint16
    dataTable = CycIF_removezero(dataTable); % Remove zeros

    % Loop through each column in the table and convert to uint16
    for col = 1:width(dataTable)
       if isnumeric(dataTable{:, col})
            dataTable{:,2:end}=uint16(dataTable{:,2:end});
       else
            warning('Column %d is not numeric and was not converted to uint16.', col);
       end
    end 
    % Assign the table to a variable named after newFileName
    % assignin('base', newFileName, dataTable);

    % Store the table in the ith cell of the cell array
    dataStruct1.(newFileName) = dataTable;
end

disp('Processing completed.');