% MATLAB script to select an HDF5 file and display information about it

% Ask the user to select an HDF5 file
[fileName, filePath] = uigetfile('*.hdf5', 'Select an HDF5 file');

% Check if the user selected a file or canceled the operation
if isequal(fileName, 0)
    disp('No file selected. Exiting...');
    return;
end

% Full path to the selected file
fullFileName = fullfile(filePath, fileName);

% Display the selected file
disp(['Selected file: ', fullFileName]);

% Open the HDF5 file
info = h5info(fullFileName);

% Display general information about the HDF5 file
disp('General Information about the HDF5 file:');
disp(['File: ', fullFileName]);
disp(['Number of Groups: ', num2str(length(info.Groups))]);
disp(['Number of Datasets: ', num2str(length(info.Datasets))]);

% Display the groups and datasets in the file
disp('Groups and Datasets:');
for i = 1:length(info.Groups)
    disp(['Group: ', info.Groups(i).Name]);
end

for i = 1:length(info.Datasets)
    disp(['Dataset: ', info.Datasets(i).Name]);
    
    % Display dataset attributes (if any)
    if isfield(info.Datasets(i), 'Attributes') && ~isempty(info.Datasets(i).Attributes)
        disp('  Attributes:');
        for j = 1:length(info.Datasets(i).Attributes)
            attrName = info.Datasets(i).Attributes(j).Name;
            attrValue = info.Datasets(i).Attributes(j).Value;
            disp(['    ', attrName, ': ', num2str(attrValue)]);
        end
    end
end
