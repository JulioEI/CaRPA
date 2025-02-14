% MATLAB Script to Select HDF5 or H5 File and Output the Number of Frames

% Open file selection dialog to choose an HDF5 or H5 file
[hdf5_file, hdf5_path] = uigetfile({'*.h5;*.hdf5', 'HDF5 Files (*.h5, *.hdf5)'}, 'Select an HDF5 File');

% Check if a file was selected (if the user cancels, hdf5_file will be 0)
if hdf5_file == 0
    disp('No file selected. Exiting...');
    return;
end

% Combine path and filename
hdf5_file_path = fullfile(hdf5_path, hdf5_file);

% Check if the selected file exists
if exist(hdf5_file_path, 'file') ~= 2
    disp('The selected file does not exist.');
    return;
end

% Get the structure of the HDF5 file
info = h5info(hdf5_file_path);

% Display the structure of the HDF5 file
disp('HDF5 File Structure:');
disp(info);

% Try to find and read the "frames" attribute (or dataset)
try
    % Search the entire file structure for a "frames" attribute or dataset
    frames = -1;  % Default value if no frames found

    % Check all datasets and groups for 'frames' data
    for i = 1:length(info.Groups)
        group = info.Groups(i);
        if isfield(group, 'Datasets')
            for j = 1:length(group.Datasets)
                dataset_name = group.Datasets(j).Name;
                if contains(dataset_name, 'frames', 'IgnoreCase', true)
                    frames_data = h5read(hdf5_file_path, ['/' group.Name '/' dataset_name]);
                    frames = numel(frames_data);  % Get the number of elements
                    break;
                end
            end
        end
    end

    % Display the number of frames (if found)
    if frames ~= -1
        disp(['Number of frames: ', num2str(frames)]);
    else
        disp('No frames dataset found.');
    end
catch ME
    % If an error occurs, display the error message
    disp(['Error: ', ME.message]);
end
