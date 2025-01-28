function loadRepoFunctions()
	% loads the necessary directories to have the batch functions present
	% biafra ahanonu
	% started: 2013.12.24 [08:46:11]

	% changelog
		% updated: 2015.08.24 [23:17:44] - removed Miji and other specific aspects of the function for the moment. Renaming function as well to be more appropriate to its function.
		% 2015.10.08 [21:01:19] - included a filter for .git folder to avoid adding that to Matlab path

	% add controller directory and subdirectories to path, DOES NOT add any folder names /private
	pathList = genpath(pwd);
	pathListArray = strsplit(pathList,';');
	pathFilter = cellfun(@isempty,regexpi(pathListArray,'\.git'));
	pathListArray = pathListArray(pathFilter);
	pathList = strjoin(pathListArray,';');
% 	addpath(pathList);

	setFigureDefaults();
end