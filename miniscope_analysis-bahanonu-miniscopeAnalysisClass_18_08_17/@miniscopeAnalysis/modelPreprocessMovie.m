function obj = modelPreprocessMovie(obj)
% preprocess movies


	% verify that there are folders present
	if isempty(obj.inputFolders)
		uiwait(msgbox('Please run modelAddNewFolders to add folders to the object','Success','modal'));
		return
	end
	% check that Miji exists, if not, have user enter information
	% try
	% 	Miji
	% 	MIJ.exit
	% catch
	if exist('Miji.m','file')==2
		display(['Miji located in: ' which('Miji.m')]);
		% Miji is loaded, continue
	else
		% pathToMiji = inputdlg('Enter path to Miji.m in Fiji (e.g. \Fiji.app\scripts):',...
		%              'Miji path', [1 100]);
		% pathToMiji = pathToMiji{1};
		pathToMiji = uigetdir('\.','Enter path to Miji.m in Fiji (e.g. \Fiji.app\scripts)');
		if ischar(pathToMiji)
			privateLoadBatchFxnsPath = 'private\privateLoadBatchFxns.m';
			if exist(privateLoadBatchFxnsPath,'file')~=0
				fid = fopen(privateLoadBatchFxnsPath,'at')
				fprintf(fid, '\npathtoMiji = ''%s'';\n', pathToMiji);
				fclose(fid);
			end
			addpath(pathToMiji);
		end
	end


	if obj.guiEnabled==1
		scnsize = get(0,'ScreenSize');
		[fileIdxArray, ok] = listdlg('ListString',obj.fileIDNameArray,'ListSize',[scnsize(3)*0.2 scnsize(4)*0.25],'Name','which folders to analyze?');
	else
		if isempty(obj.foldersToAnalyze)
			fileIdxArray = 1:length(obj.fileIDNameArray);
		else
			fileIdxArray = obj.foldersToAnalyze;
		end
	end

	options.fileFilterRegexp = 'concat_.*.h5';
	folderListInfo = {obj.inputFolders{fileIdxArray}};
	options.datasetName = obj.inputDatasetName;

	% controllerPreprocessMovie2('folderListPath',folderListInfo,'fileFilterRegexp',options.fileFilterRegexp,'datasetName',options.datasetName,'frameList',[]);
	obj.modelPreprocessMovieFunction('folderListPath',folderListInfo,'fileFilterRegexp',options.fileFilterRegexp,'datasetName',options.datasetName,'frameList',[]);
end