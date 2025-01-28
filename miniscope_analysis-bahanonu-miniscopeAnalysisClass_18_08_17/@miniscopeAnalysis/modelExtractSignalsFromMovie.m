function obj = modelExtractSignalsFromMovie(obj)
	% runs signal extraction algorithms and associated saving/other functions
	% biafra ahanonu
	% started: 2015.01.05 (date might be wrong, likely late 2014)
	% inputs
		% inputMovie - a string or a cell array of strings pointing to the movies to be analyzed (recommended). Else, [x y t] matrix where t = frames.
		% numExpectedComponents - number of expected components
	% outputs
		% cnmfAnalysisOutput - structure containing extractedImages and extractedSignals along with input parameters to the algorithm
	% READ BEFORE RUNNING
		% Get CVX from http://cvxr.com/cvx/doc/install.html
		% Run the below commands in Matlab after unzipping
		% cvx_setup
		% cvx_save_prefs (permanently stores settings)

	% changelog
		% 2016.02.19 - rewrite of code to allow non-overwrite mode, so that multiple computers can connect to the same server and process the same series of folders in parallel while automatically ignoring folders that have already been processed. Could extend to include some date-based measure for analysis re-runs
	% TODO
		%

	scnsize = get(0,'ScreenSize');
	signalExtractionMethodStr = {'PCAICA','PCAICA_old','EM','EXTRACT','CNMF','ROI'};
	currentIdx = find(strcmp(signalExtractionMethodStr,obj.signalExtractionMethod));
	signalExtractionMethodDisplayStr = {'PCAICA (Mukamel, 2009) | Hakan/Tony version','PCAICA (Mukamel, 2009) | Jerome version','CELLMax (Lacey)','EXTRACT (Hakan)','CNMF (Pnevmatikakis, 2015)','ROI - only do after running either PCAICA, CELLMax, or EXTRACT'};
	[signalIdxArray, ok] = listdlg('ListString',signalExtractionMethodDisplayStr,'ListSize',[scnsize(3)*0.4 scnsize(4)*0.4],'Name','which signal extraction method?','InitialValue',currentIdx);
	% signalIdxArray
	signalExtractionMethod = signalExtractionMethodStr(signalIdxArray);

	oldPCAICA = 0;
	if iscell(signalExtractionMethod)
		nSignalExtractMethods = length(signalExtractionMethod);
	else
		nSignalExtractMethods = 1;
	end
	for signalExtractNo = 1:nSignalExtractMethods
		switch signalExtractionMethod{signalExtractNo}
			case 'PCAICA_old'
				oldPCAICA = 1;
				signalExtractionMethod = {'PCAICA'};
			otherwise
		end
	end

	% overwriteAnalysisFileSwitch = [1 0];
	% usrIdxChoiceStr = {'overwrite analysis','DO NOT overwrite analysis'};
	% [sel, ok] = listdlg('ListString',usrIdxChoiceStr);
	% overwriteAnalysisFileSwitch = overwriteAnalysisFileSwitch(sel);

	movieSettings = inputdlg({...
			'ALL | Regular expression for movie files to extract signals from',...
			'ALL | overwrite existing analysis files? (1 = yes, 0 = no)',...
			'PCAICA | output units: std, 2norm, fl, or var',...
			'CNMF | save each parameter run to new directory? (1 = yes, 0 = no)',...
			'CNMF | iterate over parameter space? (1 = yes, 0 = no)',...
			'CNMF | only run initialization algorithm? (1 = yes, 0 = no)',...
			'CELLMax | read movie chunks from disk? (1 = yes, 0 = no)',...
			'CELLMax | fraction of total frames subset each iteration?',...
			'CELLMax | number of min iterations?',...
			'CELLMax | number of max iterations?',...
			'CELLMax | gridSpacing?',...
			'CELLMax | gridWidth?',...
			'CELLMax | elimination threshold scaled phi?',...
			'CELLMax | max square tile size?',...
			'ALL | number of parallel workers',...
			'ALL | parallel enabled',...
			'ALL | Input HDF5 dataset name',...
			'ALL | Use default options (1 = yes, 0 = no)',...
			'ALL | Runtime Matlab profiler (1 = yes, 0 = no)',...
		},...
		'Cell extraction parameters',1,...
		{...
			'downsample',...
			'1',...
			'fl',...
			'0',...
			'0',...
			'0',...
			'1',...
			'0.5',...
			'200',...
			'460',...
			'',...
			'',...
			'0.005',...
			'101',...
			num2str(feature('numCores')),...
			'1',...
			obj.inputDatasetName,...
			'0',...
			'0',...
		}...
	);
	setNo = 1;
	obj.fileFilterRegexp = movieSettings{setNo};setNo = setNo+1;
	overwriteAnalysisFileSwitch = str2num(movieSettings{setNo});setNo = setNo+1;
	pcaicaOutputUnits = movieSettings{setNo};setNo = setNo+1;
	saveEachRunNewDirSwitch = str2num(movieSettings{setNo});setNo = setNo+1;
	iterateOverParameterSpace = str2num(movieSettings{setNo});setNo = setNo+1;
	onlyRunInitialization = str2num(movieSettings{setNo});setNo = setNo+1;
	options.readMovieChunks = str2num(movieSettings{setNo});setNo = setNo+1;
	options.pctTotalFramesSubsetUse = str2num(movieSettings{setNo});setNo = setNo+1;
	options.minIters = str2num(movieSettings{setNo});setNo = setNo+1;
	options.maxIters = str2num(movieSettings{setNo});setNo = setNo+1;
	options.gridSpacing = str2num(movieSettings{setNo});setNo = setNo+1;
	options.gridWidth = str2num(movieSettings{setNo});setNo = setNo+1;
	options.threshForElim = str2num(movieSettings{setNo});setNo = setNo+1;
	options.maxSqSize = str2num(movieSettings{setNo});setNo = setNo+1;
	options.numWorkers = str2num(movieSettings{setNo});setNo = setNo+1;
	options.useParallel = str2num(movieSettings{setNo});setNo = setNo+1;
	obj.inputDatasetName = movieSettings{setNo};setNo = setNo+1;
	options.defaultOptions = str2num(movieSettings{setNo});setNo = setNo+1;
	options.profiler = str2num(movieSettings{setNo});setNo = setNo+1;

	% get files to process
	[fileIdxArray idNumIdxArray nFilesToAnalyze nFiles] = obj.getAnalysisSubsetsToAnalyze();

	if iscell(signalExtractionMethod)
		nSignalExtractMethods = length(signalExtractionMethod);
	else
		nSignalExtractMethods = 1;
	end
	% pre processing for each
	for signalExtractNo = 1:nSignalExtractMethods
		switch signalExtractionMethod{signalExtractNo}
			case 'ROI'
				%
			case 'PCAICA'
				obj.signalExtractionMethod = signalExtractionMethod{signalExtractNo};

				pcaicaPCsICsSwitchStr = {'Subject','Folder'};
				[signalIdxArray, ok] = listdlg('ListString',pcaicaPCsICsSwitchStr,'ListSize',[scnsize(3)*0.4 scnsize(4)*0.4],'Name','Select PCs/ICs by subject or folder?','InitialValue',currentIdx);
				% signalIdxArray
				pcaicaPCsICsSwitchStr = pcaicaPCsICsSwitchStr{signalIdxArray};


				% use subject or folder for PC/ICA list?
				switch pcaicaPCsICsSwitchStr
					case 'Subject'
						% create expected PC/ICs signal structure if doesn't already exist
						subjectList = unique(obj.subjectStr);
						try obj.numExpectedSignals.(obj.signalExtractionMethod).(obj.subjectStr{1});check.numExpectedSignals=1; catch; check.numExpectedSignals=0; end
						% if isempty(obj.numExpectedSignals)
						if check.numExpectedSignals==0
							for subjectNum = 1:length(subjectList)
								obj.numExpectedSignals.(obj.signalExtractionMethod).(subjectList{subjectNum}) = [];
							end
						end

						% create default [PCs ICs] list else empty
						defaultList = {};
						for subjectNum = 1:length(subjectList)
							if ~isempty(obj.numExpectedSignals.(obj.signalExtractionMethod).(subjectList{subjectNum}))
								defaultList{subjectNum} = num2str(obj.numExpectedSignals.(obj.signalExtractionMethod).(subjectList{subjectNum}));
							else
								defaultList{subjectNum} = '';
							end
						end

						% ask user for nPCs/ICs
						numExpectedSignalsArray = inputdlg(subjectList,'number of PCs/ICs to use [PCs ICs]',[1 100],defaultList);
						for subjectNum = 1:length(subjectList)
							obj.numExpectedSignals.(obj.signalExtractionMethod).(subjectList{subjectNum}) = str2num(numExpectedSignalsArray{subjectNum});
						end
					case 'Folder'
						nFolders = length(fileIdxArray);

						try obj.numExpectedSignals.(obj.signalExtractionMethod).Folders;check.numExpectedSignals=1; catch; check.numExpectedSignals=0; end
						% if isempty(obj.numExpectedSignals)
						if check.numExpectedSignals==0
							for folderNo = 1:nFolders
								obj.numExpectedSignals.(obj.signalExtractionMethod).Folders{folderNo} = [];
							end
						end

						% create default [PCs ICs] list else empty
						defaultList = {};
						for folderNo = 1:nFolders
							if ~isempty(obj.numExpectedSignals.(obj.signalExtractionMethod).Folders{folderNo})
								defaultList{folderNo} = num2str(obj.numExpectedSignals.(obj.signalExtractionMethod).Folders{folderNo});
							else
								defaultList{folderNo} = '';
							end
						end

						% ask user for nPCs/ICs and store
						numExpectedSignalsArray = inputdlg(obj.folderBasePlaneSaveStr,'number of PCs/ICs to use [PCs ICs]',[1 100],defaultList);
						for folderNo = 1:nFolders
							obj.numExpectedSignals.(obj.signalExtractionMethod).Folders{folderNo} = str2num(numExpectedSignalsArray{folderNo});
						end

					otherwise
						% body
				end

			case 'EM'
				obj.signalExtractionMethod = signalExtractionMethod{signalExtractNo};

				cellmaxIntMethod = {'grid','ica'};
				[signalIdxArray, ok] = listdlg('ListString',cellmaxIntMethod,'ListSize',[scnsize(3)*0.2 scnsize(4)*0.25],'Name','Which type of initialization method to use for CELLMax?');
				% signalIdxArray
				options.CELLMax.initMethod = cellmaxIntMethod{signalIdxArray};
				% if strcmp(options.CELLMax.initMethod,'grid')
					% get each animals grid spacing and width
					% nFiles = length(obj.rawSignals);
					subjectList = unique(obj.subjectStr(fileIdxArray));
					for thisSubjectStr=subjectList
						try
							display(repmat('=',1,21))
							thisSubjectStr = thisSubjectStr{1};
							display(thisSubjectStr);
							validFoldersIdx = find(strcmp(thisSubjectStr,obj.subjectStr));
							% filter for folders chosen by the user
							validFoldersIdx = intersect(validFoldersIdx,fileIdxArray);
							if isempty(validFoldersIdx)
								continue;
							end

							movieList = getFileList(obj.inputFolders{validFoldersIdx(1)}, obj.fileFilterRegexp);
							DFOFList.(thisSubjectStr) = loadMovieList(movieList,'convertToDouble',0,'frameList',[1:500],'inputDatasetName',obj.inputDatasetName,'treatMoviesAsContinuous',1);
						catch

						end
					end
					for thisSubjectStr=subjectList
						try
							display(repmat('=',1,21))
							thisSubjectStr = thisSubjectStr{1};
							display(thisSubjectStr);
							validFoldersIdx = find(strcmp(thisSubjectStr,obj.subjectStr));
							% filter for folders chosen by the user
							validFoldersIdx = intersect(validFoldersIdx,fileIdxArray);
							if isempty(validFoldersIdx)
								continue;
							end

							movieList = getFileList(obj.inputFolders{validFoldersIdx(1)}, obj.fileFilterRegexp);
							% DFOF = loadMovieList(movieList,'convertToDouble',0,'frameList',[1:500],'inputDatasetName',obj.inputDatasetName,'treatMoviesAsContinuous',1);
							DFOF = DFOFList.(thisSubjectStr);

							h=figure;
							maxProj=max(DFOF,[],3);
							imagesc(maxProj); hax=gca;
							title('Select a region covering one cell (best to select one near another cell). Double-click region to continue.')
							cell1=imellipse(hax);
							wait(cell1);
							img1=createMask(cell1);
							gridWidth.(thisSubjectStr)=ceil(sqrt(sum(img1(:)==1)/pi));
							title('Select the closest neighboring cell. Double-click region to continue.')
							cell2=imellipse(hax);
							wait(cell2);
							pos1=getPosition(cell1);
							pos2=getPosition(cell2);
							gridSpacing.(thisSubjectStr)=ceil(norm(pos1(1:2)-pos2(1:2)))+1;
							close(h);
							clear DFOF;
						catch

						end
					end
					gridWidth
					gridSpacing
				% else
				% 	subjectList = unique(obj.subjectStr(fileIdxArray));
				% 	for thisSubjectStr=subjectList
				% 		display(repmat('=',1,21))
				% 		thisSubjectStr = thisSubjectStr{1};
				% 		display(thisSubjectStr);
				% 		gridWidth.(thisSubjectStr) = NaN;
				% 		gridSpacing.(thisSubjectStr) = NaN;
				% 	end
				% end
			case 'EXTRACT'
				obj.signalExtractionMethod = signalExtractionMethod{signalExtractNo};
				%
			case 'CNMF'
				obj.signalExtractionMethod = signalExtractionMethod{signalExtractNo};
				% create expected components if doesn't already exist
				subjectList = unique(obj.subjectStr);
				if isempty(obj.numExpectedSignals)
					% &isfield(obj.numExpectedSignals,obj.signalExtractionMethod)
					for subjectNum = 1:length(subjectList)
						obj.numExpectedSignals.(obj.signalExtractionMethod).(subjectList{subjectNum}) = [];
					end
				end

				% create default [components] list else empty
				defaultList = {};
				for subjectNum = 1:length(subjectList)
					if isfield(obj.numExpectedSignals,obj.signalExtractionMethod)&isfield(obj.numExpectedSignals.(obj.signalExtractionMethod),subjectList{subjectNum})
						defaultList{subjectNum} = num2str(obj.numExpectedSignals.(obj.signalExtractionMethod).(subjectList{subjectNum}));
					else
						defaultList{subjectNum} = '';
					end
				end

				% ask user for nPCs/ICs
				numExpectedSignalsArray = inputdlg(subjectList,'number of components to use for each subject',[1 100],defaultList);
				for subjectNum = 1:length(subjectList)
					obj.numExpectedSignals.(obj.signalExtractionMethod).(subjectList{subjectNum}) = str2num(numExpectedSignalsArray{subjectNum});
				end

				if iterateOverParameterSpace==1

					paramSetTmp = inputdlg({...
							'merge_thr | Merging threshold (positive between 0  and 1)',...
							'noise_range | Range of normalized frequencies over which to average PSD (2 x1 vector)',...
							'nrgthr | Energy threshold (positive between 0 and 1)',...
							'beta | Weight on squared L1 norm of spatial components',...
							'maxIter | Maximum number of HALS iterations',...
							'bSiz | Expansion factor for HALS localized updates'...
						},...
						'CNMF parameters, leave field blank to ignore',1,...
						{...
							'{0.98, 0.99}',...
							'{[0.10,0.6],[0.05,0.8]}',...
							'{0.85,0.80}',...
							'{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9}',...
							'{5,10,20,30}',...
							'{1,2,3,4,5,6}'...
						}...
					);
					parameterSpaceStr = {'merge_thr','noise_range','nrgthr','beta','maxIter','bSiz'};
					% initialize based on parameter strings above
					nParams = length(paramSetTmp);
					for paramNo = 1:nParams
						if ~isempty(paramSetTmp{paramNo})
							paramSetMaster.(parameterSpaceStr{paramNo}) = eval(paramSetTmp{paramNo});
						end
					end
				end
			otherwise
				% body
		end
	end


	if options.useParallel==1
		manageParallelWorkers('setNumCores',options.numWorkers);
	end

	nFolders = length(fileIdxArray);
	for thisFileNumIdx = 1:length(fileIdxArray)
		try
			% currentDateTimeStr = char(datetime('now','TimeZone','local','Format','yyyyMMdd_HHmm'));

			fileNum = fileIdxArray(thisFileNumIdx);
			obj.fileNum = fileNum;
			display(repmat('=',1,21))
			% display([num2str(fileNum) '/' num2str(nFolders) ': ' obj.fileIDNameArray{obj.fileNum}]);
			display([num2str(thisFileNumIdx) '/' num2str(nFilesToAnalyze) ' (' num2str(fileNum) '/' num2str(nFiles) '): ' obj.fileIDNameArray{obj.fileNum}]);

			%
			fileFilterRegexp = obj.fileFilterRegexp;
			if iscell(signalExtractionMethod)
				nSignalExtractMethods = length(signalExtractionMethod);
			else
				nSignalExtractMethods = 1;
			end
			% pre processing for each
			for signalExtractNo = 1:nSignalExtractMethods
				obj.signalExtractionMethod = signalExtractionMethod{signalExtractNo};

				currentDateTimeStr = datestr(now,'yyyymmdd_HHMM','local');
				diary([obj.logSavePath filesep currentDateTimeStr '_' obj.folderBaseSaveStr{obj.fileNum} '_' signalExtractionMethod{signalExtractNo} '.log']);
				startTime = tic;
				display(repmat('-',1,7))
				display(signalExtractionMethod{signalExtractNo})

				% set the save variable names and determine whether to skip files
				thisDirSaveStr = [obj.inputFolders{obj.fileNum} filesep obj.date{obj.fileNum} '_' obj.protocol{obj.fileNum} '_' obj.fileIDArray{obj.fileNum}];
				switch signalExtractionMethod{signalExtractNo}
					case 'ROI'
						saveID = {obj.rawROItracesSaveStr};
						saveVariable = {'ROItraces'};
					case 'PCAICA'
						options.rawICfiltersSaveStr = '_ICfilters.mat';
						options.rawICtracesSaveStr = '_ICtraces.mat';
						% saveID = {options.rawICfiltersSaveStr,options.rawICtracesSaveStr,obj.rawPCAICAStructSaveStr};
						% saveVariable = {'IcaFilters','IcaTraces','pcaicaAnalysisOutput'};
						saveID = {obj.rawPCAICAStructSaveStr};
						saveVariable = {'pcaicaAnalysisOutput'};
					case 'EM'
						saveID = {obj.rawEMStructSaveStr};
						saveVariable = {'emAnalysisOutput'};
					case 'EXTRACT'
						saveID = {obj.rawEXTRACTStructSaveStr};
						saveVariable = {'extractAnalysisOutput'};
					case 'CNMF'
						saveID = {obj.rawCNMFStructSaveStr,'_paramSet.mat'};
						saveVariable = {'cnmfAnalysisOutput','saveParams'};
					otherwise
						% do nothing
				end

				% skip this analysis if files already exist
				if overwriteAnalysisFileSwitch==0
					checkSaveString = [thisDirSaveStr saveID{1}];
					if exist(checkSaveString,'file')~=0
						display('SKIPPING ANALYSIS FOR THIS FOLDER')
						continue
					end
				end

				% save temporary file to prevent file checking from starting multiple runs on the same folder
				savestring = [thisDirSaveStr saveID{1}];
				display(['saving temporary: ' savestring])
				tmpVar = 'Peace is our Profession.';
				% save(savestring,saveVariable{i},'-v7.3','emOptions');
				save(savestring,'tmpVar');

				switch signalExtractionMethod{signalExtractNo}
					case 'ROI'
						runROISignalFinder();
						saveRunTimes('roi');
					case 'PCAICA'
						runPCAICASignalFinder();
						saveRunTimes('pcaica');
					case 'EM'
						emOptions = runEMSignalFinder();
						saveRunTimes('cellmax_v3');
						clear emOptions;
					case 'EXTRACT'
						emOptions = runEXTRACTSignalFinder();
						saveRunTimes('extract');
						clear extractAnalysisOutput;
						%
					case 'CNMF'
						[cnmfOptions] = runCNMFignalFinder();
						saveRunTimes('cnmf');
						clear cnmfOptions;
					otherwise
						% body
				end
				toc(startTime)
				diary OFF;
			end
		catch err
			display(repmat('@',1,7))
			disp(getReport(err,'extended','hyperlinks','on'));
			display(repmat('@',1,7))
		end
	end

	% add information about the extracted signals to the object for later processing
	objGuiOld = obj.guiEnabled;
	obj.guiEnabled = 0;
	obj.modelVarsFromFiles();
	obj.guiEnabled = 0;
	obj.viewCreateObjmaps();
	obj.guiEnabled = objGuiOld;

	function saveRunTimes(algorithm)

		if verLessThan('matlab', '8.4.0')
			return
		end

		runtimeTablePath = [obj.dataSavePathFixed filesep 'database_processing_runtimes.csv'];
		runtimeTableExists = 0;
		if exist(runtimeTablePath,'file')
			[runtimeTable] = readExternalTable(runtimeTablePath,'delimiter',',');
			addRow = size(runtimeTable,1)+1;
			% runtimeTable = table2struct(runtimeTable);
			runtimeTable.runtime_seconds(addRow,1) = toc(startTime);
			runtimeTableExists = 1;
		else
			runtimeTable = table(0,...
				{'0000.00.00'},...
				{'00:00'},...
				{'tmp'},...
				{'tmp'},...
				0,...
				0,...
				0,...
				0,...
				0,...
				0,...
				0,...
				0,...
				0,...
				0,...
				0,...
				{'tmp'},...
				0,...
				0,...
				0,...
				{'tmp'},...
				0,...
				'VariableNames',{...
				'runtime_seconds',...
				'date',...
				'daytime',...
				'folder',...
				'algorithm',...
				'frames',...
				'width',...
				'height',...
				'parallel',...
				'workers',...
				'minIters',...
				'maxIters',...
				'maxSqSize',...
				'maxDeltaParams',...
				'gridSpacing',...
				'gridWidth',...
				'initMethod',...
				'sqSizeX',...
				'sqSizeY',...
				'numSignalsDetected',...
                'versionAlgorithm',...
                'selectRandomFrames'})
			% runtimeTable.runtime_seconds = toc(startTime);

			addRow = size(runtimeTable,1)+1;
			% runtimeTable = table2struct(runtimeTable);
			runtimeTable.runtime_seconds(addRow,1) = toc(startTime);
			runtimeTableExists = 1;
		end

		movieList = getFileList(obj.inputFolders{obj.fileNum}, fileFilterRegexp);
		movieDims = loadMovieList(movieList,'convertToDouble',0,'frameList',[],'inputDatasetName',obj.inputDatasetName,'treatMoviesAsContinuous',1,'getMovieDims',1);

		runtimeTable.date{addRow,1} = datestr(now,'yyyy.mm.dd','local');
		runtimeTable.daytime{addRow,1} = datestr(now,'HH:MM','local');
		runtimeTable.folder{addRow,1} = obj.folderBaseSaveStr{obj.fileNum};
		runtimeTable.algorithm{addRow,1} = algorithm;
		runtimeTable.frames(addRow,1) = movieDims.z;
		runtimeTable.width(addRow,1) = movieDims.x;
		runtimeTable.height(addRow,1) = movieDims.y;

		fprintf('saving %s\n',runtimeTablePath);
		writetable(runtimeTable,runtimeTablePath,'FileType','text','Delimiter',',');
		parametersToAdd = {'minIters','maxIters','maxSqSize','maxSqSize','maxDeltaParams','gridSpacing','gridWidth','initMethod','sqSizeX','sqSizeY','numSignalsDetected','selectRandomFrames'};
		if strcmp(algorithm,'cellmax_v3')
            try
                runtimeTable.parallel(addRow,1) = emOptions.useParallel;
            catch
                runtimeTable.parallel(addRow,1) = NaN;
            end
			runtimeTable.workers(addRow,1) = 7;
            fn_structdisp(emOptions);
			for parameterNo = 1:length(parametersToAdd)
				parameterStr = parametersToAdd{parameterNo};
				% check that parameter name exists
				% if any(strcmp(parameterStr,fieldnames(runtimeTable)))
				if isfield(emOptions.CELLMaxoptions,parameterStr)
				else
					continue
				end
				if ~isfield(runtimeTable,parameterStr)
					if isfield(emOptions.CELLMaxoptions,parameterStr)
						if ischar(emOptions.CELLMaxoptions.(parameterStr))
							runtimeTable.(parameterStr){addRow,1} = '';
						else
							runtimeTable.(parameterStr)(addRow,1) = NaN;
						end
					end
				end
				if isfield(emOptions.CELLMaxoptions,'maxSqSize')&&~isempty(emOptions.CELLMaxoptions.(parameterStr))
					if iscell(runtimeTable.(parameterStr))
						runtimeTable.(parameterStr){addRow,1} = emOptions.CELLMaxoptions.(parameterStr);
					else
						runtimeTable.(parameterStr)(addRow,1) = emOptions.CELLMaxoptions.(parameterStr);
					end
				else
					if iscell(runtimeTable.(parameterStr))
						runtimeTable.(parameterStr){addRow,1} = '';
					else
						runtimeTable.(parameterStr)(addRow,1) = NaN;
					end
				end
			end
		else
			runtimeTable.parallel(addRow,1) = NaN;
			runtimeTable.workers(addRow,1) = NaN;
			for parameterNo = 1:length(parametersToAdd)
				parameterStr = parametersToAdd{parameterNo};
				if isfield(runtimeTable,parameterStr)
					display([parameterStr ' | ' num2str(iscell(runtimeTable.(parameterStr)))])
					if iscell(runtimeTable.(parameterStr))
						runtimeTable.(parameterStr){addRow,1} = '';
					else
						runtimeTable.(parameterStr)(addRow,1) = NaN;
					end
				else
					try
						runtimeTable.(parameterStr)(addRow,1) = NaN;
					catch
						runtimeTable.(parameterStr){addRow,1} = '';
					end
				end
			end
		end

		% if runtimeTableExists==0
		% 	runtimeTable = struct2table(runtimeTable);
		% end
		writetable(runtimeTable,runtimeTablePath,'FileType','text','Delimiter',',');
	end
	function runROISignalFinder()

		movieList = getFileList(obj.inputFolders{obj.fileNum}, fileFilterRegexp);
		[inputMovie thisMovieSize Npixels Ntime] = loadMovieList(movieList);

		[inputSignals inputImages signalPeaks signalPeaksArray] = modelGetSignalsImages(obj,'returnType','raw_images');
		[ROItraces] = applyImagesToMovie(inputImages,inputMovie);
		clear inputMovie inputImages;
		[figHandle figNo] = openFigure(1, '');
			ROItracesTmp = ROItraces;
			ROItracesTmp(ROItracesTmp<0.1) = 0;
			imagesc(ROItracesTmp);
			ylabel('filter number');xlabel('frame');
			colormap(obj.colormap);colorbar;
			title(obj.fileIDNameArray{obj.fileNum})
			set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
			obj.modelSaveImgToFile([],'ROItraces_','current',obj.fileIDArray{obj.fileNum});

		tracesSaveDimOrder = '[signalNo frameNo]';

		% =======
		for i=1:length(saveID)
			savestring = [thisDirSaveStr saveID{i}];
			display(['saving: ' savestring])
			save(savestring,saveVariable{i},'tracesSaveDimOrder');
		end
		% =======
	end

	function [extractAnalysisOutput] = runEXTRACTSignalFinder()

		movieList = getFileList(obj.inputFolders{obj.fileNum}, fileFilterRegexp);
		[inputMovie thisMovieSize Npixels Ntime] = loadMovieList(movieList);
		inputMovie(isnan(inputMovie)) = 0;

		opts.movie_dataset = obj.inputDatasetName;
		% opts.save_to_movie_dir = 1;
		% % make larger if using 2x downsampled movie
		% opts.spat_linfilt_halfwidth = 2;
		% opts.ss_cell_size_threshold = 5;
		% opts.spat_medfilt_enabled = 0;
		% opts.trim_pixels = 0.4;
		% opts.verbos = 2;
		% opts.disableGPU = 1;

		% options.turboreg = getFxnSettings();
		% options.turboreg
		% options.datasetName = options.turboreg.datasetName;

		% settingDefaults = struct(...
		%     'movie_dataset',{{'/1','/Movie','/movie'}},...
		%     'save_to_movie_dir',  {{1,0}},...
		%     'spat_linfilt_halfwidth', {{2,5}},...
		%     'ss_cell_size_threshold', {{5,10}},...
		%     'spat_medfilt_enabled', {{0,1}},...
		%     'trim_pixels', {{0.4,0.6}},...
		%     'verbos', {{0,1}},...
		%     'disableGPU', {{1,0}}...
		% );
		% settingStr = struct(...
		%     'movie_dataset',{{'/1','/Movie','/movie'}},...
		%     'save_to_movie_dir',  {{1,0}},...
		%     'spat_linfilt_halfwidth', {{2,5}},...
		%     'ss_cell_size_threshold', {{5,10}},...
		%     'spat_medfilt_enabled', {{0,1}},...
		%     'trim_pixels', {{0.4,0.6}},...
		%     'verbos', {{0,1}},...
		%     'disableGPU', {{1,0}}...
		% );

		[h,w,t] = size(inputMovie);

		opts.max_cell_radius=30;
		opts.min_cell_spacing=5;
		opts.remove_duplicate_cells = 0;
		% Use GPU
		opts.compute_device='gpu';

		% This is how to call the function 'partition_helper()' to find out how many partitions are necessary:
		num_parts = partition_helper(h,w,t,opts.min_cell_spacing,opts.max_cell_radius);

		% Below call returned num_parts=20. We decide to partition x axis to 4, and y axis to 5. This makes 20 parititions overall.
		nPlotsRoot = sqrt(num_parts);
		if nPlotsRoot<2
			nPlotsRoot = 2;
		end
		integ = fix(nPlotsRoot);
		fract = abs(nPlotsRoot - integ);
		opts.num_partition_y = ceil(nPlotsRoot);
		opts.num_partition_x = floor(nPlotsRoot)+round(fract)

		min_cell_spacing=3;
		max_cell_radius=10;
		num_partition_x=3;
		num_partitiony=3;
		cell_keep_tolerance=5;
		subtract_background=1;

		opts.config.diffuse_F=1;
		opts.config.smooth_T = 0;
		opts.config.smooth_F = 0;
		% opts.config.cell_keep_tolerance

		% [filters,traces,info,opts] = extractor(movieList{1},opts);
		% [filters,traces,info,opts] = extractor(inputMovie,opts);
		outStruct = extractor(inputMovie,opts);

		im_dup_corr_thresh = 0.05; % Image correlation threshold
		trace_dup_corr_thresh = 0.6; % Trace correlation threshold
		outStruct = remove_duplicates(outStruct,im_dup_corr_thresh,trace_dup_corr_thresh);

		extractAnalysisOutput.filters = outStruct.spatial_weights;
		% permute so it is [nCells frames]
		extractAnalysisOutput.traces = permute(outStruct.temporal_weights, [2 1]);
		extractAnalysisOutput.info = outStruct.info;
		extractAnalysisOutput.opts = outStruct.opts;
		extractAnalysisOutput.file = movieList{1};

		% =======
		% save EXTRACT signals
		for i=1:length(saveID)
			savestring = [thisDirSaveStr saveID{i}];
			display(['saving: ' savestring])
			% save(savestring,saveVariable{i},'-v7.3','emOptions');
			save(savestring,saveVariable{i});
		end
		% =======
	end
	function runPCAICASignalFinder()
		switch pcaicaPCsICsSwitchStr
			case 'Subject'
				nPCsnICs = obj.numExpectedSignals.(obj.signalExtractionMethod).(obj.subjectStr{obj.fileNum})
			case 'Folder'
				nPCsnICs = obj.numExpectedSignals.(obj.signalExtractionMethod).Folders{obj.fileNum}
			otherwise
				% body
		end
		% return;

		nPCs = nPCsnICs(1);
		nICs = nPCsnICs(2);
		%
		movieList = getFileList(obj.inputFolders{obj.fileNum}, fileFilterRegexp);
		% [inputMovie] = loadMovieList(movieList,'convertToDouble',0,'frameList',[]);

		if oldPCAICA==1
			display('running PCA-ICA, old version...')
			[PcaFilters PcaTraces] = runPCA(movieList, '', nPCs, fileFilterRegexp);
			if isempty(PcaFilters)
				display('PCs are empty, skipping...')
				return;
			end
			[IcaFilters IcaTraces] = runICA(PcaFilters, PcaTraces, '', nICs, '');
			traceSaveDimOrder = '[nComponents frames]';
			% reorder if needed
			options.IcaSaveDimOrder = 'xyz';
			if strcmp(options.IcaSaveDimOrder,'xyz')
				IcaFilters = permute(IcaFilters,[2 3 1]);
				imageSaveDimOrder = 'xyz';
			else
				imageSaveDimOrder = 'zxy';
			end
		else
			display('running PCA-ICA, new version...')
			[PcaOutputSpatial PcaOutputTemporal PcaOutputSingularValues PcaInfo] = run_pca(movieList, nPCs, 'movie_dataset_name',obj.inputDatasetName);

			if isempty(PcaOutputTemporal)
				display('PCs are empty, skipping...')
				return;
			end

			display('+++')
			movieDims = loadMovieList(movieList,'convertToDouble',0,'frameList',[],'inputDatasetName',obj.inputDatasetName,'treatMoviesAsContinuous',1,'getMovieDims',1);

			% output_units = 'fl';
			% output_units = 'std';
			[IcaFilters, IcaTraces, IcaInfo] = run_ica(PcaOutputSpatial, PcaOutputTemporal, PcaOutputSingularValues, movieDims.x, movieDims.y, nICs, 'output_units',pcaicaOutputUnits);
			IcaTraces = permute(IcaTraces,[2 1]);
			traceSaveDimOrder = '[nComponents frames]';
			% reorder if needed
			options.IcaSaveDimOrder = 'xyz';
			if strcmp(options.IcaSaveDimOrder,'xyz')
				imageSaveDimOrder = 'xyz';
			else
				IcaFilters = permute(IcaFilters,[3 1 2]);
				imageSaveDimOrder = 'zxy';
			end
			pcaicaAnalysisOutput.IcaInfo = IcaInfo;
		end


		pcaicaAnalysisOutput.IcaFilters = IcaFilters;
		pcaicaAnalysisOutput.IcaTraces = IcaTraces;
		pcaicaAnalysisOutput.imageSaveDimOrder = imageSaveDimOrder;
		pcaicaAnalysisOutput.traceSaveDimOrder = traceSaveDimOrder;
		pcaicaAnalysisOutput.nPCs = nPCs;
		pcaicaAnalysisOutput.nICs = nICs;
		pcaicaAnalysisOutput.time.startTime = startTime;
		pcaicaAnalysisOutput.time.endTime = toc(startTime);
		pcaicaAnalysisOutput.time.dateTime = datestr(now,'yyyymmdd_HHMM','local');
		pcaicaAnalysisOutput.movieList = movieList;
		% =======
		% save ICs
		saveOutput = 1;
		if saveOutput==1
			for i=1:length(saveID)
				savestring = [thisDirSaveStr saveID{i}];
				display(['saving: ' savestring])
				save(savestring,saveVariable{i});
			end
		end
		% =======
	end
	function [emOptions] = runEMSignalFinder()
		% emOptions.dsMovieDatasetName = options.datasetName;
		% emOptions.movieDatasetName = options.datasetName;
		movieList = getFileList(obj.inputFolders{obj.fileNum}, fileFilterRegexp);
        for k = 1:length(movieList)
            % movieFilename=[];
            % upsampledMovieList = getFileList(thisDir, fileFilterRegexp);
            % mpiprofile on
            % emOptions.CELLMaxoptions = emOptions.EMoptions;
            display(['input movie: ' movieList{k}])

            % =====================
            clear emOptions;

            if strcmp(options.CELLMax.initMethod,'grid')
                emOptions.CELLMaxoptions.initMethod = 'grid';
            elseif strcmp(options.CELLMax.initMethod,'ica')
                emOptions.CELLMaxoptions.initMethod='ica';
            end
            emOptions.CELLMaxoptions.gridSpacing = gridSpacing.(obj.subjectStr{obj.fileNum});
            emOptions.CELLMaxoptions.gridWidth = gridWidth.(obj.subjectStr{obj.fileNum});
            if ~isempty(options.gridSpacing)
                emOptions.CELLMaxoptions.gridSpacing = options.gridSpacing;
                emOptions.CELLMaxoptions.gridWidth = options.gridWidth;
            end
            emOptions.useParallel = options.useParallel;
            emOptions.CELLMaxoptions.inputSizeManual = 0;

            emOptions.CELLMaxoptions.subsampleMethod = 'resampleRemaining';
            % [maxIters nMovieFrames]
            % options.subsampleFrameMatrix = [];
            % [1 nMovieFrames] - vector of frames to use in a movie
            % options.subsampleFrameVector = [];
            % options.selectRandomFrames=1;
            % options.numFramesRandom=2000;
            % 0 to 1, percentage of frames per iteration to select
            emOptions.CELLMaxoptions.percentFramesPerIteration = options.pctTotalFramesSubsetUse;
            % subsampleMethod = 'resampleRemaining', fr
            emOptions.CELLMaxoptions.percentRemainingSubsample = 0.75;%1; CUSTOM_CODE!!
            emOptions.CELLMaxoptions.maxSqSize = options.maxSqSize;
            emOptions.CELLMaxoptions.threshForElim = options.threshForElim;

            if options.defaultOptions==0
                emOptions.CELLMaxoptions.localICimgs = [];
                emOptions.CELLMaxoptions.localICtraces = [];
                emOptions.CELLMaxoptions.minIters = options.minIters;
                emOptions.CELLMaxoptions.maxIters = options.maxIters;
                emOptions.CELLMaxoptions.inputSizeManual = 0;
                emOptions.CELLMaxoptions.numSigmasThresh = 0.5;
                emOptions.CELLMaxoptions.nParallelWorkers = options.numWorkers;
                emOptions.CELLMaxoptions.generateNovelSeed = 1;
                % emOptions.CELLMaxoptions.randNumGenSeed = 2;
                movieDims = loadMovieList(movieList{k},'getMovieDims',1,'inputDatasetName',obj.inputDatasetName);
                emOptions.CELLMaxoptions.numFramesRandom = round(movieDims.z*options.pctTotalFramesSubsetUse);
                if emOptions.CELLMaxoptions.numFramesRandom<3e3
                    emOptions.CELLMaxoptions.numFramesRandom = 3e3;
                end
                emOptions.CELLMaxoptions.readMovieChunks = options.readMovieChunks;
            end
            emOptions.movieDatasetName=obj.inputDatasetName;
            emOptions.CELLMaxoptions.movieFilename = movieList{k};
            % =====================

            fn_structdisp(emOptions);

            if options.profiler==1
                currentDateTimeStr = datestr(now,'yyyymmdd_HHMM','local');
                profilerSaveLocation = [obj.inputFolders{obj.fileNum} filesep 'profilerCELLMax_' currentDateTimeStr];
                display(['Profiler will be saved to: ' profilerSaveLocation])
                profile on
            end

            [emAnalysisOutput, ~] = CELLMax_Wrapper(movieList{k},'options',emOptions);

            if options.profiler==1
                profile off
                profsave(profile('info'),profilerSaveLocation);
            end
            % [emAnalysisOutput, ~] = EM_CellFind_Wrapper(movieList{1},[],'options',emOptions);
            % emOptions.CELLMaxoptions.sqSizeX = NaN;
            % emOptions.CELLMaxoptions.sqSizeY = NaN;

            emOptions.CELLMaxoptions.sqSizeX = [];
            emOptions.CELLMaxoptions.sqSizeY = [];
            emAnalysisOutput.dsCellTraces = emAnalysisOutput.cellTraces;
            emOptions.CELLMaxoptions.numSignalsDetected = size(emAnalysisOutput.dsCellTraces,1);
            emOptions.versionCellmax = emAnalysisOutput.versionCellmax;
            % emOptions.EMoptions = emOptions.CELLMaxoptions;
            % mpiprofile off
            % mpiprofile viewer
            % pause

            % output.cellImages : images representing sources found (candidate cells). not all will be cells. Size is [x y numCells]
            % output.centroids : centroids of each cell image, x (horizontal) and then y (vertical). Size is [numCells 2]
            % output.convexHulls : convex hull (line tracing around edge) of each cell, in x then y. Cell Array, Size is [numCells 1], each entry is hull of one cell.
            % output.dsEventTimes : event timings on the down sampled probability traces.
            % output.dsScaledProbabilities : a scaled probability trace for each cell, from the downsampled movie. Can be used as a denoised fluorescence trace.
            % output.dsCellTraces : fluorescence traces for each cell, from the temporally downsampled movie. Size is [numCells numFrames] for numFrames of downsampled movie
            % output.cellTraces : fluorescence traces for each cell, from the full temporal resolution movie. Size is [numCells numFrames] for numFrames of full movie
            % output.eventTimes : event timings as output by detectEvents.
            % output.EMoptions : options that EM was run with. Good to keep for recordkeeping purposes.

            emOptions.time.startTime = startTime;
            emOptions.time.endTime = toc(startTime);
            emAnalysisOutput

            % =======
            % save output components
            
%             for i=1:length(saveID)
%                 savestring = [thisDirSaveStr saveID{i}];
%                 display(['saving: ' savestring])
%                 save(savestring,saveVariable{i},'-v7.3','emOptions');
%                 % save(savestring,saveVariable{i},'emOptions');
%             end
            for i=1:length(saveID) %CUSTOM CODE!
                splitSaveID = split(saveID{i},'.');
                saveIDK = join([splitSaveID(1),'_',num2str(k),'.',splitSaveID(2)],'');
                savestring = join([thisDirSaveStr saveIDK],'');
                savestring = savestring{1};
                display(['saving: ' savestring])
                save(savestring,saveVariable{i},'-v7.3','emOptions');
                % save(savestring,saveVariable{i},'emOptions');
            end
            % =======
        end
    end
	function [cnmfOptions] = runCNMFignalFinder()

		% if cvx is not in the path, ask user
		if isempty(which('cvx_begin'))
			[filePath,folderPath,~] = uigetfile(['*.*'],'select cvx_setup.m');
			run([folderPath filesep filePath]);
		end

		movieList = getFileList(obj.inputFolders{obj.fileNum}, fileFilterRegexp);

		% Y = loadMovieList(movieList,'convertToDouble',0,'frameList',[],'inputDatasetName',obj.inputDatasetName,'treatMoviesAsContinuous',1);

		% Merging threshold (positive between 0  and 1)
		cnmfOptions.merge_thr = 0.95;
		% Range of normalized frequencies over which to average PSD (2 x1 vector)
		cnmfOptions.noise_range = [0.10,0.6];
		% Size of 2-d median filter (2 x 1 array of positive integers)
		cnmfOptions.medw = [2,2];
		% Energy threshold (positive between 0 and 1)
		% cnmfOptions.nrgthr = 0.97;
		cnmfOptions.nrgthr = 0.85;
		% Morphological closing operator for post-processing (binary image)
		cnmfOptions.clos_op = strel('disk',3,0);
		% Flag for computing noise values sequentially for memory reasons
		cnmfOptions.split_data = 0;
		% Spatial down-sampling factor (scalar >= 1)
		cnmfOptions.ssub = 1;
		% Temporal down-sampling factor (scalar >= 1)
		cnmfOptions.tsub = 1;
		% Maximum number of sparse NMF iterations
		cnmfOptions.snmf_max_iter = 200;
		% Weight on squared L1 norm of spatial components
		cnmfOptions.beta = 0.5;
		% Weight on frobenius norm of temporal components * max(Y)^2
		cnmfOptions.eta = 1;
		% Relative change threshold for stopping sparse_NMF
		cnmfOptions.err_thr = 1e-4; %1e-5
		% Maximum number of HALS iterations
		cnmfOptions.maxIter = 5;
		% Expansion factor for HALS localized updates
		cnmfOptions.bSiz = 3; %1e-5

		% Standard deviation of Gaussian kernel for initialization
		cnmfOptions.otherCNMF.tau = 2;
		% cnmfOptions.otherCNMF.tau = [3 2 1; 3 2 1]';
		% Order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
		cnmfOptions.otherCNMF.p = 2;

		if onlyRunInitialization==1
			% run only initialization algorithm
			cnmfOptions.nonCNMF.onlyRunInitialization = 1;
		else
			% run only initialization algorithm
			cnmfOptions.nonCNMF.onlyRunInitialization = 0;
		end

		% set of parameters to vary
		if iterateOverParameterSpace==1
			paramSet = paramSetMaster;
		else
			paramSet.beta = {0.1};
			paramSet.maxIter = {5};
			paramSet.bSiz = {1};
		end

		% make a list of all possible combinations
		paramNames = fieldnames(paramSet);
		nParams = length(paramNames);
		paramStr = 'paraSpace = combvec(';
		for paramNo = 1:nParams
			paramStr = [paramStr '1:' num2str(length(paramSet.(paramNames{paramNo})))];
			if paramNo~=nParams
				paramStr = [paramStr ','];
			end
		end
		paramStr = [paramStr ');'];
		eval(paramStr);

		for paramNo = 1:nParams
			paramIdx = paraSpace(paramNo,:);
			paramSet.(paramNames{paramNo}) = paramSet.(paramNames{paramNo})(paramIdx(:));
		end
		% [p,q,r] = meshgrid(1:length(paramSet.merge_thr),1:length(paramSet.noise_range),1:length(paramSet.nrgthr));
		% paraSpace = [p(:) q(:) r(:)];
		% paraSpace = combvec(1:2,1:4,1:3);
		% paramSet.merge_thr = paramSet.merge_thr(p(:));
		% paramSet.noise_range = paramSet.noise_range(q(:));
		% paramSet.nrgthr = paramSet.nrgthr(r(:));

		% paramSet.tau = {};
		% paramSet.p = {};
		nParameterSets = size(paraSpace,2);
		% decide whether to iterate over new parameters
		iterateParameters = 0;
		saveParams.null = 1;

		if iterateParameters==0
			nParameterSets = 1;
		end
		for parameterSetNo = 1:nParameterSets
			try
				display(repmat('*',1,14))
				% display([num2str(fileNum) '/' num2str(nFolders) ': ' obj.fileIDNameArray{obj.fileNum}]);
				display(['parameter set:' num2str(parameterSetNo) '/' num2str(nParameterSets)]);

				if iterateParameters==1
					% add the iterated parameter here
					paramNames = fieldnames(paramSet);
					nParams = length(paramNames);
					for paramNo = 1:nParams
						cnmfOptions.(paramNames{paramNo}) = paramSet.(paramNames{paramNo}){parameterSetNo};
						saveParams.(paramNames{paramNo}) = paramSet.(paramNames{paramNo}){parameterSetNo};
					end
					saveParams
				end

				startTime = tic;
				cnmfAnalysisOutput = [];
				[cnmfAnalysisOutput] = computeCnmfSignalExtraction(movieList,obj.numExpectedSignals.(obj.signalExtractionMethod).(obj.subjectStr{obj.fileNum}),'options',cnmfOptions);
				% [cnmfAnalysisOutput] = computeCnmfSignalExtractionOriginal(movieList,obj.numExpectedSignals.(obj.signalExtractionMethod).(obj.subjectStr{obj.fileNum}),'options',cnmfOptions);

				obj.modelSaveImgToFile([],'initializationROIs_',1337,[obj.folderBaseSaveStr{obj.fileNum} '_run0' num2str(parameterSetNo)]);
				[figHandle figNo] = openFigure(1337, '');hold off;
				obj.modelSaveImgToFile([],'cellmapContours_',1339,[obj.folderBaseSaveStr{obj.fileNum} '_run0' num2str(parameterSetNo)]);
				[figHandle figNo] = openFigure(1339, '');hold off;

				cnmfAnalysisOutput.time.startTime = startTime;
				cnmfAnalysisOutput.time.endTime = tic;
				cnmfAnalysisOutput.time.totalTime = toc(startTime);

				% save ICs
				if saveEachRunNewDirSwitch==0
					saveID = {obj.rawCNMFStructSaveStr,'_paramSet.mat'};
					saveVariable = {'cnmfAnalysisOutput','saveParams'};
					thisDirSaveStr = [obj.inputFolders{obj.fileNum} filesep obj.date{obj.fileNum} '_' obj.protocol{obj.fileNum} '_' obj.fileIDArray{obj.fileNum}];
				else
					saveID = {obj.rawCNMFStructSaveStr,'_paramSet.mat'};
					saveVariable = {'cnmfAnalysisOutput','saveParams'};
					thisDirSaveStr = [obj.inputFolders{obj.fileNum} filesep 'run0' num2str(parameterSetNo) filesep];
					if (~exist(thisDirSaveStr,'dir')) mkdir(thisDirSaveStr); end;
					thisDirSaveStr = [thisDirSaveStr obj.date{obj.fileNum} '_' obj.protocol{obj.fileNum} '_' obj.fileIDArray{obj.fileNum}];
				end

				% =======
				for i=1:length(saveID)
					savestring = [thisDirSaveStr saveID{i}];
					display(['saving: ' savestring])
					% save(savestring,saveVariable{i},'-v7.3','emOptions');
					save(savestring,saveVariable{i});
				end
				% =======
			catch err
				display(repmat('@',1,7))
				disp(getReport(err,'extended','hyperlinks','on'));
				display(repmat('@',1,7))
			end
		end
	end
end
function [turboregSettingStruct] = getFxnSettings()

	% propertySettings = turboregSettingDefaults;

	propertyList = fieldnames(turboregSettingDefaults);
	nPropertiesToChange = size(propertyList,1);

	% add current property to the top of the list
	for propertyNo = 1:nPropertiesToChange
		property = char(propertyList(propertyNo));
		propertyOptions = turboregSettingStr.(property);
		propertySettingsStr.(property) = propertyOptions;
		% propertySettingsStr.(property);
	end

	uiListHandles = {};
	uiTextHandles = {};
	uiXIncrement = 0.05;
	uiYOffset = 0.95;
	uiTxtSize = 0.3;
	uiBoxSize = 0.4;
	[figHandle figNo] = openFigure(1337, '');
	clf
	uicontrol('Style','Text','String','processing options','Units','normalized','Position',[0.0 uiYOffset-uiXIncrement*(nPropertiesToChange+1) 0.3 0.05],'BackgroundColor','white','HorizontalAlignment','Left');
	for propertyNo = 1:nPropertiesToChange
		property = char(propertyList(propertyNo));
		uiTextHandles{propertyNo} = uicontrol('Style','Text','String',[property ': '],'Units','normalized','Position',[0.0 uiYOffset-uiXIncrement*propertyNo uiTxtSize 0.05],'BackgroundColor','white','HorizontalAlignment','Left');
		uiListHandles{propertyNo} = uicontrol('Style', 'popup','String', propertySettingsStr.(property),'Units','normalized','Position', [uiTxtSize uiYOffset-uiXIncrement*propertyNo uiBoxSize 0.05]);
	end
	uicontrol('Style','Text','String','press enter to continue','Units','normalized','Position',[0.0 uiYOffset-uiXIncrement*(nPropertiesToChange+1) 0.3 0.05],'BackgroundColor','white','HorizontalAlignment','Left');
	pause

	for propertyNo = 1:nPropertiesToChange
		property = char(propertyList(propertyNo));
		uiListHandleData = get(uiListHandles{propertyNo});
		turboregSettingStruct.(property) = turboregSettingDefaults.(property){uiListHandleData.Value};
	end
	close(1337)
end