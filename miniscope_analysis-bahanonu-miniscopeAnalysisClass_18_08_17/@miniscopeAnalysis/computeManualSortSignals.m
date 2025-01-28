function obj = computeManualSortSignals(obj)
	% compute peaks for all signals if not already input
	% biafra ahanonu
	% branched from controllerAnalysis: 2014.08.01 [16:09:16]
	% inputs
		%
	% outputs
		%

	% changelog
		% 2014.10.09 - finished re-implementing for behaviorAnalysis class
	% TODO
		% ADD PERSONS NAME TO THE FILE

	% =======
	options.emSaveRaw = '_emAnalysis.mat';
	options.emSaveSorted = obj.sortedEMStructSaveStr;
	options.cleanedICfiltersSaveStr = obj.sortedICfiltersSaveStr;
	options.cleanedICtracesSaveStr = obj.sortedICtracesSaveStr;
	options.cleanedICdecisionsSaveStr = obj.sortedICdecisionsSaveStr;
	% =======

	display(repmat('#',1,21))
	display('computing signal peaks...')
	[fileIdxArray idNumIdxArray nFilesToAnalyze nFiles] = obj.getAnalysisSubsetsToAnalyze();
	for thisFileNumIdx = 1:nFilesToAnalyze
		try
			fileNum = fileIdxArray(thisFileNumIdx);
			obj.fileNum = fileNum;
			% fileNum = obj.fileNum;
			display(repmat('#',1,21))
			display([num2str(thisFileNumIdx) '/' num2str(nFilesToAnalyze) ' (' num2str(fileNum) '/' num2str(nFiles) '): ' obj.fileIDNameArray{obj.fileNum}]);
			% =======
			% path to current folder
			currentFolderPath = obj.inputFolders{obj.fileNum};
			% process movie regular expression
			fileFilterRegexp = obj.fileFilterRegexp;
			% get list of movies
			movieList = getFileList(currentFolderPath, fileFilterRegexp);
			% subject information
			subject = obj.subjectNum{obj.fileNum};
			assay = obj.assay{obj.fileNum};
			subjAssayIDStr = obj.fileIDNameArray{obj.fileNum};
			folderBaseSaveStr = obj.folderBaseSaveStr{obj.fileNum};
			%
			currentFolderSaveStr = [currentFolderPath filesep obj.folderBaseSaveStr{obj.fileNum}];

			usrIdxChoiceSignalType = obj.signalExtractionMethod;
			% =======
			if ~exist('usrIdxChoiceSettings','var')|strcmp(usrIdxChoiceSettings,'per folder settings')
				% usrIdxChoiceStr = {'sorting','viewing'};
				% [sel, ok] = listdlg('ListString',usrIdxChoiceStr);
				% usrIdxChoiceSortType = usrIdxChoiceStr{sel};

			 %    usrIdxChoiceStr = {'load movie','do not load movie'};
			 %    [sel, ok] = listdlg('ListString',usrIdxChoiceStr);
			 %    usrIdxChoiceMovie = usrIdxChoiceStr{sel};

			 %    usrIdxChoiceStr = {'do not classify','classify'};
			 %    [sel, ok] = listdlg('ListString',usrIdxChoiceStr);
			 %    usrIdxChoiceClassification = usrIdxChoiceStr{sel};

			 %    usrIdxChoiceStr = {'DO NOT show ROI trace','show ROI trace'};
			 %    usrIdxChoiceStrNum = [0 1];
			 %    [sel, ok] = listdlg('ListString',usrIdxChoiceStr);
			 %    usrIdxChoiceROI = usrIdxChoiceStrNum(sel);

			 %    usrIdxChoiceStr = {'start with blank','start with auto classify'};
			 %    usrIdxChoiceStrNum = [0 1];
			 %    [sel, ok] = listdlg('ListString',usrIdxChoiceStr);
			 %    usrIdxChoiceAutoValid = usrIdxChoiceStrNum(sel);
				[settingStruct] = subfxnGetSettings('Signal sorting settings',fileFilterRegexp);
				fn=fieldnames(settingStruct);
				for i=1:length(fn)
				  eval([fn{i} '=settingStruct.' fn{i} ';']);
				end
				% usrIdxChoiceFileFilterRegexp
				obj.fileFilterRegexp = usrIdxChoiceFileFilterRegexp;
				fileFilterRegexp = obj.fileFilterRegexp;

				% usrIdxChoiceSettings = settingStruct.usrIdxChoiceSettings;
			end

			% get list of movies
			movieList = getFileList(currentFolderPath, fileFilterRegexp);

			% if ~exist('usrIdxChoiceSettings','var')
			% 	usrIdxChoiceStr = {'settings across all folders','per folder settings'};
			% 	[sel, ok] = listdlg('ListString',usrIdxChoiceStr);
			% 	usrIdxChoiceSettings = usrIdxChoiceStr{sel};
			% end
		    % =======
		    if usrIdxChoiceAutoValid==2
		    	dlgBoxMsg = sprintf('select previous decisions to load for %s',obj.fileIDNameArray{obj.fileNum});
			    [filePathDecisions,folderPathDecisions,~] = uigetfile(['.\private\tmp' filesep '*.*'],dlgBoxMsg);
			    % exit if user picks nothing
			    % if folderListInfo==0; return; end
			    tmpDecisionList = {[folderPathDecisions filesep filePathDecisions]};
			    display(['loading temp decisions: ' tmpDecisionList{1}])
			    load(tmpDecisionList{1});
			end

			switch usrIdxChoiceSignalType
				case 'PCAICA'
					[rawSignals rawImages signalPeaks signalPeaksArray] = modelGetSignalsImages(obj,'returnType','raw');
					% check if the folder has temporary decisions to load (e.g. if a crash occured)
					if usrIdxChoiceAutoValid==3
						previousDecisionList = getFileList(currentFolderPath, options.cleanedICdecisionsSaveStr);
						if ~isempty(previousDecisionList)
							display(['loading previous decisions: ' previousDecisionList{1}])
							load(previousDecisionList{1});
						end
					end
					% ioptions.minValConstant = -1;
					ioptions.minValConstant = -0.1;
					ioptions.threshold = 0.5;
				case 'EM'
					[rawSignals, rawImages, signalPeaks, signalPeaksArray] = modelGetSignalsImages(obj,'returnType','raw_CellMax');
					if usrIdxChoiceAutoValid==3
						previousDecisionList = getFileList(currentFolderPath, options.emSaveSorted);
						if ~isempty(previousDecisionList)
							display(['loading previous decisions: ' previousDecisionList{1}])
							load(previousDecisionList{1});
	                        valid = validCellMax;
	                    end
                 	end
					ioptions.minValConstant = -400;
					ioptions.threshold = 0.5;
				case 'EXTRACT'
					[rawSignals, rawImages, signalPeaks, signalPeaksArray] = modelGetSignalsImages(obj,'returnType','raw');
					if usrIdxChoiceAutoValid==3
						previousDecisionList = getFileList(currentFolderPath, obj.sortedEXTRACTStructSaveStr);
						if ~isempty(previousDecisionList)
							display(['loading previous decisions: ' previousDecisionList{1}])
							load(previousDecisionList{1});
	                        valid = validEXTRACT;
	                    end
	                end
					% valid = obj.validAuto{obj.fileNum};
					ioptions.minValConstant = -10;
					ioptions.threshold = 0.5;
				case 'CNMF'
					[rawSignals, rawImages, signalPeaks, signalPeaksArray] = modelGetSignalsImages(obj,'returnType','raw');
					if usrIdxChoiceAutoValid==3
						previousDecisionList = getFileList(currentFolderPath, obj.sortedCNMFStructSaveStr);
						if ~isempty(previousDecisionList)
							display(['loading previous decisions: ' previousDecisionList{1}])
							load(previousDecisionList{1});
	                        valid = validCNMF;
	                    end
	                end
					% valid = obj.validAuto{obj.fileNum};
					ioptions.minValConstant = -200;
					ioptions.threshold = 0.3;
				otherwise
					% body
			end

			if usrIdxChoiceAutoValid==0
				display('starting with blank decisions...')
				% valid = valid*0;valid = valid+3;
				% valid = 3*ones([1 size(rawSignals,1)]);
				valid = obj.valid{obj.fileNum}.(obj.signalExtractionMethod).auto;
				valid = valid*0;valid = valid+3;
			elseif usrIdxChoiceAutoValid==1
				display('starting with automatically sorted decisions...')
				valid = obj.valid{obj.fileNum}.(obj.signalExtractionMethod).auto;
				% [~, ~, ~, ~, valid] = modelGetSignalsImages(obj,'returnType','filtered','returnOnlyValid',1);
				% valid = obj.validAuto{obj.fileNum};
			end

			% valid
			% =======
		    % load movie?
		    if strcmp(usrIdxChoiceMovie,'load movie')
		        % load movies
		        if isempty(movieList)
		        	[filePath,folderPath,~] = uigetfile([currentFolderPath filesep '*.*'],'select movie to load');
		        	% exit if user picks nothing
		        	% if folderListInfo==0; return; end
		        	movieList = [folderPath filesep filePath];
		        end
		        [ioptions.inputMovie o m n] = loadMovieList(movieList);
		        % 'frameList',1:1000
		        tmpMovie = ioptions.inputMovie(1:10,1:10,:);
		        if nanmean(tmpMovie(:))>0.9
		            display('setting mean to zero')
		            ioptions.inputMovie = ioptions.inputMovie-1;
		        end
		    else

		    end
		    % =======
		    ioptions.signalPeaks = signalPeaks;
		    ioptions.signalPeaksArray = signalPeaksArray;
		    % ioptions.inputStr = subjAssayIDStr;
		    ioptions.valid = valid;
		    ioptions.sessionID = [folderBaseSaveStr '_' num2str(java.lang.System.currentTimeMillis)];
		    ioptions.inputStr = [num2str(thisFileNumIdx) '/' num2str(nFilesToAnalyze) ' (' num2str(fileNum) '/' num2str(nFiles) '): ' obj.folderBaseSaveStr{obj.fileNum}];
		    ioptions.showROITrace = usrIdxChoiceROI;
		    ioptions.cropSizeLength = usrIdxCropSizeLength;
		    ioptions.threshold = userIdxImageThreshold;
		    ioptions.colormap = obj.colormap;
		    % ioptions.classifierFilepath = options.classifierFilepath;
		    % ioptions.classifierType = options.classifierType;
			[rawImages rawSignals valid] = signalSorter(rawImages, rawSignals,'options',ioptions);
			% rawImages = rawImages(valid,:,:);
			% rawSignals = rawSignals(valid,:);

			% add manual sorting to object
			obj.validManual{obj.fileNum} = valid;
			% commandwindow;

			% save sorted ICs
			if strcmp(usrIdxChoiceSortType,'sorting')
				switch usrIdxChoiceSignalType
					case 'PCAICA'
						IcaFilters = rawImages;
						IcaTraces = rawSignals;
						% saveID = {options.cleanedICfiltersSaveStr,options.cleanedICtracesSaveStr,options.cleanedICdecisionsSaveStr}
						% saveVariable = {'IcaFilters','IcaTraces','valid'}
						saveID = {options.cleanedICdecisionsSaveStr}
						saveVariable = {'valid'}
						for i=1:length(saveID)
							savestring = [currentFolderSaveStr saveID{i}];
							display(['saving: ' savestring])
							save(savestring,saveVariable{i});
						end
					case 'EM'
						validCellMax = valid;
						saveID = {options.emSaveSorted}
						saveVariable = {'validCellMax'}
						for i=1:length(saveID)
							savestring = [currentFolderSaveStr saveID{i}];
							display(['saving: ' savestring])
							save(savestring,saveVariable{i});
						end
					case 'EXTRACT'
						validEXTRACT = valid;
						saveID = {obj.sortedEXTRACTStructSaveStr}
						saveVariable = {'validEXTRACT'}
						for i=1:length(saveID)
							savestring = [currentFolderSaveStr saveID{i}];
							display(['saving: ' savestring])
							save(savestring,saveVariable{i});
						end
					case 'CNMF'
						validCNMF = valid;
						saveID = {obj.sortedCNMFStructSaveStr}
						saveVariable = {obj.validCNMFStructVarname}
						for i=1:length(saveID)
							savestring = [currentFolderSaveStr saveID{i}];
							display(['saving: ' savestring])
							save(savestring,saveVariable{i});
						end
					otherwise
						% body
				end
			end

			clear rawImages rawSignals valid
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
end
function [settingStruct] = subfxnGetSettings(inputTitleStr,fileFilterRegexp)

	regSettingDefaults = struct(...
		'usrIdxChoiceSortType', {{'sorting','viewing'}},...
		'usrIdxChoiceMovie',  {{'load movie','do not load movie'}},...
		'usrIdxChoiceClassification', {{'do not classify','classify'}},...
		'usrIdxChoiceROI', {{0,1}},...
		'usrIdxChoiceAutoValid',{{0,1,2,3}},...
		'usrIdxChoiceSettings',{{'settings across all folders','per folder settings'}},...
		'usrIdxChoiceFileFilterRegexp',{{fileFilterRegexp,'turboreg','crop','manualCut','dfof','downsample','other'}},...
		'usrIdxCropSizeLength',{{15,5,10,15,20,25,30,35}},...
		'userIdxImageThreshold',{{0.5,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95}}...
	);
	regSettingStr = struct(...
	    'usrIdxChoiceSortType', {{'sorting (decisions are saved)','viewing (decisions are NOT saved)'}},...
	    'usrIdxChoiceMovie',  {{'load movie (show images cut to signal peaks)','do not load movie'}},...
	    'usrIdxChoiceClassification', {{'do not classify','classify'}},...
	    'usrIdxChoiceROI', {{'DO NOT show ROI trace','show ROI trace'}},...
	    'usrIdxChoiceAutoValid',{{'start with blank','start with auto classify','start with TEMP manually chosen classifications (e.g. backups)','start with FINISHED manually chosen classifications'}},...
	    'usrIdxChoiceSettings',{{'settings across all folders','per folder settings'}},...
	    'usrIdxChoiceFileFilterRegexp',{{fileFilterRegexp,'turboreg','crop','manualCut','dfof','downsample','other (manually enter name)'}},...
	    'usrIdxCropSizeLength',{{15,5,10,15,20,25,30,35}},...
	    'userIdxImageThreshold',{{0.5,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95}}...
	);

	% propertySettings = regSettingDefaults;

	propertyList = fieldnames(regSettingDefaults);
	nPropertiesToChange = size(propertyList,1);

	% add current property to the top of the list
	for propertyNo = 1:nPropertiesToChange
		property = char(propertyList(propertyNo));
		propertyOptions = regSettingStr.(property);
		propertySettingsStr.(property) = propertyOptions;
		% propertySettingsStr.(property);
	end

	uiListHandles = {};
	uiTextHandles = {};
	uiXIncrement = 0.03;
	uiYOffset = 0.90;
	uiTxtSize = 0.3;
	uiBoxSize = 0.4;
	[figHandle figNo] = openFigure(1337, '');
	clf
	uicontrol('Style','Text','String',inputTitleStr,'Units','normalized','Position',[0.0 uiYOffset-uiXIncrement*(0) 0.3 0.05],'BackgroundColor','white','HorizontalAlignment','Left');
	for propertyNo = 1:nPropertiesToChange
		property = char(propertyList(propertyNo));
		uiTextHandles{propertyNo} = uicontrol('Style','text','String',[property ': ' 10],'Units','normalized','Position',[0.0 uiYOffset-uiXIncrement*propertyNo+0.03 uiTxtSize 0.02],'BackgroundColor',[0.9 0.9 0.9],'ForegroundColor','black','HorizontalAlignment','Left');
		% uiTextHandles{propertyNo}.Enable = 'Inactive';
		uiListHandles{propertyNo} = uicontrol('Style', 'popup','String', propertySettingsStr.(property),'Units','normalized','Position', [uiTxtSize uiYOffset-uiXIncrement*propertyNo uiBoxSize 0.05],'Callback',@(hObject,callbackdata){set(hObject, 'Backgroundcolor', [208,229,180]/255)});
	end
	uicontrol('Style','Text','String','press enter to continue','Units','normalized','Position',[0.0 uiYOffset-uiXIncrement*(nPropertiesToChange+2) 0.3 0.05],'BackgroundColor','white','HorizontalAlignment','Left');
	% uicontrol('Style','Text','String',inputTitleStr,'Units','normalized','Position',[0.0 uiYOffset 0.15 0.05],'BackgroundColor','white','HorizontalAlignment','Left');
	pause

	for propertyNo = 1:nPropertiesToChange
		property = char(propertyList(propertyNo));
		uiListHandleData = get(uiListHandles{propertyNo});
		settingStruct.(property) = regSettingDefaults.(property){uiListHandleData.Value};
	end
	close(1337)
end
    % ostruct.inputImages{ostruct.counter} = IcaFilters;
    % ostruct.inputSignals{ostruct.counter} = IcaTraces;
    % ostruct.validArray{ostruct.counter} = valid;

    % if exist(options.classifierFilepath, 'file')&strcmp(usrIdxChoiceClassification,'classify')&0
    %     display(['loading: ' options.classifierFilepath]);
    %     load(options.classifierFilepath)
    %     options.trainingOrClassify = 'classify';
    %     ioption.classifierType = options.classifierType;
    %     ioption.trainingOrClassify = options.trainingOrClassify;
    %     ioption.inputTargets = {ostruct.validArray{ostruct.counter}};
    %     ioption.inputStruct = classifierStruct;
    %     [ostruct.classifier] = classifySignals({ostruct.inputImages{ostruct.counter}},{ostruct.inputSignals{ostruct.counter}},'options',ioption);
    %     valid = ostruct.classifier.classifications;
    %     % originalValid = valid;
    %     validNorm = normalizeVector(valid,'normRange','oneToOne');
    %     validDiff = [0 diff(valid')];
    %     %
    %     figure(100020);close(100020);figure(100020);
    %     plot(valid);hold on;
    %     plot(validDiff,'g');
    %     %
    %     % validQuantiles = quantile(valid,[0.4 0.3]);
    %     % validHigh = validQuantiles(1);
    %     % validLow = validQuantiles(2);
    %     validHigh = 0.7;
    %     validLow = 0.5;
    %     %
    %     valid(valid>=validHigh) = 1;
    %     valid(valid<=validLow) = 0;
    %     valid(isnan(valid)) = 0;
    %     % questionable classification
    %     valid(validDiff<-0.3) = 2;
    %     valid(valid<validHigh&valid>validLow) = 2;
    %     %
    %     plot(valid,'r');
    %     plot(validNorm,'k');box off;
    %     legend({'scores','diff(scores)','classification','normalized scores'})
    %     % valid
    % else
    %     display(['no classifier at: ' options.classifierFilepath])
    % end

