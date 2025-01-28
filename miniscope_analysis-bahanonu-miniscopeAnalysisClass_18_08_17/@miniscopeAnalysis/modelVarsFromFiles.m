function obj = modelVarsFromFiles(obj)
	% get signals and images from input folders
	% biafra ahanonu
	% branched from controllerAnalysis: 2014.08.01 [16:09:16]
	% inputs
		%
	% outputs
		%

	% changelog
		% 2017.01.14 [20:06:04] - support switched from [nSignals x y] to [x y nSignals]
	% TODO
		% ADD SUPPORT FOR EM ANALYSIS

	%========================
	% options.populationDistanceType = 'mahal';

	% get options
	% options = getOptions(options,varargin);
	% display(options)
	% unpack options into current workspace
	% fn=fieldnames(options);
	% for i=1:length(fn)
	% 	eval([fn{i} '=options.' fn{i} ';']);
	% end
	%========================

	display(repmat('#',1,21))
	display('loading files...')

	signalExtractionMethod = obj.signalExtractionMethod;
	% usrIdxChoiceStr = {'PCAICA','EM'};
	% [sel, ok] = listdlg('ListString',usrIdxChoiceStr);
	% usrIdxChoiceList = {2,1};
	% signalExtractionMethod = usrIdxChoiceStr{sel};

	optFieldnames = fieldnames(obj.filterImageOptions);
	if obj.guiEnabled==1
		usrInput = inputdlg([optFieldnames; {'numStdsForThresh'}; {'loadVarsToRam'}],...
			'Automatic filtering parameters',1,...
			[cellfun(@num2str,struct2cell(obj.filterImageOptions),'UniformOutput',false); {'3'}; {num2str(obj.loadVarsToRam)}]...
		);
		for fieldnameNo = 1:length(optFieldnames)
			obj.filterImageOptions.(optFieldnames{fieldnameNo}) = str2num(usrInput{fieldnameNo});
		end
		numStdsForThresh = str2num(usrInput{fieldnameNo+1});
		obj.loadVarsToRam = str2num(usrInput{fieldnameNo+2});
	else
		% only update the threshold
		numStdsForThresh = 3;
	end

	[fileIdxArray idNumIdxArray nFilesToAnalyze nFiles] = obj.getAnalysisSubsetsToAnalyze();
	for thisFileNumIdx = 1:nFilesToAnalyze
		try
			fileNum = fileIdxArray(thisFileNumIdx);
			obj.fileNum = fileNum;
			display(repmat('=',1,21))
			% display([num2str(thisFileNumIdx) '/' num2str(nFilesToAnalyze) ' (' num2str(fileNum) '/' num2str(nFiles) '): ' obj.fileIDNameArray{obj.fileNum}]);
	% nFolders = length(obj.dataPath);
	% for fileNum = 1:nFolders
	% 	display(repmat('-',1,7))
	% 	try
			obj.rawSignals{fileNum} = [];
			obj.rawImages{fileNum} = [];
			obj.signalPeaks{fileNum} = [];
			obj.signalPeaksArray{fileNum} = [];
			obj.nSignals{fileNum} = [];
			obj.nFrames{fileNum} = [];
			obj.objLocations{fileNum} = [];
			obj.validManual{fileNum} = [];
			obj.validAuto{fileNum} = [];
			if strmatch('#',obj.dataPath{fileNum})
				% display([num2str(fileNum) '/' num2str(nFolders) ' | skipping: ' obj.dataPath{fileNum}]);
				display([num2str(thisFileNumIdx) '/' num2str(nFilesToAnalyze) ' (' num2str(fileNum) '/' num2str(nFiles) ') | skipping: ' obj.fileIDNameArray{obj.fileNum}]);
				obj.rawSignals{fileNum} = [];
				obj.rawImages{fileNum} = [];
				continue;
			else
				% display([num2str(fileNum) '/' num2str(nFolders) ': ' obj.dataPath{fileNum}]);
				display([num2str(thisFileNumIdx) '/' num2str(nFilesToAnalyze) ' (' num2str(fileNum) '/' num2str(nFiles) '): ' obj.fileIDNameArray{obj.fileNum}]);
			end

			switch signalExtractionMethod
				case 'PCAICA'
					% [signalTraces signalImages signalPeaks signalPeaksArray] = modelGetSignalsImages(obj,'returnType','sorted');

					regexPairs = {...
						% {'_ICfilters_sorted.mat','_ICtraces_sorted.mat'},...
						% {'holding.mat','holding.mat'},..
						{obj.rawPCAICAStructSaveStr,obj.sortedICdecisionsSaveStr,obj.classifierICdecisionsSaveStr,obj.rawICfiltersSaveStr,obj.rawICtracesSaveStr},...
						{obj.sortedICdecisionsSaveStr,obj.classifierICdecisionsSaveStr,obj.rawICfiltersSaveStr,obj.rawICtracesSaveStr},...
						{obj.rawICfiltersSaveStr,obj.rawICtracesSaveStr},...
						{obj.sortedICfiltersSaveStr,obj.sortedICtracesSaveStr}...
						% {obj.rawEMStructSaveStr},...
					};
					% get list of files to load
					filesToLoad = getFileList(obj.dataPath{fileNum},regexPairs{1});
					rawFiles = 0;
					filesToLoad = [];
					fileToLoadNo = 1;
					nRegExps = length(regexPairs);
					while isempty(filesToLoad)
						filesToLoad = getFileList(obj.dataPath{obj.fileNum},regexPairs{fileToLoadNo});
						fileToLoadNo = fileToLoadNo+1;
						if fileToLoadNo>nRegExps
							break;
						end
					end
					if fileToLoadNo==2
						rawFiles = 1;
					end
					cellfun(@display,filesToLoad)
					if isempty(filesToLoad)
					% if(~exist(filesToLoad{1}, 'file'))
						display('no files!');
					    continue
					end

					% get secondary list of files to load
					% if isempty(filesToLoad)|length(filesToLoad)<3
					%     filesToLoad = getFileList(obj.dataPath{fileNum},regexPairs{2});
					%     rawFiles = 1;
					%     if isempty(filesToLoad)
					%     % if(~exist(filesToLoad{1}, 'file'))
					%     	display('no files!');
					%         continue
					%     end
					% end
					% load files in order
					for i=1:length(filesToLoad)
					    display(['loading: ' filesToLoad{i}]);
					    load(filesToLoad{i});
					end
					if exist('pcaicaAnalysisOutput','var')
						signalTraces = double(pcaicaAnalysisOutput.IcaTraces);
						% signalImages = permute(double(pcaicaAnalysisOutput.IcaFilters),[3 1 2]);
						if strcmp(pcaicaAnalysisOutput.imageSaveDimOrder,'xyz')
							signalImages = double(pcaicaAnalysisOutput.IcaFilters);
						elseif strcmp(pcaicaAnalysisOutput.imageSaveDimOrder,'zxy')
							signalImages = permute(double(pcaicaAnalysisOutput.IcaFilters),[2 3 1]);
							% inputImages = pcaicaAnalysisOutput.IcaFilters;
						else
							% inputImages = permute(double(pcaicaAnalysisOutput.IcaFilters));
							signalImages = pcaicaAnalysisOutput.IcaFilters;
						end

						clear pcaicaAnalysisOutput;
					else
						signalImages = permute(IcaFilters,[2 3 1]);
						signalTraces = IcaTraces;
						clear IcaFilters IcaTraces;
					end
					rawFiles = 1;
				case 'EM'
					regexPairs = {...
						{obj.rawEMStructSaveStr,obj.sortedEMStructSaveStr,obj.classifierEMStructSaveStr}...
					};
					% get list of files to load
					filesToLoad = getFileList(obj.dataPath{fileNum},regexPairs{1});
					if isempty(filesToLoad)
				    	display('no files!');
				        continue
					end
					% load files in order
					for i=1:length(filesToLoad)
					    display(['loading: ' filesToLoad{i}]);
					    load(filesToLoad{i});
					end
					% signalImages = permute(emAnalysisOutput.cellImages,[3 1 2]);
					signalImages = emAnalysisOutput.cellImages;
					if isfield(emAnalysisOutput,'dsCellTraces')
						if length(emAnalysisOutput.dsCellTraces)==1
							signalTraces = emAnalysisOutput.cellTraces;
						else
							signalTraces = emAnalysisOutput.dsCellTraces;
						end
					else
						signalTraces = emAnalysisOutput.cellTraces;
					end
					% size(signalTraces)
					rawFiles = 1;
				case 'EXTRACT'
					regexPairs = {...
						{obj.rawEXTRACTStructSaveStr,obj.sortedEXTRACTStructSaveStr,obj.classifierEXTRACTStructSaveStr}...
					};
					% get list of files to load
					filesToLoad = getFileList(obj.dataPath{fileNum},regexPairs{1});
					if isempty(filesToLoad)
				    	display('no files!');
				        continue
					end
					% load files in order
					for i=1:length(filesToLoad)
					    display(['loading: ' filesToLoad{i}]);
					    load(filesToLoad{i});
					end
					signalImages = double(extractAnalysisOutput.filters);
					% signalImages = double(permute(extractAnalysisOutput.filters,[3 1 2]));
					% signalTraces = double(permute(extractAnalysisOutput.traces, [2 1]));
					signalTraces = double(extractAnalysisOutput.traces);
					% size(signalTraces)
					% size(signalImages)
					% class(signalTraces)
					% class(signalImages)
					rawFiles = 1;
				case 'CNMF'
					regexPairs = {...
						{obj.rawCNMFStructSaveStr,obj.sortedCNMFStructSaveStr,obj.classifierCNMFStructSaveStr}...
					};
					% get list of files to load
					filesToLoad = getFileList(obj.dataPath{fileNum},regexPairs{1});
					if isempty(filesToLoad)
				    	display('no files!');
				        continue
					end
					% load files in order
					for i=1:length(filesToLoad)
					    display(['loading: ' filesToLoad{i}]);
					    load(filesToLoad{i});
					end
					% signalImages = double(permute(cnmfAnalysisOutput.extractedImages,[3 1 2]));
					signalImages = double(cnmfAnalysisOutput.extractedImages);
					% signalTraces = double(cnmfAnalysisOutput.extractedSignals);
					signalTraces = double(cnmfAnalysisOutput.extractedSignalsEst);
					rawFiles = 1;
				otherwise
					% body
			end
			% if manually sorted signals, add
			if exist('valid','var')|exist('validCellMax','var')|exist('validEXTRACT','var')
				display('adding manually sorted values...')
				if exist('valid','var')
					obj.validManual{fileNum} = valid;
					obj.valid{fileNum}.(obj.signalExtractionMethod).manual = valid;
					clear valid;
				end
				if exist('validCellMax','var')
					obj.validManual{fileNum} = validCellMax;
					obj.valid{fileNum}.(obj.signalExtractionMethod).manual = validCellMax;
					clear validCellMax;
				end
				if exist(obj.validEXTRACTStructVarname,'var')
					obj.validManual{fileNum} = validEXTRACT;
					obj.valid{fileNum}.(obj.signalExtractionMethod).manual = validEXTRACT;
					clear validEXTRACT;
				end
				if exist(obj.validCNMFStructVarname,'var')
					obj.validManual{fileNum} = validEXTRACT;
					obj.valid{fileNum}.(obj.signalExtractionMethod).manual = validEXTRACT;
					clear validEXTRACT;
				end
				display('clearing manual variable...')
			end
			if exist('validClassifier','var')
				display('adding classifier annotation for signals...')
				obj.valid{fileNum}.(obj.signalExtractionMethod).classifier = validClassifier;
				display('clearing manual variable...')
				clear validClassifier;
			end

			% compute peaks
			if exist('signalTraces','var')
				% [obj.signalPeaks{fileNum}, obj.signalPeaksArray{fileNum}] = computeSignalPeaks(signalTraces, 'makePlots', 0,'makeSummaryPlots',0);
				[~, obj.signalPeaksArray{fileNum}] = computeSignalPeaks(signalTraces, 'makePlots', 0,'makeSummaryPlots',0,'numStdsForThresh',numStdsForThresh);
				obj.nSignals{fileNum} = size(signalTraces,1);
				obj.nFrames{fileNum} = size(signalTraces,2);
			end

			% rawFiles
			if rawFiles==1
				if exist('signalTraces','var')
					signalTracesTmp = signalTraces;
				else
					signalTracesTmp = [];
				end

				% % get list of movies
				% movieList = getFileList(obj.dataPath{fileNum}, obj.fileFilterRegexp);
				% [inputMovie o m n] = loadMovieList(movieList);
				% signalImagesTmp = NaN([size(signalImages,1) 20 20]);
				% for imageNo = 1:size(signalImages,1)
				% 	[signalImagesTmp(imageNo,:,:)] = viewMontage(inputMovie,signalImages(imageNo,:,:),signalTraces(imageNo,:),obj.signalPeaksArray{fileNum});
				% end
				% signalImages = signalImagesTmp;
				display(['traces dims: ' num2str(size(signalTracesTmp))])
				display(['images dims: ' num2str(size(signalImages))])
				[~, ~, validAuto, imageSizes, imgFeatures] = filterImages(signalImages, signalTracesTmp,'featureList',obj.classifierImageFeaturesNames,'options',obj.filterImageOptions);

				obj.classifierFeatures{fileNum}.(obj.signalExtractionMethod).imageFeatures = imgFeatures;
				% obj.classifierFeatures{fileNum}.signalFeatures = ;

				    [figHandle figNo] = openFigure(98, '');
					set(figHandle,'PaperUnits','inches','PaperPosition',[0 0 16 9])
				    obj.modelSaveImgToFile([],'objSize_','current',[]);
				    [figHandle figNo] = openFigure(1997, '');
					set(figHandle,'PaperUnits','inches','PaperPosition',[0 0 16 9])
				    obj.modelSaveImgToFile([],'objFeatures_','current',[]);

				% classify signals with a classifier

				% [filterImageGroups] = groupImagesByColor(signalImages,validAuto+1);
				% obj.rawImagesFiltered{fileNum} = createObjMap(filterImageGroups);
				size(validAuto)
				% validAuto
				display(['adding valid{' num2str(fileNum) '}.' obj.signalExtractionMethod '.auto identifications...'])
				obj.validAuto{fileNum} = validAuto;
				obj.valid{fileNum}.(obj.signalExtractionMethod).auto = validAuto;
				clear validAuto
				% [figHandle figNo] = openFigure(2014+round(rand(1)*100), '');
				%     imagesc(filterImageGroups);
				%     colormap(customColormap([]));
				%     box off; axis off;
				%     % colorbar
			end

			% get the x/y coordinates
			if isempty(signalImages);continue;end;
			[xCoords yCoords] = findCentroid(signalImages,'thresholdValue',0.8,'imageThreshold',0.3);
			obj.objLocations{fileNum}.(obj.signalExtractionMethod) = [xCoords(:) yCoords(:)];

			% add files
			if obj.loadVarsToRam == 1
				display('Loading variables into ram.')
			    if exist('signalTraces','var')
			    	obj.rawSignals{fileNum} = signalTraces;
			    end
			    if exist('signalImages','var')
			    	obj.rawImages{fileNum} = signalImages;
		    	end
		    else
		    	obj.rawSignals{fileNum} = [];
		    	obj.rawImages{fileNum} = [];
		    end
		    clear signalTraces signalImages
		catch err
			display(repmat('@',1,7))
			disp(getReport(err,'extended','hyperlinks','on'));
			display(repmat('@',1,7))
		end
	end
	obj.guiEnabled = 0;
	obj.modelModifyRegionAnalysis();
	obj.guiEnabled = 1;
end
function [croppedPeakImages2] = viewMontage(inputMovie,inputImage,thisTrace,signalPeakArray)

    if isempty(signalPeakArray)
        imagesc(inputImage);
        colormap(customColormap([]));
        axis off;
        croppedPeakImages2 = inputImage;
        return
    end
    % signalPeakArray
    maxSignalsToShow = 20;
    peakSignalAmplitude = thisTrace(signalPeakArray(:));
    % peakSignalAmplitude
    [peakSignalAmplitude peakIdx] = sort(peakSignalAmplitude,'descend');
    % peakSignalAmplitude
    signalPeakArray = signalPeakArray(peakIdx);
    if length(signalPeakArray)>maxSignalsToShow
        % choose a random subset
        signalPeakArray = signalPeakArray(1:maxSignalsToShow);
    end
    signalPeakArray = {signalPeakArray};
    % signalPeakArray
    croppedPeakImages = compareSignalToMovie(inputMovie, inputImage, thisTrace,'getOnlyPeakImages',1,'waitbarOn',0,'extendedCrosshairs',0,'crosshairs',0,'signalPeakArray',signalkkPeakArray);
	croppedPeakImages = squeeze(nanmean(croppedPeakImages(:,:,2:end),3));
end