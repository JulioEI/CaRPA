classdef miniscopeAnalysis < dynamicprops
	% Performs analysis on behavior (response signals) compared to stimulus or other continuous signals during a trial.
	% biafra ahanonu
	% started: 2014.07.31
	% updated: 2017.01.15 [01:31:54]
	% This is a re-write of old code from controllerAnalysis. Encapsulating the functions and variables as methods and properties in a class allows easier maintenance/flexibility.
	% inputs
		%
	% outputs
		%

	% changelog
		%
	% TODO
		%

	% dynamicprops is a subclass of handle, allowing addition of properties

	properties(GetAccess = 'public', SetAccess = 'public')
		% public read and write access.

		defaultObjDir = pwd;
		% 0 = load variables from disk, reduce RAM usage. 1 = load from disk to ram, faster for analysis.
		loadVarsToRam = 0;
		% show GUI for view functions?
		guiEnabled = 1;
		% indices for folders to analyze, [] = all
		foldersToAnalyze = [];
		% indices for stimuli to analyze, [] = all
		discreteStimuliToAnalyze = [];
		% io settings
		fileFilterRegexp = 'crop';
		% behavior video regexp
		behaviorVideoRegexp = '';
		% loop over all files during analysis? 'individual' or 'group'
		analysisType  = 'group';
		% 1 = perform certain analysis on dF/F instead of peaks
		dfofAnalysis = 0;
		% 'filtered' returns auto/manually filtered signals/images, 'raw' returns raw
		modelGetSignalsImagesReturnType = 'filtered'
		% name of input dataset name for preprocessing
		inputDatasetName = '/1';
		%
		stimTriggerOnset = 0;
		% paths for specific types of files
		currentDateTimeStr = datestr(now,'yyyymmdd','local');
		picsSavePath = ['private' filesep 'pics' filesep datestr(now,'yyyymmdd','local') filesep];
		dataSavePath = ['private' filesep 'data' filesep datestr(now,'yyyymmdd','local') filesep];
		logSavePath = ['private' filesep 'logs' filesep datestr(now,'yyyymmdd','local') filesep];
		%
		dataSaveFilenameModifier = '';
		% table save
		delimiter = ',';
		% name of i/o HDF5 dataset names
		hdf5Datasetname = '/1';
		% type of images to save analysis as '-dpng','-dmeta','-depsc2'
		imgSaveTypes = {'-dpng'};
		% colormap to be used
		% colormap = customColormap([]);
		colormap = customColormap({[0 0 1],[1 1 1],[0.5 0 0],[1 0 0]});
		% colormap = diverging_map(linspace(0,1,100),[0 0 1],[1 0 0]);
		% use for stimulus related viewing functions
		% frames before/after stimulus to look
		timeSequence = [-50:50];
		postStimulusTimeSeq = [0:10];
		% bin analysis, integer only
		binDownsampleAmount = 1;
		%
		stimulusTableValueName = 'frameSessionDownsampled';
		stimulusTableFrameName = 'frameSessionDownsampled';
		stimulusTableTimeName = 'time';
		stimulusTableSessionName = 'trial';
		%


		% methods
		currentMethod = 'modelAddNewFolders';

		% Region analysis
		regionModSaveStr = '_regionModSelectUser.mat'

		% PCAICA names
		rawPCAICAStructSaveStr = '_pcaicaAnalysis.mat';
		rawICfiltersSaveStr = '_ICfilters.mat';
		rawICtracesSaveStr = '_ICtraces.mat';
		sortedICfiltersSaveStr = '_ICfilters_sorted.mat';
		sortedICtracesSaveStr = '_ICtraces_sorted.mat';
		sortedICdecisionsSaveStr = '_ICdecisions.mat';
		classifierICdecisionsSaveStr = '_ICclassifierDecisions.mat';
		structPCAICAVarname = 'pcaicaAnalysisOutput';
		% ROI names
		rawROItracesSaveStr = '_ROItraces.mat';
		% EM names
		rawEMStructSaveStr = '_emAnalysis.mat';
		sortedEMStructSaveStr = '_emAnalysisSorted.mat';
		classifierEMStructSaveStr = '_emAnalysisClassifierDecisions.mat';
		structEMVarname = 'emAnalysisOutput';
		% validEMStructVarname = 'validCellMax';
		% EXTRACT names
		rawEXTRACTStructSaveStr = '_extractAnalysis.mat';
		sortedEXTRACTStructSaveStr = '_extractAnalysisSorted.mat';
		classifierEXTRACTStructSaveStr = '_extractAnalysisClassifierDecisions.mat';
		validEXTRACTStructVarname = 'validEXTRACT';
		structEXTRACTVarname = '';
		% EXTRACT names
		rawCNMFStructSaveStr = '_cnmfAnalysis.mat';
		sortedCNMFStructSaveStr = '_cnmfAnalysisSorted.mat';
		classifierCNMFStructSaveStr = '_cnmfAnalysisClassifierDecisions.mat';
		validCNMFStructVarname = 'validCNMF';
		structCNMRVarname = 'cnmfAnalysisOutput';
		% PCAICA, EM, or EXTRACT
		signalExtractionMethod = 'PCAICA';%EM

		settingOptions = struct(...
			'analysisType',  {{'group','individual'}},...
			'loadVarsToRam', {{0,1}},...
			'guiEnabled', {{0,1}},...
			'dfofAnalysis', {{0,1}},...
			'picsSavePath', {{['private' filesep 'pics' filesep datestr(now,'yyyymmdd','local') filesep]}},...
			'delimiter', {{',','tab'}},...
			'imgSaveTypes', {{'-dpng','-dmeta','-depsc2'}}...
		);

		filterImageOptions = struct(...
			'minNumPixels', 10,...
			'maxNumPixels', 100,...
			'SNRthreshold', 1.45,...
			'minPerimeter', 5,...
			'maxPerimeter', 50,...
			'minSolidity', 0.8,...
			'minEquivDiameter', 3,...
			'maxEquivDiameter', 30,...
			'slopeRatioThreshold', 0.04...
		);

		downsampleRawOptions = struct(...
			'folderListInfo','A:\data\processing\',...
			'downsampleSaveFolder','B:\data\processing\',...
			'downsampleSrcFolder','E:\data\raw\',...
			'downsampleFactor','4',...
			'fileFilterRegexp','recording.*.hdf5',...
			'datasetName','/images',...
			'maxChunkSize','25000',...
			'srcFolderFilterRegexp','201\d',...
			'srcSubfolderFileFilterRegexp','recording.*.(txt|xml)',...
			'srcSubfolderFileFilterRegexpExt','(.txt|.xml)',...
			'downsampleSaveFolderTwo','E:\data\raw\',...
			'downsampleFactorTwo','2',...
			'outputDatasetName','/1'...
		);

		% io folders
		inputFolders = {};
		videoDir = '';
		videoSaveDir = '';
		trackingDir = '';
		stimulusDir = '';
		% if want to automatically save object to a specific location.
		objSaveLocation = [];

		% signal related
		% either the raw signals (traces) or
		rawSignals = {};
		%
		rawImages = {};
		% computed signal peaks/locations, to reduce computation in functions
		signalPeaks = {};
		%
		signalPeaksArray = {};
		% computed centroid locations {[x y],...}
		objLocations = {};
		% mean correlation coefficient between image and movie
		imageMovieCorr = {};
		% cellmaps indicating which cells were filtered
		rawImagesFiltered = {};
		% structure of classifier structures, each field in the property should be named after the folder's subject, e.g. classifierStructs.m667 is the classification structure for subject m667
		classifierStructs = {};
		% structure of classifier structures for each folder, e.g. after running classification on each
		classifierFolderStructs = {};
		% cell array with {signalNo}.signalFeatures, {signalNo}.imageFeatures
		classifierFeatures = {};
		% cell array with {signalNo}.signalFeatures, {signalNo}.imageFeatures
		classifierImageFeaturesNames = {'EquivDiameter','Area','Perimeter','Solidity'};
		% structure for all valid classifications to go
		valid = {};
		% Automated or manual classification
		validManual = {};
		% from automated classification
		validAuto = {};
		% valid cells based on a regional modification
		validRegionMod = {};
		% polygon vertices from previously selected regions
		validRegionModPoly = {};
		% whether or not rawSignals/rawImages have been replaced by only valid signals, hence, ignore validManual/validAuto
		validPurge = 0;
		% ROI to use for exclusion analysis
		analysisROIArray = {};
		% number of expected [PCs ICs] for PCA-ICA, alter for other procedures
		numExpectedSignals = {};

		% subject info
		% all are cell array of strings or numbers as specified in the name
		dataPath = {};
		subjectNum = {};
		subjectStr = {};
		assay = {};
		protocol = {};
		assayType = {};
		assayNum = {};
		imagingPlane = {};
		imagingPlaneNum = {};
		date = {};
		fileIDArray = {};
		fileIDNameArray = {};
		folderBaseSaveStr = {};
		folderBasePlaneSaveStr = {};
		folderBaseDisplayStr = {};

		% path to CSV/TAB file or matlab table containing trial information and frames when stimuli occur
		discreteStimulusTable = {};
		% cell array of strings
		stimulusNameArray = {};
		% cell array of strings, used for saving pictures/etc.
		stimulusSaveNameArray = {};
		% cell array of numbered values for stimulus, e.g. {65,10}
		stimulusIdArray = {};
		% [1 numTrialFrames] vectors with 1 for when stimulus occurs
		stimulusVectorArray = {};
		% vector sequence before/after to analyze stimulus
		stimulusTimeSeq = {};

		% path to a CSV/TAB file
		continuousStimulusTable = {};
		% cell array of strings
		continuousStimulusNameArray = {};
		% cell array of strings, used for saving pictures/etc.
		continuousStimulusSaveNameArray = {};
		% cell array of numbered values for stimulus, e.g. {65,10}
		continuousStimulusIdArray = {};
		% [1 numTrialFrames] vectors with 1 for when stimulus occurs
		continuousStimulusVectorArray = {};
		% vector sequence before/after to analyze stimulus
		continuousStimulusTimeSeq = {};

		% session trial names
		sessionTrialNames = '';

		% behavior metrics
		behaviorMetricTable = {};
		behaviorMetricNameArray = {};
		behaviorMetricIdArray = {};
	end
	properties(GetAccess = 'public', SetAccess = 'private')
		% public read access, but private write access.

		% summary statistics data save stores
		sumStats = {};
		detailStats = {};
		saveDataType = {};

		% counters and stores
		% index of current folder
		fileNum = 1;
		% same as fileNum, will transfer to this since more clear
		folderNum = 1;
		% number to current stimulus index
		stimNum = 1;
		figNames = {};
		figNo = {};
		figNoAll = 777;

		% signal related
		nSignals = {};
		nFrames = {};
		signalPeaksCopy = {};
		alignedSignalArray = {};
		alignedSignalShuffledMeanArray = {};
		alignedSignalShuffledStdArray = {};

		% stimulus
		% reorganize discreteStimulusTable into stimulus structures to reduce memory footprint
		discreteStimulusArray = {};
		discreteStimMetrics = {};

		% stimulus
		% reorganize discreteStimulusTable into stimulus structures to reduce memory footprint
		continuousStimulusArray = {};
		continuousStimMetrics = {};

		% behavior metric
		% reorganize behaviorMetricTable into stimulus structures to reduce memory footprint/provide common IO
		behaviorMetricArray = {};

		% distance metrics
		distanceMetric = {};
		distanceMetricShuffleMean = {};
		distanceMetricShuffleStd = {};

		% correlation metrics
		corrMatrix = {};

		% significant signals, different variables for controlling which signals are statistically significant, given some test
		currentSignificantArray = [];
		significantArray = {};
		sigModSignals = {};
		sigModSignalsAll = {};
		ttestSignSignals = {};

		% cross session alignment
		globalIDs = [];
		globalIDCoords = {};
		globalIDFolders = {};
		globalIDImages = {};
		globalRegistrationCoords = {};
		globalObjectMapTurboreg = [];
		globalStimMetric = [];
	end
	properties(GetAccess = 'private', SetAccess = 'private')
		% private read and write access
	end
	properties(Constant = true)
		% cannot be changed after object is created
		FRAMES_PER_SECOND =  5;
		DOWNSAMPLE_FACTOR =  4;
		MICRON_PER_PIXEL =  2.37;
	end

	methods
		% methods, including the constructor are defined in this block
		function obj = miniscopeAnalysis(varargin)
			% CLASS CONSTRUCTOR
			warning on;
% 			clc
			display([...
			'S-Lab Imaging Analysis Class v3.20170731' 10 ...
			'Biafra Ahanonu <<a href="emailto:bahanonu@gmail.com">bahanonu@gmail.com</a>>' 10 10 ...
			'Made in USA' 10 ...
			'* * * * * * * * * * =========================' 10 ...
			'* * * * * * * * * * :::::::::::::::::::::::::' 10 ...
			'* * * * * * * * * * =========================' 10 ...
			'* * * * * * * * * * :::::::::::::::::::::::::' 10 ...
			'* * * * * * * * * * =========================' 10 ...
			':::::::::::::::::::::::::::::::::::::::::::::' 10 ...
			'=============================================' 10 ...
			':::::::::::::::::::::::::::::::::::::::::::::' 10 ...
			'=============================================' 10 ...
			':::::::::::::::::::::::::::::::::::::::::::::' 10 ...
			'=============================================' 10 ...
			':::::::::::::::::::::::::::::::::::::::::::::' 10 ...
			'=============================================' 10])
			display(repmat('#',1,7))
			display('Constructing imaging analysis object...')

			% Because the obj
			%========================
			% obj.exampleOption = '';
			% get options
			obj = getOptions(obj,varargin);
			% display(options)
			% unpack options into current workspace
			% fn=fieldnames(options);
			% for i=1:length(fn)
			%	 eval([fn{i} '=options.' fn{i} ';']);
			% end
			%========================

			obj = initializeObj(obj);

			display('Done with setup!')
			display(repmat('#',1,7))

			display([...
			'Run processing pipeline by typing into command window:' 10 ...
			'<a href="">obj.runPipelineProcessing</a>' 10 ...
			'or for advanced features: ' 10 ...
			'<a href="">obj.runPipeline</a>' 10])
		end
		% getter and setter functions
		function dataPath = get.dataPath(obj)
			dataPath = obj.dataPath;
		end
	end
	methods(Static = true)
		% functions that are related but not dependent on instances of the class
		% function obj = loadObj(oldObj)
		%	  [filePath,folderPath,~] = uigetfile('*.*','select text file that points to analysis folders','example.txt');
		%	  % exit if user picks nothing
		%	  % if folderListInfo==0; return; end
		%	  load([folderPath filesep filePath]);
		%	  oldObj = obj;
		%	  obj = miniscopeAnalysis;
		%	  obj = getOptions(obj,oldObj);
		% end
	end
	methods(Access = private)
		% methods only executed by other class methods

		% model methods, usually for input-output like saving information to files
		% for obtaining the current stim from tables
		[behaviorMetric] = modelGetBehaviorMetric(obj,inputID)
	end
	methods(Access = protected)
		% methods only executed by other class methods, also available to subclasses
	end
	methods(Access = public)
		% these are in separate M-files

		[output] = modelGetStim(obj,idNum,varargin)

		% view help about the object
		obj = help(obj)

		% specific experiments
		obj = behaviorAntipsychotic(obj)
		obj = behaviorProtocolLoad(obj)

		% compute methods, performs some computation and returns calculation to class property
		obj = computeParameterSweepBehaviorNeuroMetrics(obj)
		obj = computeDiscreteAlignedSignal(obj,varargin)
		obj = computeSpatioTemporalClustMetric(obj)
		obj = computeSignalPeaksFxn(obj)
		obj = computeMatchObjBtwnTrials(obj)
		obj = computeAcrossTrialSignalStimMetric(obj)
		obj = computeSaveSignalPeak(obj)
		%
		obj = computeContinuousAlignedSignal(obj)
		%
		obj = computeClassifyTrainSignals(obj)
		obj = computeManualSortSignals(obj)
		% just need stimulus files
		obj = computePopulationDistance(obj,varargin)
		obj = computeDiscreteDimReduction(obj,varargin)
		obj = computeDiscreteStimulusDecoder(obj,varargin)
		obj = computeDiscreteRateStats(obj)
		obj = computeTrialSpecificActivity(obj)
		obj = modelEditStimTable(obj)

		% view methods, for displaying charts
		% no prior computation
		obj = viewSignalPeakDetection(obj)
		obj = viewStimTrigTraces(obj)
		obj = viewCorr(obj)
		obj = viewCreateObjmaps(obj)
		obj = viewContinuousSignalVideo(obj)

		% pre-processing checking
		obj = viewMovie(obj)
		obj = viewMovieFiltering(obj)
		obj = viewMovieRegistrationTest(obj)
		obj = modelModifyMovies(obj)

		% require pre-computation, individual
		obj = viewStimTrig(obj)
		obj = viewObjmapStimTrig(obj)
		obj = viewChartsPieStimTrig(obj)
		obj = viewObjmapSignificant(obj)
		obj = viewSpatioTemporalMetric(obj)
		% require pre-computation, group
		obj = viewPlotSignificantPairwise(obj)
		obj = viewObjmapSignificantPairwise(obj)
		obj = viewObjmapSignificantAllStims(obj)
		% require pre-computation, group, global alignment
		obj = viewMatchObjBtwnSessions(obj)
		% movies
		obj = viewMovieCreateSignalBasedStimTrig(obj)
		% require pre-computation and behavior metrics, individual
		obj = viewSignalBehaviorCompare(obj)
		% requires manual sorting
		obj = viewCompareSignalExtractionMethods(obj)
		% need tracking
		obj = viewOverlayTrackingToVideo(obj)

		% model methods, usually for input-output like saving information to files
		obj = modelReadTable(obj,varargin)
		obj = modelTableToStimArray(obj,varargin)
		obj = modelGetFileInfo(obj)
		obj = modelVerifyDataIntegrity(obj)
		obj = modelSaveImgToFile(obj,saveFile,thisFigName,thisFigNo,thisFileID)
		obj = modelSaveSummaryStats(obj)
		obj = modelSaveDetailedStats(obj)
		obj = modelVarsFromFiles(obj)
		obj = modelModifyRegionAnalysis(obj,varargin)
		obj = modelExtractSignalsFromMovie(obj)
		obj = modelAddNewFolders(obj)
		obj = modelPreprocessMovie(obj)
		obj = modelDownsampleRawMovies(obj)
		obj = modelBatchCopyFiles(obj)
		% helps clean and load tracking data
		obj = modelTrackingData(obj)


		% helper
		[inputSignals inputImages signalPeaks signalPeaksArray valid] = modelGetSignalsImages(obj,varargin)
		[fileIdxArray idNumIdxArray nFilesToAnalyze nFiles] = getAnalysisSubsetsToAnalyze(obj)
		[turboregSettingStruct] = getRegistrationSettings(obj,inputTitleStr)

		% set methods, for IO to specific variables in a controlled manner
		obj = setMainSettings(obj)
		obj = setStimulusSettings(obj)

		function obj = runPipelineProcessing(obj)
		  fxnsToRun = {...
		  '=======setup=======',
		  'modelAddNewFolders',
		  'saveObj',
		  'initializeObj',
		  'setMainSettings',
		  'setStimulusSettings',
		  'modelEditStimTable',
		  '=======preprocess=======',
		  'modelGetFileInfo',
		  'modelVerifyDataIntegrity',
		  'modelBatchCopyFiles',
		  '===',
		  'modelDownsampleRawMovies',
		  'viewMovieFiltering',
		  'viewMovieRegistrationTest',
		  'modelPreprocessMovie',
		  'modelModifyMovies',
		  'modelExtractSignalsFromMovie',
		  '===',
		  'modelVarsFromFiles',
		  '=======signal sorting=======',
		  'computeManualSortSignals',
		  'modelModifyRegionAnalysis',
		  'computeClassifyTrainSignals',
		  'viewCompareSignalExtractionMethods',
		  '=======preprocess verification=======',
		  'viewMovie',
		  'viewSignalPeakDetection',
		  'viewSubjectMovieFrames'
		  'viewMovieCreateSideBySide',
		  'viewCreateObjmaps',
		  'viewSignalStats',
		  'computeSaveSignalPeak',
		  '=======tracking=======',
		  'modelTrackingData',
		  'viewOverlayTrackingToVideo',
		  '=======across session analysis: compute/view=======',
		  'computeMatchObjBtwnTrials',
		  'viewMatchObjBtwnSessions',
		  '=======discrete analysis: save=======',
		  'modelSaveSummaryStats',
		  'modelSaveDetailedStats'
		  };
		  obj.runPipeline('fxnsToRun',fxnsToRun);
		end
		function obj = display(obj)
			obj.runPipeline;
			% display('hello');
		end
		function obj = showVars(obj)
			obj.disp;
		end
		function obj = runPipeline(obj,varargin)
			setFigureDefaults();
			set(0, 'DefaultUICOntrolFontSize', 14)
% 			close all;clc;

			fxnsToRun = {...
			'=======setup=======',
			'showVars',
			'modelAddNewFolders',
			'saveObj',
			'initializeObj',
			'setMainSettings',
			'setStimulusSettings',
			'modelEditStimTable',
			'behaviorProtocolLoad',
			'=======preprocess=======',
			'modelGetFileInfo',
			'modelVerifyDataIntegrity',
			'modelBatchCopyFiles',
			'===',
			'modelDownsampleRawMovies',
			'viewMovieFiltering',
			'viewMovieRegistrationTest',
			'modelPreprocessMovie',
			'modelModifyMovies',
			'modelExtractSignalsFromMovie',
			'===',
			'modelVarsFromFiles',
			'=======signal sorting=======',
			'computeManualSortSignals',
			'modelModifyRegionAnalysis',
			'computeClassifyTrainSignals',
			'viewCompareSignalExtractionMethods',
			'=======preprocess verification=======',
			'viewMovie',
			'viewSignalPeakDetection',
			'viewSubjectMovieFrames'
			'viewMovieCreateSideBySide',
			'viewMovieCreateSignalBasedStimTrig',
			'viewCreateObjmaps',
			'viewSignalStats',
			'computeSaveSignalPeak',
			'=======tracking=======',
			'modelTrackingData',
			'viewOverlayTrackingToVideo',
			'=======across session analysis: compute/view=======',
			'computeMatchObjBtwnTrials',
			'computeAcrossTrialSignalStimMetric',
			'viewMatchObjBtwnSessions',
			'=======discrete analysis: save=======',
			'modelSaveSummaryStats',
			'modelSaveDetailedStats'
			};
			%========================
			options.fxnsToRun = fxnsToRun;
			% get options
			options = getOptions(options,varargin);
			% display(options)
			% unpack options into current workspace
			% fn=fieldnames(options);
			% for i=1:length(fn)
			%	eval([fn{i} '=options.' fn{i} ';']);
			% end
			%========================
			fxnsToRun = options.fxnsToRun;
			% initialDir = pwd;
			% set back to initial directory in case exited early
			% restoredefaultpath;
			% loadRepoFunctions();
			if strcmp(obj.defaultObjDir,pwd)~=1
				cd(obj.defaultObjDir);
			end

			if ischar(obj.videoDir)
				obj.videoDir = {obj.videoDir};
			end

			% ensure private folders are set
			if ~exist(obj.picsSavePath,'dir');mkdir(obj.picsSavePath);end
			if ~exist(obj.dataSavePath,'dir');mkdir(obj.dataSavePath);end
			if ~exist(obj.logSavePath,'dir');mkdir(obj.logSavePath);end

			props = properties(obj);
			totSize = 0;
			% for ii=1:length(props)
			%	  currentProperty = getfield(obj, char(props(ii)));
			%	  s = whos('currentProperty');
			%	  totSize = totSize + s.bytes;
			% end
			% sprintf('%.f',totSize*1.0e-6)
			% fprintf(1, '%d bytes\n', totSize*1.0e-6);
			% runs all currently implemented view functions

			% turn off gui elements, run in batch
			obj.guiEnabled = 0;

			scnsize = get(0,'ScreenSize');
			dlgSize = [scnsize(3)*0.7 scnsize(4)*0.8];

			currentIdx = find(strcmp(fxnsToRun,obj.currentMethod));
			[idNumIdxArray, ok] = listdlg('ListString',fxnsToRun,'InitialValue',currentIdx(1),'ListSize',dlgSize,'Name','Sir! I have a plan!');

			if ok==0
				return
			end

			[guiIdx, ok] = listdlg('ListString',{'Yes','No'},'InitialValue',1,'ListSize',dlgSize,'Name','Gui Enabled?');
			% idNumIdxArray
			% turn off gui elements, run in batch
			obj.guiEnabled = guiIdx==1;

			excludeList = {'showVars','setMainSettings','modelAddNewFolders','saveObj','setStimulusSettings','modelDownsampleRawMovies'};

			fxnsToRun = {fxnsToRun{idNumIdxArray}};
			obj.currentMethod = fxnsToRun{1};

			if isempty(intersect(fxnsToRun,excludeList))
				scnsize = get(0,'ScreenSize');
				usrIdxChoiceStr = {'PCAICA','EM','EXTRACT','CNMF','ROI'};
				usrIdxChoiceDisplay = {'PCAICA (Mukamel, 2009)','CELLMax (Lacey)','EXTRACT (Hakan)','CNMF (Pnevmatikakis, 2015)','ROI'};
				% use current string as default
				currentIdx = find(strcmp(usrIdxChoiceStr,obj.signalExtractionMethod));
				[sel, ok] = listdlg('ListString',usrIdxChoiceDisplay,'InitialValue',currentIdx,'ListSize',dlgSize,'Name','Cell extraction algorithm to use for analysis');
				% (Americans love a winner)
				usrIdxChoiceList = {2,1};
				obj.signalExtractionMethod = usrIdxChoiceStr{sel};
			end
			if ~isempty(obj.inputFolders)&isempty(intersect(fxnsToRun,excludeList))
				if isempty(obj.protocol)
					obj.modelGetFileInfo();
				end
				folderNumList = strsplit(num2str(1:length(obj.inputFolders)),' ');
				selectList = strcat(folderNumList(:),'/',num2str(length(obj.inputFolders)),' | ',obj.date(:),' _ ',obj.protocol(:),' _ ',obj.fileIDArray(:),' | ',obj.inputFolders(:));
				% set(0, 'DefaultUICOntrolFontSize', 16)
				% select subjects to analyze
				subjectStrUnique = unique(obj.subjectStr);
				[subjIdxArray, ok] = listdlg('ListString',subjectStrUnique,'ListSize',dlgSize,'Name','which subjects to analyze?');
				subjToAnalyze = subjectStrUnique(subjIdxArray);
				subjToAnalyze = find(ismember(obj.subjectStr,subjToAnalyze));
				% get assays to analyze
				assayStrUnique = unique(obj.assay(subjToAnalyze));
				[assayIdxArray, ok] = listdlg('ListString',assayStrUnique,'ListSize',dlgSize,'Name','which assays to analyze?');
				assayToAnalyze = assayStrUnique(assayIdxArray);
				assayToAnalyze = find(ismember(obj.assay,assayToAnalyze));
				% filter for folders chosen by the user
				validFoldersIdx = intersect(subjToAnalyze,assayToAnalyze);
				% if isempty(validFoldersIdx)
				%	  continue;
				% end
				useAltValid = {'no additional filter','manual index entry','manually sorted folders','not manually sorted folders','manual classification already in obj',['has ' obj.signalExtractionMethod ' extracted cells'],['missing ' obj.signalExtractionMethod ' extracted cells'],'fileFilterRegexp','valid auto'};
				useAltValidStr = {'no additional filter','manual index entry','manually sorted folders','not manually sorted folders','manual classification already in obj',['has extracted cells'],'missing extracted cells','fileFilterRegexp','valid auto'};
				[choiceIdx, ok] = listdlg('ListString',useAltValid,'ListSize',dlgSize,'Name','Choose additional folder sorting filters');
				if ok==1
					useAltValid = useAltValidStr{choiceIdx};
				else
					useAltValid = 0;
				end
				% useAltValid = 0;
				switch useAltValid
					case 'manual index entry'
					 theseSettings = inputdlg({...
							 'list (separated by commas) of indexes'
						 },...
						 'Folders to process',1,...
						 {...
							 '1'
						 }...
					 );
					 validFoldersIdx = str2num(theseSettings{1});
					case 'missing extracted cells'
						switch obj.signalExtractionMethod
							case 'PCAICA'
								missingRegexp = obj.rawICfiltersSaveStr;
							case 'EM'
								missingRegexp = obj.rawEMStructSaveStr;
							case 'EXTRACT'
								missingRegexp = obj.rawEXTRACTStructSaveStr;
							otherwise
								missingRegexp = obj.rawICfiltersSaveStr;
						end
						validFoldersIdx2 = [];
						for folderNo = 1:length(obj.dataPath)
							filesToLoad = getFileList(obj.dataPath{folderNo},missingRegexp);
							if isempty(filesToLoad)
								display(['no extracted signals: ' obj.dataPath{folderNo}])
								validFoldersIdx2(end+1) = folderNo;
							end
						end
						validFoldersIdx = intersect(validFoldersIdx,validFoldersIdx2)
					case 'has extracted cells'
						switch obj.signalExtractionMethod
							case 'PCAICA'
								cellRegexp = obj.rawICfiltersSaveStr;
							case 'EM'
								cellRegexp = obj.rawEMStructSaveStr;
							case 'EXTRACT'
								cellRegexp = obj.rawEXTRACTStructSaveStr;
							otherwise
								cellRegexp = obj.rawICfiltersSaveStr;
						end
						validFoldersIdx2 = [];
						for folderNo = 1:length(obj.dataPath)
							filesToLoad = getFileList(obj.dataPath{folderNo},cellRegexp);
							if ~isempty(filesToLoad)
								display(['has extracted signals: ' obj.dataPath{folderNo}])
								validFoldersIdx2(end+1) = folderNo;
							end
						end
						validFoldersIdx = intersect(validFoldersIdx,validFoldersIdx2)
					case 'fileFilterRegexp'
						validFoldersIdx2 = [];
						for folderNo = 1:length(obj.dataPath)
							filesToLoad = getFileList(obj.dataPath{folderNo},obj.fileFilterRegexp);
							if isempty(filesToLoad)
								validFoldersIdx2(end+1) = folderNo;
								display(['missing dfof: ' obj.dataPath{folderNo}])
							end
						end
						validFoldersIdx = intersect(validFoldersIdx,validFoldersIdx2)
					case 'valid auto'
						validFoldersIdx = find(cell2mat(cellfun(@isempty,obj.validAuto,'UniformOutput',0)));
					case 'not manually sorted folders'
						switch obj.signalExtractionMethod
							case 'PCAICA'
								missingRegexp = obj.sortedICdecisionsSaveStr;
							case 'EM'
								missingRegexp = obj.sortedEMStructSaveStr;
							case 'EXTRACT'
								missingRegexp = obj.sortedEXTRACTStructSaveStr;
							otherwise
								missingRegexp = obj.sortedICdecisionsSaveStr;
						end
						validFoldersIdx = [];
						display(['missingRegexp: ' missingRegexp])
						for folderNo = 1:length(obj.inputFolders)
							filesToLoad = getFileList(obj.inputFolders{folderNo},missingRegexp);
							% filesToLoad
							%filesToLoad
							if isempty(filesToLoad)
								validFoldersIdx(end+1) = folderNo;
								display(['not manually sorted: ' obj.dataPath{folderNo}])
							else
								display(['manually sorted: ' obj.dataPath{folderNo}])
							end
						end
					case 'manually sorted folders'
						switch obj.signalExtractionMethod
							case 'PCAICA'
								missingRegexp = obj.sortedICdecisionsSaveStr;
							case 'EM'
								missingRegexp = obj.sortedEMStructSaveStr;
							case 'EXTRACT'
								missingRegexp = obj.sortedEXTRACTStructSaveStr;
							otherwise
								missingRegexp = obj.sortedICdecisionsSaveStr;
						end
						validFoldersIdx = [];
						missingRegexp
						for folderNo = 1:length(obj.inputFolders)
							filesToLoad = getFileList(obj.inputFolders{folderNo},missingRegexp);
							%filesToLoad
							if ~isempty(filesToLoad)
								validFoldersIdx(end+1) = folderNo;
								display(['manually sorted: ' obj.dataPath{folderNo}])
							end
						end
					case 'manual classification already in obj'
						validFoldersIdx = find(arrayfun(@(x) ~isempty(x{1}),obj.validManual));
					otherwise
						% body
				end
				[fileIdxArray, ok] = listdlg('ListString',selectList,'ListSize',dlgSize,'Name','which folders to analyze?','InitialValue',validFoldersIdx);
				if ok==0
					return
				end
				obj.foldersToAnalyze = fileIdxArray;
				if isempty(obj.stimulusNameArray)
					obj.discreteStimuliToAnalyze = [];
				else
					[idNumIdxArray, ok] = listdlg('ListString',obj.stimulusNameArray,'ListSize',dlgSize,'Name','which stimuli to analyze?');
					if ok==0
						return
					end
					obj.discreteStimuliToAnalyze = idNumIdxArray;
				end
			end
			for thisFxn=fxnsToRun
				try
					display(repmat('!',1,21))
					display(['Running: obj.' thisFxn{1}]);
					obj.(thisFxn{1});
				catch err
					display(repmat('@',1,7))
					disp(getReport(err,'extended','hyperlinks','on'));
					display(repmat('@',1,7))
					if strcmp(obj.defaultObjDir,pwd)~=1
						restoredefaultpath;
						cd(obj.defaultObjDir);
						loadRepoFunctions();
					end
				end
			end
			obj.guiEnabled = 1;
			obj.foldersToAnalyze = [];
			% set back to initial directory in case exited early
			% restoredefaultpath;
			% loadRepoFunctions();
			if strcmp(obj.defaultObjDir,pwd)~=1
				cd(obj.defaultObjDir);
			end
		end

		function GetSize(obj)
			props = properties(obj);
			totSize = 0;
			for ii=1:length(props)
				currentProperty = getfield(obj, char(props(ii)));
				s = whos('currentProperty');
				totSize = totSize + s.bytes;
			end
			fprintf(1, '%d bytes\n', totSize);
		end
		% save the current object instance
		function obj = saveObj(obj)

			if isempty(obj.objSaveLocation)
				[filePath,folderPath,~] = uiputfile('*.*','select folder to save object mat file to','miniscopeAnalysis_properties.mat');
				% exit if user picks nothing
				% if folderListInfo==0; return; end
				savePath = [folderPath filesep filePath];
				% tmpObj = obj;
				% obj = struct(obj);
				obj.objSaveLocation = savePath;
			else
				savePath = obj.objSaveLocation;
			end
			display(['saving to: ' savePath])
			try
			  save(savePath,'obj','-v7.3');
			catch
			  display('Problem saving, choose new location...')
			  obj.objSaveLocation = [];
			  obj.saveObj();
			end
			% obj = tmpObj;
		end

		function obj = initializeObj(obj)
			% load dependencies.
			loadRepoFunctions();
			% ensure private folders are set
			if ~exist(obj.picsSavePath,'dir');mkdir(obj.picsSavePath);end
			if ~exist(obj.dataSavePath,'dir');mkdir(obj.dataSavePath);end
			if ~exist(obj.logSavePath,'dir');mkdir(obj.logSavePath);end

			% if use puts in a single folder or a path to a txt file with folders
			if ~isempty(obj.rawSignals)&strcmp(class(obj.rawSignals),'char')
				if isempty(regexp(obj.rawSignals,'.txt'))&exist(obj.rawSignals,'dir')==7
					% user just inputs a single directory
					obj.rawSignals = {obj.rawSignals};
				else
					% user input a file linking to directories
					fid = fopen(obj.rawSignals, 'r');
					tmpData = textscan(fid,'%s','Delimiter','\n');
					obj.rawSignals = tmpData{1,1};
					fclose(fid);
				end
				obj.inputFolders = obj.rawSignals;
				obj.dataPath = obj.rawSignals;
			end
			% add subject information to object given datapath
			if ~isempty(obj.dataPath)
				obj.modelGetFileInfo();
			else
				display('No folder paths input, run <a href="">modelAddNewFolders</a> method.');
				% warning('Input data paths for all files!!! option: dataPath')
			end
			if ~isempty(obj.discreteStimulusTable)&~strcmp(class(obj.discreteStimulusTable),'table')
				obj.modelReadTable('table','discreteStimulusTable');
				obj.modelTableToStimArray('table','discreteStimulusTable','tableArray','discreteStimulusArray','nameArray','stimulusNameArray','idArray','stimulusIdArray','valueName',obj.stimulusTableValueName,'frameName',obj.stimulusTableFrameName);
			end
			if ~isempty(obj.continuousStimulusTable)&~strcmp(class(obj.continuousStimulusTable),'table')
				obj.delimiter = ',';
				obj.modelReadTable('table','continuousStimulusTable','addFileInfoToTable',1);
				obj.delimiter = ',';
				obj.modelTableToStimArray('table','continuousStimulusTable','tableArray','continuousStimulusArray','nameArray','continuousStimulusNameArray','idArray','continuousStimulusIdArray','valueName',obj.stimulusTableValueName,'frameName',obj.stimulusTableFrameName,'grabStimulusColumnFromTable',1);
			end
			% load behavior tables
			if ~isempty(obj.behaviorMetricTable)&~strcmp(class(obj.behaviorMetricTable),'table')
				obj.modelReadTable('table','behaviorMetricTable');
				obj.modelTableToStimArray('table','behaviorMetricTable','tableArray','behaviorMetricArray','nameArray','behaviorMetricNameArray','idArray','behaviorMetricIdArray','valueName','value');
			end
			% modify stimulus naming scheme
			if ~isempty(obj.stimulusNameArray)
				obj.stimulusSaveNameArray = obj.stimulusNameArray;
				obj.stimulusNameArray = strrep(obj.stimulusNameArray,'_',' ');
			end
			% load all the data
			if ~isempty(obj.rawSignals)&strcmp(class(obj.rawSignals{1}),'char')
				display('paths input, going to load files')
				obj.guiEnabled = 0;
				obj = modelVarsFromFiles(obj);
				obj.guiEnabled = 1;
			end
			% check if signal peaks have already been calculated
			if isempty(obj.signalPeaks)&~isempty(obj.rawSignals)
				% obj.computeSignalPeaksFxn();
			else
				display('No folder data specified, load data with <a href="">modelVarsFromFiles</a> method.');
				% warning('no signal data input!!!')
			end
			% load stimulus tables
		end
	end
end