function obj = computeClassifyTrainSignals(obj)
	% compute peaks for all signals if not already input
	% biafra ahanonu
	% branched from controllerAnalysis: 2014.08.01 [16:09:16]
	% inputs
		%
	% outputs
		%

	% changelog
		% 2017.01.14 [20:06:04] - support switched from [nSignals x y] to [x y nSignals]
	% TODO
		% Add in option for different types of classification

	% display('NOT FULLY CONVERTED TO CLASS METHOD YET...')
	% return

	scnsize = get(0,'ScreenSize');
	classifyOrTrain = {'training','classify','classify_output_table','classify to valid'};
	'export classifier parameters'
	[fileIdxArray, ok] = listdlg('ListString',classifyOrTrain,'ListSize',[scnsize(3)*0.2 scnsize(4)*0.25],'Name','train or classify?');
	classifyOrTrain = classifyOrTrain(fileIdxArray);

	trainingSet = {'all','using only 0%, 50%, and 100% trials','use first 3 trials','use first 1 trials','manual'};
	[fileIdxArray, ok] = listdlg('ListString',trainingSet,'ListSize',[scnsize(3)*0.2 scnsize(4)*0.25],'Name','train or classify?');
	trainingSet = trainingSet{fileIdxArray};

    [fileIdxArray idNumIdxArray nFilesToAnalyze nFiles] = obj.getAnalysisSubsetsToAnalyze();

    movieSettings = inputdlg({...
            'processed imaging movie regexp:',...
        },...
        'Settings for HDF5 data checking',[1 100],...
        {...
            obj.fileFilterRegexp,...
        }...
    );
    fileFilterRegexp = movieSettings{1};obj.fileFilterRegexp = fileFilterRegexp;
    for classifyOrTrainStr = classifyOrTrain
    	classifyOrTrainStr = classifyOrTrainStr{1};
	    switch classifyOrTrainStr
	    	case 'classify to valid'
			    saveAllClassifiedSignals = {'Yes','No'};
			    [tmpArray, ok] = listdlg('ListString',saveAllClassifiedSignals,'ListSize',[scnsize(3)*0.2 scnsize(4)*0.25],'Name','Save classifications to file?');
			    saveAllClassifiedSignals = tmpArray==1;
	    		% body
	    	otherwise
	    		% body
	    end
	end

	switch trainingSet
		case 'manual'
			theseSettings = inputdlg({...
					'list (separated by commas) of indexes'
				},...
					'Folders to process',1,...
				{...
					'1'
				}...
			);
			tmpFolderIdx = str2num(theseSettings{1});
		otherwise
			% do nothing
	end

	display(repmat('#',1,21))
	display('computing signal peaks...')
	nFiles = length(obj.rawSignals);
	subjectList = unique(obj.subjectStr(fileIdxArray));
	for thisSubjectStr=subjectList
		display(repmat('=',1,21))
		thisSubjectStr = thisSubjectStr{1};
		display(thisSubjectStr);
		useImageMovieCorr = 1;
		for classifyOrTrainStr = classifyOrTrain
			classifyOrTrainStr = classifyOrTrainStr{1};
			switch classifyOrTrainStr
				case 'training'
					try
						% get folders with correct subject and that have already been manually classified
						validFoldersIdx = find(strcmp(thisSubjectStr,obj.subjectStr));
						% fileIdxArray
						validFoldersIdx = intersect(validFoldersIdx,fileIdxArray);
						validManualIdx = find(arrayfun(@(x) ~isempty(x{1}),obj.validManual));
						trainingFoldersIdx = intersect(validFoldersIdx,validManualIdx);
						% use only three trials
						switch trainingSet
							case 'all'
								% body
							case 'using only 0%, 50%, and 100% trials'
								display('using only 0%, 50%, and 100% trials')
								trainingFoldersIdx = trainingFoldersIdx(floor(quantile(1:length(trainingFoldersIdx),[0 0.5 1])));
								trainingFoldersIdx = unique(trainingFoldersIdx);
							case 'use first 3 trials'
								trainingFoldersIdx = trainingFoldersIdx(1:3);
							case 'use first 1 trials'
								trainingFoldersIdx = trainingFoldersIdx(1);
							case 'manual'
								trainingFoldersIdx = trainingFoldersIdx(tmpFolderIdx);
							otherwise
								% body
						end
						%
						ioption.classifierType = 'all';
						ioption.trainingOrClassify = 'training';
						ioption.inputTargets = {obj.validManual{trainingFoldersIdx}};
						% for idx = 1:nTrainingFolders
						trainingRawImages = {};
						trainingRawSignals = {};
						additionalFeatures = {};
						obj.folderBaseSaveStr{trainingFoldersIdx}
						nTrainingFolders = length(trainingFoldersIdx);
	                    nTrainingFolders
						for idx = 1:nTrainingFolders
							ioption.inputTargets{idx} = logical(ioption.inputTargets{idx}==1);
							ioption.inputTargets{idx}(isnan(ioption.inputTargets{idx})) = 0;
							folderNo = trainingFoldersIdx(idx);
							obj.fileNum = folderNo;
							display(repmat('-',1,21))
							display([num2str(idx) '/' num2str(nTrainingFolders) ' (' num2str(obj.fileNum) '/' num2str(nFiles) '): ' obj.fileIDNameArray{obj.fileNum}]);
							[rawSignals, rawImages, signalPeaks, signalPeaksArray] = modelGetSignalsImages(obj,'returnType','raw');

							% get image corr with movie if available
							if useImageMovieCorr==1
								try imageMovieCorr=obj.imageMovieCorr{obj.fileNum}.(obj.signalExtractionMethod); catch; imageMovieCorr=[]; end
									% [imageMovieCorr useImageMovieCorr additionalFeatures additionalFeatureList]
								[imageMovieCorr, useImageMovieCorr, additionalFeatures, additionalFeatureList] = getImageCorr(obj,idx,'train',obj.inputFolders{obj.fileNum},obj.fileFilterRegexp,imageMovieCorr,rawImages,rawSignals,signalPeaks,signalPeaksArray,additionalFeatures);
								% [imageMovieCorr ioption{folderNo} useImageMovieCorr] = getImageCorr(trainOrClassify,obj.inputFolders{obj.fileNum},obj.fileFilterRegexp,imageMovieCorr);
								% ,ioption{folderNo}
								obj.imageMovieCorr{obj.fileNum}.(obj.signalExtractionMethod) = imageMovieCorr;

								% getImageCorr('train');
							end

							trainingRawImages{idx} = rawImages;
							trainingRawSignals{idx} = rawSignals;

							try obj.classifierFeatures{obj.fileNum}.(obj.signalExtractionMethod).imageFeatures;imageFeaturesCheck=1; catch; imageFeaturesCheck=0; end
							if imageFeaturesCheck==1
								display('using pre-computed image features...')
								ioption.inputImageFeatures{idx} = obj.classifierFeatures{obj.fileNum}.(obj.signalExtractionMethod).imageFeatures;
							else
								ioption.inputImageFeatures{idx} = [];
							end
						end
						ioption.additionalFeatureList = additionalFeatureList;
						ioption.additionalFeatures = additionalFeatures;
						ioption.featureList = obj.classifierImageFeaturesNames;

						% classify signals
						[outputStruct] = classifySignals(trainingRawImages,trainingRawSignals,'options',ioption);
						obj.classifierStructs.(thisSubjectStr) = outputStruct;
						figure(756)
						set(756,'PaperUnits','inches','PaperPosition',[0 0 16 9])
						obj.modelSaveImgToFile([],'classifyFeatures_','current',[]);

						figure(757)
						set(757,'PaperUnits','inches','PaperPosition',[0 0 16 9])
						obj.modelSaveImgToFile([],'classifyFeaturesSide_','current',[]);
					catch err
						display(repmat('@',1,7))
						disp(getReport(err,'extended','hyperlinks','on'));
						display(repmat('@',1,7))
					end
				case 'classify'
					validFoldersIdx = find(strcmp(thisSubjectStr,obj.subjectStr));
                    validFoldersIdx = intersect(validFoldersIdx,fileIdxArray);
                    if isempty(validFoldersIdx)
                        continue
                    end
					nValidFolders = length(validFoldersIdx);

					% ========================
					manageParallelWorkers('parallel',1,'maxCores',7,'openCloseParallelPool','open');
					% ========================
					imageMovieCorrAll = obj.imageMovieCorr;
					signalExtractionMethod = obj.signalExtractionMethod;
					validManual = obj.validManual;
					% get subject classifier structure
					classifierStruct = obj.classifierStructs.(thisSubjectStr);
                    imgCorrSwitch = sum(strcmp('imageMovieCorr',classifierStruct.options.additionalFeatureList)>0);
					% folderBaseSaveStr = obj.folderBaseSaveStr(obj.fileNum);
					folderBaseSaveStr = obj.folderBaseSaveStr;
					for idx = 1:nValidFolders
						try
							obj.fileNum = validFoldersIdx(idx);
							folderNo = validFoldersIdx(idx);
							display(repmat('-',1,21))
							display([num2str(idx) '/' num2str(nValidFolders) ' (' num2str(folderNo) '/' num2str(nFiles) '): ' obj.fileIDNameArray{obj.fileNum}]);

							display(folderBaseSaveStr{folderNo})
							[rawSignals{folderNo}, rawImages{folderNo}, signalPeaks{folderNo}, signalPeaksArray{folderNo}] = modelGetSignalsImages(obj,'returnType','raw');
						catch err
							display(repmat('@',1,7))
							disp(getReport(err,'extended','hyperlinks','on'));
							display(repmat('@',1,7))
						end
                    end

                    additionalFeatures = {};
					for idx = 1:nValidFolders
						try
							obj.fileNum = validFoldersIdx(idx);
							folderNo = validFoldersIdx(idx);
							display(repmat('-',1,21))
							display([num2str(idx) '/' num2str(nValidFolders) ' (' num2str(folderNo) '/' num2str(nFiles) '): ' obj.fileIDNameArray{obj.fileNum}]);

							display(folderBaseSaveStr{folderNo})
							% [rawSignals, rawImages, signalPeaks, signalPeaksArray] = modelGetSignalsImages(obj,'returnType','raw');
							if isempty(rawSignals{folderNo})
								continue
							end
							% ioption{folderNo}.classifierType = 'all';
							% ioption{folderNo}.trainingOrClassify = 'classify';
							if ~isempty(validManual{folderNo})
								inputTargets = {validManual{folderNo}};
								inputTargets{1} = logical(inputTargets{1}==1);
								inputTargets{1}(isnan(inputTargets{1})) = 0;
							else
								inputTargets = {};
							end
							% ioption{folderNo}.inputStruct = classifierStruct;

							% get image corr with movie if available
							if imgCorrSwitch>0
                                % z = imageMovieCorrAll{folderNo};
                                %z.(signalExtractionMethod)
								try imageMovieCorr=imageMovieCorrAll{folderNo}.(signalExtractionMethod); catch; imageMovieCorr=[]; end
								[imageMovieCorr, ~,  additionalFeatures, additionalFeatureList] = getImageCorr(obj,idx,'classify',obj.inputFolders{obj.fileNum},obj.fileFilterRegexp,imageMovieCorr,rawImages{folderNo},rawSignals{folderNo},signalPeaks{folderNo},signalPeaksArray{folderNo},additionalFeatures)
								imageMovieCorrAll{folderNo}.(signalExtractionMethod) = imageMovieCorr;

								% getImageCorr('classify');
                            else
								additionalFeatures = [];
								additionalFeatureList = {};
							end

							[obj.classifierFolderStructs{obj.fileNum}] = classifySignals({rawImages{folderNo}},{rawSignals{folderNo}},...
								'classifierType','all',...
								'trainingOrClassify','classify',...
								'inputTargets',inputTargets,...
								'inputStruct',classifierStruct,...
								'additionalFeatures',additionalFeatures,...
								'additionalFeatureList',additionalFeatureList);
							% clear ioption

							plotValid(obj)
						catch err
							display(repmat('@',1,7))
							disp(getReport(err,'extended','hyperlinks','on'));
							display(repmat('@',1,7))
						end
					end
					% obj.classifierFolderStructs{obj.fileNum}
					for idx = 1:nValidFolders
						try
							obj.fileNum = validFoldersIdx(idx);
							if ~isempty(obj.validManual{obj.fileNum})
								% if ~any(strcmp('summaryStats',fieldnames(ostruct)))
								if ~exist('summaryStats','var')
									summaryStats.subject{1,1} = nan;
									summaryStats.assay{1,1} = nan;
									summaryStats.assayType{1,1} = nan;
									summaryStats.assayNum{1,1} = nan;
									summaryStats.imagingPlane{1,1} = nan;
									summaryStats.FN{1,1} = nan;
									summaryStats.FP{1,1} = nan;
									summaryStats.TP{1,1} = nan;
									summaryStats.TN{1,1} = nan;
									summaryStats.classifierType{1,1} = nan;
								end
								classifierTypes = {'svm','nnet','glm','svm_nnet_glm'};
								% obj.classifierFolderStructs{obj.fileNum}.confusionPct{classifierNo}(1);
								for classifierNo = 1:length(classifierTypes)
									% summaryStats.confusionPctFN{end+1,1} = obj.classifierFolderStructs{obj.fileNum}.confusionPct{classifierNo}(1);
									% summaryStats.confusionPctFP{end+1,1} = obj.classifierFolderStructs{obj.fileNum}.confusionPct{classifierNo}(2);
									% summaryStats.confusionPctTP{end+1,1} = obj.classifierFolderStructs{obj.fileNum}.confusionPct{classifierNo}(3);
									% summaryStats.confusionPctTN{end+1,1} = obj.classifierFolderStructs{obj.fileNum}.confusionPct{classifierNo}(4);

									summaryStats.FN{end+1,1} = obj.classifierFolderStructs{obj.fileNum}.confusionPct{classifierNo}(3);
									summaryStats.FP{end+1,1} = obj.classifierFolderStructs{obj.fileNum}.confusionPct{classifierNo}(2);
									summaryStats.TP{end+1,1} = obj.classifierFolderStructs{obj.fileNum}.confusionPct{classifierNo}(4);
									summaryStats.TN{end+1,1} = obj.classifierFolderStructs{obj.fileNum}.confusionPct{classifierNo}(1);
									summaryStats.classifierType{end+1,1} = classifierTypes{classifierNo};
									summaryStats.subject{end+1,1} = obj.subjectStr{obj.fileNum};
									summaryStats.assay{end+1,1} = obj.assay{obj.fileNum};
									summaryStats.assayType{end+1,1} = obj.assayType{obj.fileNum};
									summaryStats.assayNum{end+1,1} = obj.assayNum{obj.fileNum};
									summaryStats.imagingPlane{end+1,1} = obj.imagingPlane{obj.fileNum};
								end
							end

							if exist('summaryStats','var')
								obj.sumStats = summaryStats;
								% write out summary statistics
							    savePath = [obj.dataSavePath obj.protocol{obj.fileNum} '_classifySummary.tab'];
							    display(['saving data to: ' savePath])
								writetable(struct2table(obj.sumStats),savePath,'FileType','text','Delimiter','\t');
							end
						catch err
							display(repmat('@',1,7))
							disp(getReport(err,'extended','hyperlinks','on'));
							display(repmat('@',1,7))
						end
					end
				case 'classify_output_table'
					validFoldersIdx = find(strcmp(thisSubjectStr,obj.subjectStr));
	                validFoldersIdx = intersect(validFoldersIdx,fileIdxArray);
	                nValidFolders = length(validFoldersIdx);
	                validManual = obj.validManual;
					for idx = 1:nValidFolders
						obj.fileNum = validFoldersIdx(idx);
						folderNo = validFoldersIdx(idx);
						display(repmat('-',1,21))
							display([num2str(idx) '/' num2str(nValidFolders) ' (' num2str(obj.fileNum) '/' num2str(nFiles) '): ' obj.folderBaseSaveStr{obj.fileNum}]);
						if ~isempty(validManual{folderNo})
							inputTargets = validManual{folderNo};
						else
							continue
						end
						classifierTypes = {'svm','nnet','glm','svm_nnet_glm'};
						% obj.classifierFolderStructs{obj.fileNum}.confusionPct{classifierNo}(1);
						if ~exist('summaryStats','var')
							summaryStats.subject{1,1} = nan;
							summaryStats.assay{1,1} = nan;
							summaryStats.assayType{1,1} = nan;
							summaryStats.assayNum{1,1} = nan;
							summaryStats.imagingPlane{1,1} = nan;
							summaryStats.FN{1,1} = nan;
							summaryStats.FP{1,1} = nan;
							summaryStats.TP{1,1} = nan;
							summaryStats.TN{1,1} = nan;
							summaryStats.classifierType{1,1} = nan;
						end
						for classifierNo = 1:length(classifierTypes)
							switch classifierTypes{classifierNo}
								case 'svm'
									classifierPredictions = obj.classifierFolderStructs{obj.fileNum}.svmGroups;
								case 'nnet'
									classifierPredictions = obj.classifierFolderStructs{obj.fileNum}.nnetGroups;
								case 'glm'
									classifierPredictions = obj.classifierFolderStructs{obj.fileNum}.glmGroups;
								case 'svm_nnet_glm'
									classifierPredictions = obj.classifierFolderStructs{obj.fileNum}.classifications;
								otherwise
									% body
							end

							perTab = crosstab(logical(inputTargets(:)),logical(classifierPredictions(:)>0.5));
							perTab = perTab(:);
							summaryStats.FN{end+1,1} = perTab(3);
							summaryStats.FP{end+1,1} = perTab(2);
							summaryStats.TP{end+1,1} = perTab(4);
							summaryStats.TN{end+1,1} = perTab(1);
							summaryStats.classifierType{end+1,1} = classifierTypes{classifierNo};
							summaryStats.subject{end+1,1} = obj.subjectStr{obj.fileNum};
							summaryStats.assay{end+1,1} = obj.assay{obj.fileNum};
							summaryStats.assayType{end+1,1} = obj.assayType{obj.fileNum};
							summaryStats.assayNum{end+1,1} = obj.assayNum{obj.fileNum};
							summaryStats.imagingPlane{end+1,1} = obj.imagingPlane{obj.fileNum};
						end
					end
					if exist('summaryStats','var')
						obj.sumStats = summaryStats;
						% write out summary statistics
					    savePath = [obj.dataSavePath obj.protocol{obj.fileNum} '_classifySummary.tab'];
					    display(['saving data to: ' savePath])
						writetable(struct2table(obj.sumStats),savePath,'FileType','text','Delimiter','\t');
					end

				case 'classify to valid'
					validFoldersIdx = find(strcmp(thisSubjectStr,obj.subjectStr));
                    validFoldersIdx = intersect(validFoldersIdx,fileIdxArray);
                    nValidFolders = length(validFoldersIdx);
					for idx = 1:nValidFolders
						try
							obj.fileNum = validFoldersIdx(idx);
							display(repmat('-',1,21))
								display([num2str(idx) '/' num2str(nValidFolders) ' (' num2str(obj.fileNum) '/' num2str(nFiles) '): ' obj.folderBaseSaveStr{obj.fileNum}]);
							% obj.folderBaseSaveStr{obj.fileNum}
							thisStruct = obj.classifierFolderStructs{obj.fileNum};
							fprintf('Adding to: obj.valid{%d}.(%s).classifier\n',obj.fileNum,obj.signalExtractionMethod);
							obj.valid{obj.fileNum}.(obj.signalExtractionMethod).classifier = thisStruct.classifications>0.5;
                        catch err
							display(repmat('@',1,7))
							disp(getReport(err,'extended','hyperlinks','on'));
							display(repmat('@',1,7))
                        end
					end

					% saveAllClassifiedSignals = 1;
					if saveAllClassifiedSignals==1
						% [fileIdxArray idNumIdxArray nFilesToAnalyze nFiles] = obj.getAnalysisSubsetsToAnalyze();
						for thisFileNumIdx = 1:nValidFolders
							try
								% thisFileNum = fileIdxArray(thisFileNumIdx);
								thisFileNum = validFoldersIdx(idx);
								obj.fileNum = thisFileNum;
								display(repmat('=',1,21))
								display([num2str(thisFileNumIdx) '/' num2str(nFilesToAnalyze) ' (' num2str(thisFileNum) '/' num2str(nFiles) '): ' obj.fileIDNameArray{obj.fileNum}]);


								currentFolderSaveStr = [obj.inputFolders{obj.fileNum} filesep obj.folderBaseSaveStr{obj.fileNum}];
								validClassifier = obj.valid{obj.fileNum}.(obj.signalExtractionMethod).classifier;
								if ~isempty(validClassifier)
									switch obj.signalExtractionMethod
										case 'PCAICA'
											saveID = {obj.classifierICdecisionsSaveStr}
										case 'EM'
											saveID = {obj.classifierEMStructSaveStr}
										otherwise
												% body
									end
									saveID = {obj.classifierICdecisionsSaveStr}
									saveVariable = {'validClassifier'}
									for i=1:length(saveID)
										savestring = [currentFolderSaveStr saveID{i}];
										display(['saving: ' savestring])
										save(savestring,saveVariable{i});
									end
								end
	                        catch err
								display(repmat('@',1,7))
								disp(getReport(err,'extended','hyperlinks','on'));
								display(repmat('@',1,7))
	                        end
						end
						% return
					end
				case 'compare to manual'


				otherwise
					% body
			end
		end
	end
	if exist('summaryStats','var')
		obj.sumStats = summaryStats;
		% write out summary statistics
	    savePath = [obj.dataSavePath obj.protocol{obj.fileNum} '_classifySummary.tab'];
	    display(['saving data to: ' savePath])
		writetable(struct2table(obj.sumStats),savePath,'FileType','text','Delimiter','\t');
	end
end

function [imageMovieCorr, useImageMovieCorr, additionalFeatures, additionalFeatureList] = getImageCorr(obj,idx,trainOrClassify,inputFolders,fileFilterRegexp,imageMovieCorr,rawImages,rawSignals,signalPeaks,signalPeaksArray,additionalFeatures)

	imageMovieCorr = imageMovieCorr;
	% get image correlation from movies
	movieList = getFileList(inputFolders, fileFilterRegexp);
	if isempty(movieList)
		movieList = getFileList(inputFolders, '.h5');
		if isempty(movieList)
			movieList = getFileList(inputFolders, '.h5');
		end
	end
	cellfun(@disp,movieList);
	useImageMovieCorr = 1;
	if ~isempty(movieList)
		if ~isempty(imageMovieCorr)
			outputMeanImageCorrs = obj.imageMovieCorr{obj.fileNum}.(obj.signalExtractionMethod);
			display('using pre-computed image-movie correlations')
		else
			% movieList = getFileList(inputFolders, obj.fileFilterRegexp);
			[inputMovie] = loadMovieList(movieList);
			[outputImages outputMeanImageCorrs] = createPeakTriggeredImages(inputMovie, rawImages, rawSignals,'signalPeaksArray',signalPeaksArray,'normalizeOutput',0);
			if isempty(outputMeanImageCorrs)
				display('no longer using image corr for this subject')
				useImageMovieCorr = 0;
				return;
			end
			outputMeanImageCorrs(isnan(outputMeanImageCorrs)) = 0;
			rawImages = outputImages;
			obj.imageMovieCorr{obj.fileNum}.(obj.signalExtractionMethod) = outputMeanImageCorrs;

		end
		additionalFeatureList = {'imageMovieCorr'};
		switch trainOrClassify
			case 'train'
				additionalFeatures{idx} = [outputMeanImageCorrs(:)];
				% body
			case 'classify'
				additionalFeatures{1} = [outputMeanImageCorrs(:)];
			otherwise
				% body
		end
	else
		display('no longer using image corr for this subject')
		useImageMovieCorr = 0;
		additionalFeatureList = {};
		additionalFeatures = {};
	end
end
function plotValid(obj)
	valid = obj.classifierFolderStructs{obj.fileNum}.classifications;
	% originalValid = valid;
	validNorm = normalizeVector(valid,'normRange','oneToOne');
	validDiff = [0 diff(valid')];
	%
	[figHandle figNo] = openFigure(10000, '');
	clf
	plot(valid);hold on;
	plot(validDiff,'g');
	%
	% validQuantiles = quantile(valid,[0.4 0.3]);
	% validHigh = validQuantiles(1);
	% validLow = validQuantiles(2);
	validHigh = 0.7;
	validLow = 0.5;
	%
	valid(valid>=validHigh) = 1;
	valid(valid<=validLow) = 0;
	valid(isnan(valid)) = 0;
	% questionable classification
	valid(validDiff<-0.3) = 2;
	valid(valid<validHigh&valid>validLow) = 2;
	%
	plot(valid,'r');
	plot(validNorm,'k');box off;
	legend({'scores','diff(scores)','classification','normalized scores'})
end
	% if strcmp('classify',options.trainingOrClassify)
	% 	ioption.classifierType = options.classifierType;
	% 	ioption.trainingOrClassify = options.trainingOrClassify;
	% 	ioption.inputTargets = {ostruct.validArray{ostruct.counter}};
	% 	ioption.inputStruct = classifierStruct
	% 	[ostruct.classifier] = classifySignals({ostruct.inputImages{ostruct.counter}},{ostruct.inputSignals{ostruct.counter}},'options',ioption);
	% 	% ostruct.data.confusionPct
	% 	% ostruct.classifier.confusionPct
	% 	if ~any(strcmp('summaryStats',fieldnames(ostruct)))
	% 		ostruct.summaryStats.subject{1,1} = nan;
	% 		ostruct.summaryStats.assay{1,1} = nan;
	% 		ostruct.summaryStats.assayType{1,1} = nan;
	% 		ostruct.summaryStats.assayNum{1,1} = nan;
	% 		ostruct.summaryStats.confusionPctFN{1,1} = nan;
	% 		ostruct.summaryStats.confusionPctFP{1,1} = nan;
	% 		ostruct.summaryStats.confusionPctTP{1,1} = nan;
	% 		ostruct.summaryStats.confusionPctTN{1,1} = nan;
	% 	end
	% 	ostruct.summaryStats.subject{end+1,1} = ostruct.subject{fileNum};
	% 	ostruct.summaryStats.assay{end+1,1} = ostruct.assay{fileNum};
	% 	ostruct.summaryStats.assayType{end+1,1} = ostruct.info.assayType{fileNum};
	% 	ostruct.summaryStats.assayNum{end+1,1} = ostruct.info.assayNum{fileNum};
	% 	ostruct.summaryStats.confusionPctFN{end+1,1} = ostruct.classifier.confusionPct(1);
	% 	ostruct.summaryStats.confusionPctFP{end+1,1} = ostruct.classifier.confusionPct(2);
	% 	ostruct.summaryStats.confusionPctTP{end+1,1} = ostruct.classifier.confusionPct(3);
	% 	ostruct.summaryStats.confusionPctTN{end+1,1} = ostruct.classifier.confusionPct(4);
	% end