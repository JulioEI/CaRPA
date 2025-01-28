function [outStruct] = classifySignals(inputImages,inputSignals,varargin)
    % runs a specified model on the data and attempts to train a classifier (svm, nnet, etc.) to it.
    % biafra ahanonu
    % started: 2013.08.09
    % adapted from biafra ahanonu nn_cell_classification.m
    % inputs
        % inputImages - cell array of [x y nSignals] matrices containing each set of images corresponding to inputSignals objects.
        % inputSignals - cell array of [nSignals frames] matrices containing each set of inputImages signals
    % outputs
        % outStruct.inputFeatures
        % outStruct.classifications - [1 nSignals] vector containing classification from an algorithm or the consensus classification (mode) of all algorithms if using all_classify as input.
    % options
        % inputStruct - can input previously output structure (outStruct), will automatically search for classifier and output structure won't overwrite the old one.
        % classifierType - char string | 'nnet' 'svm' 'glm'
        % trainingOrClassify - char string | 'train' or 'classify'
        % classifier - struct | input a previously trained classifier or will automatically grab from outStruct.svmClassifier, outStruct.nnetClassifier, or outStruct.glmCoeffs depending on the algorithm.
        % inputTargets - cell | cell array of [1 nSignals] vectors with 1/0 classification of good/bad
        % featureList - list of image features to use from regionprops
        % additionalFeatureList - string cell array of names for additional features to input to the classifier.
        % additionalFeatures - numeric matrix [nSignals nFeatures], nFeatures should equal length(additionalFeatureList)
    % example
        % % train dataset
        % [classifyStruct] = classifySignals({inputImages},{inputSignals},'inputTargets',valid,'classifierType','glm');
        % [classifyStruct] = classifySignals({inputImages},{inputSignals},'classifierType','all','trainingOrClassify','training','inputMovieList',{'pathToMovieFile'},'inputTargets',{inputTargets});
        % % classify dataset
        % [classifyStruct] = classifySignals({inputImages},{inputSignals},'inputStruct',classifyStruct,'classifierType','nnet','trainingOrClassify','classify');
        % [classifyStruct] = classifySignals({inputImages},{inputSignals},'classifierType','all','trainingOrClassify','classify','inputMovieList',{'pathToMovieFile'},'inputTargets',{inputTargets},'inputStruct',classifyStruct);

    % changelog
        % 2014.02.02 - generalized nnet classifier to include other classification schemes, finished svm classifier
        % 2014.02.10 - added glm classifier, included image feature list option, allow input of outStruct as inputStruct.
        % 2014.05.02 - fixed svm classify, needed to transpose input features matrix
        % 2014.06.20 - improved nnet, set NaN in input features to zero to avoid problems with classification.
        % 2015.12.22 - Added additional comments, slight restructuring.
        % 2016.01.02 - Integrating movie image correlation into function directly from computeClassifyTrainSignals
        % 2017.01.14 [20:06:04] - support switched from [nSignals x y] to [x y nSignals]
    % TODO
        % Add ability to iterate over different parameter space for nnet, svm, and glm to find the optimal set.

    %========================
    % input a previously created structure to prevent overwriting.
    options.inputStruct = struct;
    % 'nnet' 'svm' 'glm' 'all'
    options.classifierType = 'svm';
    % 'training' or 'classify'
    options.trainingOrClassify = 'training';
    % previously trained classifier
    options.classifier = [];
    % known targets to classify with (or compare classifier results to)
    options.inputTargets = [];
    % input a list of movies to use for image correlation metrics
    options.inputMovieList = {};
    % hierarchy name in hdf5 where movie is
    options.inputDatasetName = '/1';
    % list of image features to calculate using regionprops
    % options.featureList = {'Eccentricity','EquivDiameter','Area','Orientation','Perimeter','Solidity'};
    options.featureList = {'EquivDiameter','Area','Perimeter','Solidity'};
    % pre-computed image features
    options.inputImageFeatures = {};
    % full feature list
    options.fullFeatureList = ['meanPeakSnr',options.featureList,'nPeaks3Std','nPeaks1Std','meanPeakSlopeRatio','meanPeakFwhm','meanPeakAmplitude'];
    % input additional feature, should match additionalFeatures nFeatures
    options.additionalFeatureList = {};
    % matrix [nSignals nFeatures], numeric
    options.additionalFeatures = [];
    % get options
    options = getOptions(options,varargin,'recursiveStructs',0);
    display(options)
    % unpack options into current workspace
    % fn=fieldnames(options);
    % for i=1:length(fn)
    %   eval([fn{i} '=options.' fn{i} ';']);
    % end
    %========================
    try
        % obtain features from each set of images
        numInputs = length(inputImages);
        if iscell(inputImages)&isempty(inputImages{1})
            display('EMPTY INPUT TO CLASSIFIER! STOPPING.')
            return
        end
        if ~isempty(options.inputMovieList)
            options.fullFeatureList{end+1} = 'imageMovieCorr';
        end
        % get features
        for inputNo = 1:numInputs
            display(['getting features for ' num2str(inputNo) '/' num2str(numInputs)])
            if length(options.inputImageFeatures)>=inputNo&~isempty(options.inputImageFeatures{inputNo})
                display('using pre-computed image features...')
                inputFeatures{inputNo} = computeFeatures(inputImages{inputNo},inputSignals{inputNo},options.inputImageFeatures{inputNo},options,inputNo);
            else
                inputFeatures{inputNo} = computeFeatures(inputImages{inputNo},inputSignals{inputNo},[],options,inputNo);
            end
            if ~isempty(options.additionalFeatures)
                tmpAdditionalFeatures = options.additionalFeatures{inputNo};
                inputFeatures{inputNo} = horzcat(inputFeatures{inputNo},tmpAdditionalFeatures(:));
            end
        end
        % inputFeatures
        inputFeatures = cat(1,inputFeatures{:})';
        outStruct.inputFeatures = inputFeatures;

        % if input targets present, plot information on distribution of features based on input classifications
        if ~isempty(options.inputTargets)
            [options outStruct] = plotTargetVsFeatures(options,outStruct);
        end

        % replace output structure with user's input structure
        if ~isempty(options.inputStruct)
            outStruct = options.inputStruct;
        else
            outStruct.null = nan;
        end

        % add targets to output structure or get them from inputted structure
        % if strcmp('training',options.trainingOrClassify)
        %     outStruct.inputTargets = options.inputTargets;
        % elseif strcmp('classify',options.trainingOrClassify)
        %     if ~any(strcmp('inputTargets',fieldnames(outStruct)))
        %         options.inputTargets = outStruct.inputTargets;
        %     end
        % end

        %========================
        % start either training or classification of input signals/images
        classificationOption = strcat(options.classifierType,'_',options.trainingOrClassify);
        display(classificationOption)
        [figHandle figNo] = openFigure(564, '');
        switch classificationOption
            %========================
            case 'all_training'
                % train all three classifiers at once
                [outStruct.svmClassifier] = svmTrainerFxn(inputFeatures,options.inputTargets(:));
                [outStruct.nnetClassifier] = nnetTrainerFxn(inputFeatures,options.inputTargets(:));
                [outStruct.glmCoeffs] = glmTrainerFxn(inputFeatures,options.inputTargets(:));
                % break;
            case 'all_classify'
                % classify using each algorithm
                [outStruct.svmGroups] = svmClassifyFxn(outStruct.svmClassifier,inputFeatures);
                [outStruct.nnetGroups] = nnetClassifyFxn(outStruct.nnetClassifier,inputFeatures);
                [outStruct.glmGroups] = glmClassifyFxn(outStruct.glmCoeffs,inputFeatures);
                % get a consensus score
                svmClassification = logical(outStruct.svmGroups>0.5);
                nnetClassification = logical(outStruct.nnetGroups>0.5);
                glmClassification = logical(outStruct.glmGroups>0.5);
                outStruct.classifications = [svmClassification(:) nnetClassification(:) glmClassification(:)];
                outStruct.classifications = mode(outStruct.classifications,2);
                classifierGroups = {outStruct.svmGroups,outStruct.nnetGroups,outStruct.glmGroups,outStruct.classifications};
                classifierNames = {'svm','nnet','glm','consensus - svm nnet glm'};
                for classifierNo = 1:length(classifierGroups)
                    if ~isempty(options.inputTargets)
                        warning off;
                        figure(classifierNo);
                        plotconfusion(options.inputTargets(:)',classifierGroups{classifierNo}(:)')
                        title(classifierNames{classifierNo})
                        [c,cm,ind,per] = confusion(options.inputTargets(:)',classifierGroups{classifierNo}(:)');
                        % outStruct.confusionPct{classifierNo} = mean(per,1);
                        perTab = crosstab(options.inputTargets(:),logical(classifierGroups{classifierNo}(:)>0.5));
                        % outStruct.confusionPct{classifierNo} = perTab/length(options.inputTargets(:));
                        outStruct.confusionPct{classifierNo} = perTab(:);
                        outStruct.confusionStr{classifierNo} = ['TN FP FN TP'];
                        warning on;
                    end
                end
            case 'svm_training'
                [outStruct.svmClassifier] = svmTrainerFxn(inputFeatures,options.inputTargets(:));
                % break;
            case 'svm_classify'
                if isempty(options.classifier)
                    options.classifier = outStruct.svmClassifier;
                end
                [outStruct.svmGroups] = svmClassifyFxn(options.classifier,inputFeatures);
                outStruct.classifications = outStruct.svmGroups;
                % plot confusion matrix
                if ~isempty(options.inputTargets)
                    plotconfusion(options.inputTargets(:)',outStruct.svmGroups(:)')
                    [c,cm,ind,per] = confusion(options.inputTargets(:)',outStruct.svmGroups(:)');
                    outStruct.confusionPct = mean(per,1);
                end
                % break;
            %========================
            case 'nnet_training'
                [outStruct.nnetClassifier] = nnetTrainerFxn(inputFeatures,options.inputTargets(:));
                % break;
            case 'nnet_classify'
                if isempty(options.classifier)
                    options.classifier = outStruct.nnetClassifier
                end
               [outStruct.nnetGroups] = nnetClassifyFxn(options.classifier,inputFeatures);
               outStruct.classifications = outStruct.nnetGroups;
               % size(outStruct.nnetGroups)
               % figure(99920)
               % plot(options.inputTargets(:)'+2,'k')
               % hold on
               % plot(outStruct.nnetGroups(:)','r')
               % legend({'targets','output'})
               % hold off
               % pause
               % plot confusion matrix
               if ~isempty(options.inputTargets)
                   plotconfusion(options.inputTargets(:)',outStruct.nnetGroups(:)')
                   [c,cm,ind,per] = confusion(options.inputTargets(:)',outStruct.nnetGroups(:)');
                   outStruct.confusionPct = mean(per,1);
                   mean(per,1)
               end
               % break;
            %========================
            case 'glm_training'
                [outStruct.glmCoeffs] = glmTrainerFxn(inputFeatures,options.inputTargets(:));
                % break;
            case 'glm_classify'
                if isempty(options.classifier)
                    options.classifier = outStruct.glmCoeffs;
                end
                [outStruct.glmGroups] = glmClassifyFxn(options.classifier,inputFeatures);
                outStruct.classifications = outStruct.glmGroups;
                % plot confusion matrix
                if ~isempty(options.inputTargets)
                    plotconfusion(options.inputTargets(:)',outStruct.glmGroups(:)')
                    [c,cm,ind,per] = confusion(options.inputTargets(:)',outStruct.glmGroups(:)');
                    outStruct.confusionPct = mean(per,1);
                end
                % break;
            %========================
            otherwise
                display('invalid choice specified');
                return
        end
        outStruct.error = 0;
        outStruct.options = options;
    catch err
        outStruct.error = 1;
        display(repmat('@',1,7))
        disp(getReport(err,'extended','hyperlinks','on'));
        display(repmat('@',1,7))
    end
end
function [options outStruct] = plotTargetVsFeatures(options,outStruct)
    options.inputTargets = cat(2,options.inputTargets{:});
    outStruct.inputTargets = options.inputTargets;
    inputFeatures = outStruct.inputFeatures;

    % valid = options.inputTargets;
    pointColors = ['g','r'];
    options.fullFeatureList = ['meanPeakSnr',options.featureList,'nPeaks3Std','nPeaks1Std','meanPeakSlopeRatio','meanPeakFwhm','meanPeakAmplitude'];
    if ~isempty(options.inputMovieList)
        options.fullFeatureList{end+1} = 'imageMovieCorr';
    end
    nameList = options.fullFeatureList;
    if ~isempty(options.additionalFeatureList)
        nameList = [nameList options.additionalFeatureList];
    end
    [figHandle figNo] = openFigure(756, '');
        clf
    [figHandle figNo] = openFigure(757, '');
        clf
        warning off
        for pointNum = 1:2
            pointColor = pointColors(pointNum);
            if pointNum==1
                valid = logical(options.inputTargets);
            else
                valid = logical(~options.inputTargets);
            end
            for i=1:length(nameList)
                [figHandle figNo] = openFigure(756, '');
                    subplot(3,ceil(length(nameList)/3),i)
                        % eval(['iStat=imgStats.' nameList{i} ';']);
                        iStat = inputFeatures(i,:);
                        plot(find(valid),iStat(valid),[pointColor '.'])
                        % title(nameList{i})
                        hold on;box off;
                        xlabel('rank'); ylabel(nameList{i})
                        xlim([1 length(iStat)]);
                        % hold off
                [figHandle figNo] = openFigure(757, '');
                    subplot(3,ceil(length(nameList)/3),i)
                        % eval(['iStat=imgStats.' nameList{i} ';']);
                        iStat = inputFeatures(i,:);
                        % plot(rand([1 length(find(valid))])/3+pointNum,iStat(valid),[pointColor '.'])
                        iStat = iStat(valid);
                        switch pointNum
                            case 1
                                distributionPlot(iStat(:),'histOri','left','color','g','widthDiv',[2 1],'showMM',0,'histOpt',0)
                            case 2
                                distributionPlot(iStat(:),'histOri','right','color','r','widthDiv',[2 2],'showMM',0,'histOpt',0)
                            otherwise
                                % body
                        end
                        title(nameList{i})
                        hold on;box off;
                        % xlabel('rank'); ylabel(nameList{i})
                        % xlim([1 length(iStat)]);
                        % hold off
            end
        end
        warning on
    drawnow
end
function [inputFeatures] = computeFeatures(inputImages,inputSignals,inputImageFeatures,options,inputNo)
    % obtains the training features using several sub-routines

    % number of 3 std peaks
    [signalPeaks, signalPeaksArray] = computeSignalPeaks(inputSignals,'numStdsForThresh',3);
    numOfPeaks3std = sum(signalPeaks,2);
    % get the SNR for traces
    [signalSnr ~] = computeSignalSnr(inputSignals,'testpeaks',signalPeaks,'testpeaksArray',signalPeaksArray);
    % get the peak statistics
    [peakOutputStat] = computePeakStatistics(inputSignals,'waitbarOn',1,'testpeaks',signalPeaks,'testpeaksArray',signalPeaksArray);
    slopeRatio = peakOutputStat.slopeRatio;
    avgFwhm = peakOutputStat.avgFwhm;
    avgPeakAmplitude = peakOutputStat.avgPeakAmplitude;
    % number of 1 std peaks
    [signalPeaks, signalPeaksArray] = computeSignalPeaks(inputSignals,'numStdsForThresh',1);
    numOfPeaks1std = sum(signalPeaks,2);
    % Eccentricity,EquivDiameter,Area,Orientation,Perimeter,Solidity,
    if isempty(inputImageFeatures)
        [imgStats] = computeImageFeatures(inputImages, 'thresholdImages',1,'featureList',options.featureList);
        imgFeatures = cell2mat(struct2cell(imgStats))';
    else
        imgFeatures = inputImageFeatures;
    end

    % get image correlation metrics
    if ~isempty(options.inputMovieList)
        [inputMovie] = loadMovieList(options.inputMovieList{inputNo},'inputDatasetName',options.inputDatasetName);
        [outputImages outputMeanImageCorrs] = createPeakTriggeredImages(inputMovie, inputImages, inputSignals,'signalPeaksArray',signalPeaksArray,'normalizeOutput',0);
        outputMeanImageCorrs(isnan(outputMeanImageCorrs)) = 0;
    else
        outputMeanImageCorrs = [];
    end

    % concatenate all the features together
    % inputFeatures = horzcat(signalSnr(:),imgFeatures,numOfPeaks(:),slopeRatio(:),avgFwhm(:),avgPeakAmplitude(:));
    % nameList = {'SNR','Eccentricity','EquivDiameter','Area','Orientation','Perimeter','Solidity','nPeaks','slopeRatio','avgFwhm','avgPeakAmplitude'};
    if isempty(outputMeanImageCorrs)
        inputFeatures = horzcat(signalSnr(:),imgFeatures,numOfPeaks3std(:),numOfPeaks1std(:),slopeRatio(:),avgFwhm(:),avgPeakAmplitude(:));
    else
        inputFeatures = horzcat(signalSnr(:),imgFeatures,numOfPeaks3std(:),numOfPeaks1std(:),slopeRatio(:),avgFwhm(:),avgPeakAmplitude(:),outputMeanImageCorrs(:));
    end

    % figure(900)
    % plot(inputFeatures)
    % legend({'signalSnr','sizeImageObj','numOfPeaks','slopeRatio','avgFwhm','avgPeakAmplitude'});
    % pause
end
function [svmClassifier] = svmTrainerFxn(inputFeatures,inputTargets)
    % train SVM
    optionsStruct = statset('Display','iter');
    size(inputFeatures)
    size(inputTargets(:))
    inputFeatures(isnan(inputFeatures)) = 0;
    svmClassifier = svmtrain(inputFeatures,inputTargets,'options',optionsStruct);
end
function [svmGroups] = svmClassifyFxn(svmClassifier,inputFeatures)
    % classify SVM based on svmtrain support vectors
    size(inputFeatures')
    size(svmClassifier.SupportVectors)
    inputFeatures(isnan(inputFeatures)) = 0;
    svmGroups = svmclassify(svmClassifier,inputFeatures');%'showplot',true
end
function [nnetClassifier] = nnetTrainerFxn(inputFeatures,inputTargets)
    % trains a neural network classifier
    % create a network
    nnetClassifier = feedforwardnet([100 50 25 10 8 4]);
    nnetClassifier.trainFcn = 'trainscg';
    % manageParallelWorkers('parallel',options.parallel);
    % train network
    % myCluster = parcluster('local');
    % myCluster.NumWorkers = 6;
    % saveProfile(myCluster);
    % matlabpool open
    %net = train(net,trainingSet,targets,'useGPU','yes');
    % [nnetClassifier tr] = train(nnetClassifier,inputFeatures,inputTargets,'useParallel','yes','showResources','yes');
    inputFeatures(isnan(inputFeatures)) = 0;
    [nnetClassifier tr] = train(nnetClassifier,inputFeatures,inputTargets(:)','showResources','yes');
    % matlabpool close

    % ploterrhist(targets - net(trainingSet))
end
function [nnetGroups] = nnetClassifyFxn(nnetClassifier,inputFeatures)
    % scores inputs based on train and feedforwardnet
    inputFeatures(isnan(inputFeatures)) = 0;
    nnetGroups = sim(nnetClassifier,inputFeatures);
end
function [glmCoeffs] = glmTrainerFxn(inputFeatures,inputTargets)
    % trains a glm using a normal distribution
    inputFeatures(isnan(inputFeatures)) = 0;
    glmCoeffs = glmfit(inputFeatures',inputTargets,'normal');
end
function [glmGroups] = glmClassifyFxn(glmCoeffs,inputFeatures)
    % classifies features into groups based on coefficients from glmfit
    inputFeatures(isnan(inputFeatures)) = 0;
    glmGroups = glmval(glmCoeffs,inputFeatures','identity');
end