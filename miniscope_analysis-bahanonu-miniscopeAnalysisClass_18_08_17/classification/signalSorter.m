function [inputImages inputSignals choices] = signalSorter(inputImages,inputSignals,varargin)
    % displays a GUI for sorting images and their associated signals, also does preliminary sorting based on image/signal properties
    % biafra ahanonu
    % started: 2013.10.08
    % dependent code
        % getOptions.m, createObjMap.m, removeSmallICs.m, identifySpikes.m
    % inputs
        % inputImages - [N x y] matrix where N = number of images, x/y are dimensions. Use permute(inputImages,[3 1 2]) if you use [x y N] for matrix indexing.
        % inputSignals - [N time] matrix where N = number of signals (traces) and time = frames.
        % inputID - obsolete, kept for compatibility, just input empty []
        % nSignals - obsolete, kept for compatibility
    % outputs
        % inputImages - [N x y] matrix where N = number of images, x/y are dimensions with only manual choices kept.
        % inputSignals
        % choices

    % changelog
        % 2013.10.xx changed to ginput and altered UI to show more relevant information, now shows a objMap overlayed with the current filter, etc.
        % 2013.11.01 [15:48:56]
            % Finished removing all cell array indexing by day, increase maintainability.
            % Input is now filters and traces instead of loading a directory inside fxn (which is cryptic). Output is filtered traces.
            % Can now move forward AND back, 21st century stuff. Also changed some of the other controls to make UI more user friendly.
        % 2013.11.03 [12:45:03] added a panel so that you can see the average trace around all spikes in an IC filter's trace along with several other improvements.
        % 2013.11.04 [10:30:40] changed invalid subscripting to valid, previous way involved negating choices, prone to error.
        % 2013.11.13 [09:25:24] added the ability to loop around and pre-maturely exit
        % 2013.11.19 [09:19:07] auto-saves decisions in case of a crash or other problem
        % 2013.12.07 [16:30:32] added more option (e.g. 's' key to mark rest of signals as bad)
        % 2013.12.10 [09:38:57] refactored a bit to make code more clear
        % 2013.12.15 [22:48:56] now overlays the good and bad images onto the entire image cell map, good for determining whether you've hit all the 'relevant' images
        % 2014.01.05 [09:23:54] small amount of
        % 2014.01.27 - started better integration auto-detecting based on SNR, etc.
        % 2014.03.06 - integrated support for manual scoring of automatic classification via abstraction (not explicitly loading classifier, but scoring pre-defined questionable input signals)
        % 2014.03.12 - sort by SNR or random, view montage of movie frames at peak or compare the signal to the movie directly
        % 2014.05.19 - improved SNR sort for NaNs, montage handles traces with no peaks, etc.
        % 2015.11.22 - refactored to make createStimCutMovieMontage faster
        % 2016.08.06 [20:51:02] - updated to make obj cut movie montage the default instead of static images, more useful.

    % TODO
        % DONE: allow option to mark rest as bad signals
        % c should make an obj cut movie for 10 or so signals
        % set viewMontageso it uses minValTraces maxValTraces

    % ============================
    % set default options
    options.nSignals = size(inputImages,3);
    % string to display over the cell map
    options.inputStr = '';
    % can pre-load choices, 1 = good, 0 = bad, 2 = questionable
    options.valid = [];
    % directory to store temporary decisions
    options.tmpDir = ['private' filesep 'tmp'];
    % id for the current session, use system time since it'll be unique
    options.sessionID = num2str(java.lang.System.currentTimeMillis);
    % threshold for SNR auto-annotate
    options.SnrThreshold = 1.2;
    % TBD
    options.slopeRatioThreshold = 0;
    % location of classifier
    options.classifierFilepath = [];
    % type of classifier that was used
    options.classifierType = 'nnet';
    % upper range pct score to manually sort
    options.upperClassifierThres = 0.6;
    % lower range pct score to manually sort
    options.lowerClassifierThres = 0.3;
    % movie matching inputImages/inputSignals src, used to find movie frames at peaks
    options.inputMovie = [];
    % sort by the SNR
    options.sortBySNR = 0;
    % randomize order
    options.randomizeOrder = 0;
    % show ROI traces in addition to input traces
    options.showROITrace = 0;
    % pre-compute signal peaks
    options.signalPeaks = [];
    options.signalPeaksArray = [];
    % ROI for peak signal plotting
    options.peakROI = [-20:20];
    % movie stats
    options.movieMax = NaN;
    options.movieMin = NaN;
    options.movieMean = NaN;
    % min value to display
    options.minValConstant = -0.1;
    % max movie cut images to show
    options.maxSignalsToShow = 24;
    % size in pixels to crop around the movie
    options.cropSize = 5;
    options.cropSizeLength = 10;
    % what percent of movie max to have crosshair values
    options.crossHairPercent = 0.07;
    % threshold for thresholding images
    options.threshold = 0.3;
    % threshold for thresholding images
    options.thresholdOutline = 0.3;
    % enter in rgb (range 0 to 1) for background for each color
    options.backgroundGood = [208,229,180]/255;
    options.backgroundBad = [244,166,166]/255;
    options.backgroundNeutral = repmat(230,[1 3])/255;

    options.colormap = customColormap([]);
    % get options
    options = getOptions(options,varargin);
    % unpack options into current workspace
    fn=fieldnames(options);
    for i=1:length(fn)
        eval([fn{i} '=options.' fn{i} ';']);
    end
    % ============================
    options

    for figNoFake = [1996 1997 1776 1777 1778 1779 1]
        [~, ~] = openFigure(figNoFake, '');
        clf
    end

    if ~isempty(options.inputMovie)
        display('calculating movie max...')
        options.movieMax = nanmax(options.inputMovie(:));
        display('calculating movie min...')
        options.movieMin = nanmin(options.inputMovie(:));
        display('calculating movie mean...')
        try
            options.movieMean = nanmean(options.inputMovie(1:10,1:10,:));
        catch
            options.movieMean = nanmean(options.inputMovie(:));
        end
    end

    % for manual classification of automated signals
    if ~isempty(valid)&~isempty(find(valid==2))
        inputImagesBackup = inputImages;
        inputSignalsBackup = inputSignals;
        questionableSignalIdx = find(valid==2);
        inputImages = inputImages(:,:,questionableSignalIdx);
        inputSignals = inputSignals(questionableSignalIdx,:);
        validBackup = valid;
        valid = zeros(1,length(questionableSignalIdx));
    else
        validBackup = [];
    end


    % get the SNR for traces and sort traces by this if asked
    [signalSnr ~] = computeSignalSnr(inputSignals,'testpeaks',options.signalPeaks,'testpeaksArray',options.signalPeaksArray);
    if options.sortBySNR==1
        signalSnr(isnan(signalSnr)) = -Inf;
        [signalSnr newIdx] = sort(signalSnr,'descend');
        signalSnr(isinf(signalSnr)) = NaN;
        inputSignals = inputSignals(newIdx,:);
        inputImages = inputImages(:,:,newIdx);
        if ~isempty(valid)
            valid = valid(newIdx);
        end
    end

    % randomize the order if asked
    if options.randomizeOrder==1
        randIdx = randperm(options.nSignals);
        inputSignals = inputSignals(randIdx,:);
        inputImages = inputImages(:,:,randIdx);
        if ~isempty(valid)
            valid = valid(randIdx);
        end
    end
    % =======
    % create a cell map to overlay current IC filter onto
    objMap = createObjMap(inputImages);

    % =======
    % Pre-compute values
    % threshold images
    inputImagesThres = thresholdImages(inputImages,'waitbarOn',1,'binary',0,'normalizationType','zeroToOne','threshold',options.threshold);
    inputImagesThresBinary = inputImagesThres>0;
    % inputImagesThresBinary = thresholdImages(inputImages,'waitbarOn',1,'binary',1,'threshold',options.threshold);

    [imgStats] = computeImageFeatures(inputImagesThres,'thresholdImages',1);

    if isempty(options.signalPeaks)
        [signalPeaks, signalPeakIdx] = computeSignalPeaks(inputSignals,'makePlots', 0,'makeSummaryPlots',0,'waitbarOn',1);
    else
        signalPeaks = options.signalPeaks;
        signalPeakIdx = options.signalPeaksArray;
    end
    % imagesc(options.signalPeaks)
    % add max peak for those signals that don't otherwise have any
    % peaksNoneIdx = sum(signalPeaks,2)>0;
    % peaksNoneIdx = find(~nPeaksAll);
    minimumEvents = 4;
    % figure(223123)
    % plot(cellfun(@length,signalPeakIdx))
    peaksNoneIdx = cellfun(@length,signalPeakIdx)>minimumEvents;
    peaksNoneIdx = find(~peaksNoneIdx);
    % peaksNoneIdx
    if ~isempty(peaksNoneIdx)
        fprintf('Adding pseudo-peaks for %d low event signals...',length(peaksNoneIdx))
        for idxHere = peaksNoneIdx
            tmpS = inputSignals(idxHere,:);
            % plot(tmpS)
            % drawnow
            % size(tmpS)
            % [yy ii] = max(tmpS(:));
            [sortedX,sortingIndices] = sort(tmpS(:),'descend');
            % maxValues = sortedX(1:N);
            maxValueIndices = sortingIndices(1:minimumEvents);

            [~, peaksTmp] = computeSignalPeaks(tmpS,'makePlots', 0,'makeSummaryPlots',0,'waitbarOn',0,'numStdsForThresh',2.5,'outputInfo',0);
            if length(peaksTmp{1})>=minimumEvents
                maxValueIndices = [peaksTmp{1}(1:minimumEvents)];
            else
                maxValueIndices = [peaksTmp{1}];
            end
            tmpIdx = randperm(length(tmpS(:)),2);
            maxValueIndices = [maxValueIndices(:); tmpIdx(:)];
            signalPeakIdx{idxHere} = unique([maxValueIndices(:); signalPeakIdx{idxHere}(:)]);
            % signalPeakIdx{idxHere}
            signalPeaks(idxHere,maxValueIndices) = 1;
            % figure(idxHere)
            % plot(inputSignals(idxHere,:));
            % hold on
            % plot(signalPeakIdx{idxHere},inputSignals(idxHere,signalPeakIdx{idxHere}),'r+','MarkerSize',10);
        end
    end

    % get the peak statistics
    [peakOutputStat] = computePeakStatistics(inputSignals,'waitbarOn',1,'testpeaks',options.signalPeaks,'testpeaksArray',options.signalPeaksArray,'spikeROI',options.peakROI);

    %
    if ~isempty(options.inputMovie)
        [~, outputMeanImageCorrs] = createPeakTriggeredImages(options.inputMovie, inputImages, inputSignals,'cropSize',options.cropSize);
        outputMeanImageCorrs(isnan(outputMeanImageCorrs)) = 0;
        peakOutputStat.outputMeanImageCorrs = outputMeanImageCorrs(:);
        % get ROI traces
        if options.showROITrace==1
            [ROItraces] = applyImagesToMovie(inputImagesThres,options.inputMovie,'alreadyThreshold',1,'waitbarOn',1);
        else
            ROItraces = [];
        end
    else
        peakOutputStat.outputMeanImageCorrs = [];
    end
    % =======

    % display histogram of movie/trace
    histogramSwitch = 0;
    if ~isempty(options.inputMovie)&histogramSwitch==1
        figure(45684)
        subplot(3,1,1)
        hist(inputSignals(:),100);xlabel('input values');ylabel('counts')
        title(['input | min: ' num2str(nanmin(inputSignals(:))) ' | max: ' num2str(nanmax(inputSignals(:)))])
        subplot(3,1,2)
        hist(options.inputMovie(:),100);xlabel('movie values');ylabel('counts')
        title(['movie | min: ' num2str(nanmin(options.inputMovie(:))) ' | max: ' num2str(nanmax(options.inputMovie(:)))])
        subplot(3,1,3)
        hist(ROItraces(:),100);xlabel('ROI values');ylabel('counts')
        title(['ROI | min: ' num2str(nanmin(ROItraces(:))) ' | max: ' num2str(nanmax(ROItraces(:)))])
    end

    % remove small ICs unless a pre-list is loaded in
    if isempty(valid)
        [~, ~, valid, inputImageSizes] = filterImages(inputImages, inputSignals,'thresholdImages',1);
        %
        validPre = valid;
        % pre-select as valid if SNR is above a certain threshold
        validSNR = signalSnr>options.SnrThreshold;
        validPre = valid | validSNR;
        validSlope = peakOutputStat.slopeRatio>options.slopeRatioThreshold;
        validPre = validPre & validSlope;
        % Since 0=invalid, 1=valid, -1=unknown, set all '1' to unknown
        valid(find(valid==1)) = -1;
    else
        % [~, ~, ~, inputImageSizes] = filterImages(inputImagesThres, inputSignals,'thresholdImages',0);
        inputImageSizes = sum(sum(inputImagesThres,1),2);
        inputImageSizes = inputImageSizes(:);
        validPre = valid;
    end

    % =======
    % plot information about the traces
    % plotSignalStatistics(inputSignals,inputImageSizes,inputStr,'r','hold off',signalSnr,peakOutputStat.slopeRatio)
    plotSignalStatisticsWrapper(inputSignals,inputImages,validPre,inputImageSizes,inputStr,signalSnr,peakOutputStat);

    % =======
    % loop over choices
    nSignals = size(inputImages,3);
    display(['# signals: ' num2str(nSignals)]);
    % valid = ones(1,size(inputImages,1))*-1;

    choices = chooseSignals(options,1:nSignals, inputImages,inputSignals,objMap, valid, inputStr,tmpDir,sessionID,signalPeakIdx,signalSnr,inputImagesThres,inputImageSizes,peakOutputStat,ROItraces,imgStats);
    % assume all skips were good ICs that user forgot to enter
    validChoices = choices;
    validChoices(find(validChoices==-1))=1;
    validChoices = logical(validChoices);

    % =======
    % plotSignalStatisticsWrapper(inputSignals,inputImages,validChoices,inputImageSizes,inputStr,signalSnr,peakOutputStat);

    % if manually scoring automatic, combine manual classification with automatic
    if ~isempty(validBackup)&~isempty(find(validBackup==2))
        valid = validBackup;
        % add the manual scores for the questionable signals into the valid input vector
        % validChoices
        % questionableSignalIdx
        valid(questionableSignalIdx) = validChoices;
        validChoices = logical(valid);
        choices = validChoices;
        % restore original input data
        inputImages = inputImagesBackup;
        inputSignals = inputSignalsBackup;
    end

    % =======
    % filter input for valid signals
    inputImages = inputImages(:,:,validChoices);
    inputSignals = inputSignals(validChoices,:);
end
function plotSignalStatisticsWrapper(inputSignals,inputImages,validChoices,inputImageSizes,inputStr,signalSnr,peakOutputStat)
    % plot good and bad signals with different colors

    % determine number of IC filters to investigate
    pointColor = ['r','g'];
    for pointNum = 1:2
        if pointNum==1
            valid = logical(~validChoices);
        else
            valid = logical(validChoices);
        end
        % plot information about the traces
        plotSignalStatistics(inputSignals(valid,:),inputImageSizes(valid),inputStr,pointColor(pointNum),'hold on',signalSnr(valid),peakOutputStat.slopeRatio(valid))
    end
end
function plotSignalStatistics(inputSignals,inputImageSizes,inputStr,pointColor, holdState,signalSnr,slopeRatio)
    % plot statistics for input signal

    % get best fit line SNR v slopeRatio
    p = polyfit(signalSnr,slopeRatio,1);   % p returns 2 coefficients fitting r = a_1 * x + a_2
    r = p(1) .* signalSnr + p(2); % compute a new vector r that has matching datapoints in x
    if ~isempty(slopeRatio)&~isempty(signalSnr)
        % start plotting!
        figNo = 1776;%AMERICA
        [figHandle figNo] = openFigure(figNo, '');
        hold off;
        plot(normalizeVector(slopeRatio),'Color',[4 4 4]/5);hold on;
        plot(normalizeVector(signalSnr),'r');
        title(['SNR in trace signal for ' inputStr])
        hleg1 = legend('S-ratio','SNR');
        xlabel('ic rank');ylabel('SNR');box off;hold off;

        [figHandle figNo] = openFigure(figNo, '');
        hold off;
        plot(slopeRatio,'Color',[4 4 4]/5);hold on;
        plot(signalSnr,'r');
        title(['SNR in trace signal for ' inputStr])
        hleg1 = legend('S-ratio','SNR');
        xlabel('ic rank');ylabel('SNR');box off;hold off;

        [figHandle figNo] = openFigure(figNo, '');
        scatter(signalSnr,slopeRatio,[pointColor '.']);hold on;
        plot(signalSnr, r, 'k-');
        title(['SNR v S-ratio for ' inputStr])
        xlabel('SNR');ylabel('S-ratio');box off;
        eval(holdState);

        [figHandle figNo] = openFigure(figNo, '');
        scatter3(signalSnr,slopeRatio,inputImageSizes,[pointColor '.'])
        title(['SNR, S-ratio, filter size for ' inputStr])
        xlabel('SNR');ylabel('S-ratio');zlabel('ic size');
        legend({'bad','good'});rotate3d on;
        eval(holdState);
    end
end
function [valid] = chooseSignals(options,signalList, inputImages,inputSignals,objMap, valid, inputStr,tmpDir,sessionID,signalPeakIdx,signalSnr,inputImagesThres,inputImageSizes,peakOutputStat,ROItraces,imgStats)
    % manually decide which signals are good or bad, pre-computed values input to speed up movement through signals

    warning('off','all')
    warning('query','all')

    if ~exist(tmpDir,'file')
        mkdir(tmpDir);
    end

    % mainFig = openFigure(1,'full');
    mainFig = figure(1);
    % prevent matlab from giving command window focus
    set(mainFig,'KeyPressFcn', '1;');

    % location of each subplot
    if ~isempty(options.inputMovie)
        objMapPlotLoc = 4;
        tracePlotLoc = [5:6];
    else
        objMapPlotLoc = 1;
        tracePlotLoc = [4:6];
    end
    % inputMoviePlotLoc = 1:2;
    inputMoviePlotLoc = 2;
    inputMoviePlotLoc2 = 1;
    filterPlotLoc = inputMoviePlotLoc2;
    avgSpikeTracePlot = 3;
    subplotX = 3;
    subplotY = 2;

    if ~isempty(options.inputMovie)
        objMapPlotLoc = [1 2 7 8];
        tracePlotLoc = [10 11 12];
    else
        objMapPlotLoc = 1;
        tracePlotLoc = [4:6];
    end
    % inputMoviePlotLoc = 1:2;
    inputMoviePlotLoc = [3 4];
    inputMoviePlotLoc2 = [5 6];
    filterPlotLoc = inputMoviePlotLoc2;
    avgSpikeTracePlot = 9;
    subplotX = 6;
    subplotY = 2;

    % subplot(subplotY,subplotX,objMapPlotLoc);
    % subplot(subplotY,subplotX,inputMoviePlotLoc);
    % subplot(subplotY,subplotX,avgSpikeTracePlot);
    % title('tmp')
    % subplot(subplotY,subplotX,tracePlotLoc);
    % title('tmp')

    % instructions
    instructionStr =  ['up/down:good/bad | left/right: forward/back | m:(montage peak images) | c:(compare signal to movie) |',10,' f:finished | g:(goto signal) | s:(set remaining signals to bad) | signals are assumed good',10,10,10];
    instructionStr =  ['controls',10,10,'up/down:good/bad',10,'left/right: forward/back',10,'m: peak images',10,'c: movie signal',10,'f:finished',10,'g: goto signal',10,'s: remaining bad',10,'signals assumed good'];
    % sepStrNum = ' | ';
    sepStrNum = 10;
    instructionStr =  ['up/down: good/bad',sepStrNum,'left/right: forward/back',sepStrNum,'m: peak images',sepStrNum,'c: movie signal',sepStrNum,'x: movie signal array',10,'g: goto signal',sepStrNum,'s: set remaining bad',sepStrNum,'f: finished/save',10,...
    'Eccentricity>0.4',10,...
    'imageSizes>10,<100',10,...
    'Perimeter<50,>5',10,...
    'Solidity>0.8',10,...
    'EquivDiameter>3,<30',10,...
    'signalSnr>1.45',10,...
    'slopeRatio>0.02',10];
    suptitleHandle = suptitle(instructionStr);

    set(suptitleHandle,'FontSize',10,'FontWeight','normal');
    set(suptitleHandle, 'horizontalAlignment', 'left');
    set(suptitleHandle, 'units', 'normalized');
    h1 = get(suptitleHandle, 'position');
    set(suptitleHandle, 'position', [0.007 -0.25 h1(3)]);

    % plot the cell map to provide context
    subplot(subplotY,subplotX,objMapPlotLoc);
    imagesc(objMap); axis off; colormap gray;
    title(['objMap' inputStr]);hold on;

    % make color image overlays
    zeroMap = zeros(size(objMap));
    oneMap = ones(size(objMap))*0.5;
    green = cat(3, zeroMap, oneMap, zeroMap);
    blue = cat(3, zeroMap, zeroMap, oneMap);
    red = cat(3, oneMap, zeroMap, zeroMap);
    warning off
    imageOverlay = imshow(blue);
    goodFilterOverlay = imshow(green);
    badFilterOverlay = imshow(red);
    warning on
    hold off

    % get values for plotting
    % options.peakROI = [-40:40]
    peakROI = options.peakROI;
    minValTraces = nanmin(inputSignals(:));
    if minValTraces<options.minValConstant
        minValTraces = options.minValConstant;
    end
    minValTraces
    maxValTraces = nanmax(inputSignals(:));
    if maxValTraces>0.4|maxValTraces<0.3
        % maxValTraces = 0.35;
    end

    % filter based on the list
    inputImages = inputImages(:,:,signalList);
    inputSignals = inputSignals(signalList,:);

    % loop over chosen filters
    nImages = size(inputImages,3);

    % initialize loop variables
    saveData=0;
    i = 1;
    reply = 0;
    loopCount = 1;
    warning off

    if ~isempty(options.inputMovie)
        display('calculating movie min/max...')
        % minValMovie = nanmin(options.inputMovie(:));
        % maxValMovie = nanmax(options.inputMovie(:));
        maxValMovie = options.movieMax;
        minValMovie = options.movieMin;
    end

    [xCoords yCoords] = findCentroid(inputImages,'thresholdValue',0.8,'imageThreshold',options.threshold);

    % pre-calculate
    % if ~isempty(options.inputMovie)
    %     croppedPeakImages = {};
    %     reverseStr = '';
    %     for j=1:nImages
    %         figure(79879);
    %         thisTrace = inputSignals(j,:);
    %         [croppedPeakImages{j}] = viewMontage(options.inputMovie,inputImages(j,:,:),thisTrace);
    %         reverseStr = cmdWaitbar(j,nImages,reverseStr,'inputStr','getting montages','waitbarOn',1,'displayEvery',5);
    %     end
    % end

    % display('pre-loading preview images')
    objCutMovieCollection = {};
    reverseStr = '';
    nSignalsHere = size(inputImages,3);
    for signalNo = 1:nSignalsHere
        thisTrace = inputSignals(signalNo,:);
        testpeaks = signalPeakIdx{signalNo};
        try
            [objCutMovieCollection{signalNo}] = createObjCutMovieSignalSorter(options,testpeaks,thisTrace,inputImages,signalNo,options.cropSizeLength,maxValMovie);
        catch err
            objCutMovieCollection{signalNo} = {};
            display(repmat('@',1,7))
            disp(getReport(err,'extended','hyperlinks','on'));
            display(repmat('@',1,7))
        end
        if mod(signalNo,20)==0;reverseStr = cmdWaitbar(signalNo,nSignalsHere,reverseStr,'inputStr','pre-loading preview images','waitbarOn',1,'displayEvery',1);end
    end

    % only exit if user clicks options that calls for saving the data
    while saveData==0
        figure(mainFig);
        % change figure color based on nature of current choice
        % valid
        if valid(i)==1
            set(mainFig,'Color',options.backgroundGood);
        elseif valid(i)==0
            set(mainFig,'Color',options.backgroundBad);
        else
            set(mainFig,'Color',options.backgroundNeutral);
        end

        % get loop specific values
        directionOfNextChoice=0;
        thisImage = squeeze(inputImages(:,:,i));
        thisTrace = inputSignals(i,:);
        cellIDStr = ['#' num2str(i) '/' num2str(nImages)];

        if ~isempty(options.inputMovie)
            subplot(subplotY,subplotX,inputMoviePlotLoc2)
                % tic
                testpeaks = signalPeakIdx{i};
                if(~isempty(testpeaks))
                    try
                        [croppedPeakImages] = viewMontage(options.inputMovie,inputImages(:,:,i),options,thisTrace,[signalPeakIdx{i}],minValMovie,maxValMovie,options.cropSizeLength);
                    catch err
                        imagesc(thisImage);
                        display(repmat('@',1,7))
                        disp(getReport(err,'extended','hyperlinks','on'));
                        display(repmat('@',1,7))
                    end
                    try
                        sigDig = 100;
                        title([...
                            'imageCorr = ' num2str(round(peakOutputStat.outputMeanImageCorrs(i)*sigDig)/sigDig)...
                            ' | SNR = ' num2str(round(signalSnr(i)*sigDig)/sigDig)...
                            ' | S-ratio = ' num2str(round(peakOutputStat.slopeRatio(i)*sigDig)/sigDig)]);

                    catch err
                        % title([' (' num2str(sum(valid==1)) ' good)']);
                        display(repmat('@',1,7))
                        disp(getReport(err,'extended','hyperlinks','on'));
                        display(repmat('@',1,7))
                    end
                else
                    imagesc(thisImage);
                    % title([' (' num2str(sum(valid==1)) ' good)']);
                    % colormap(customColormap([]));
                end
                % toc
                % imagesc(croppedPeakImages{i});
                % colormap(customColormap([]));
                % axis off;
        else
            % show the current image
            subplot(subplotY,subplotX,filterPlotLoc)
                imagesc(thisImage);
                colormap gray
                axis off; % ij square
                title(['signal ' cellIDStr 10 '(' num2str(sum(valid==1)) ' good)']);
        end

        % use thresholded image as AlphaData to overlay on cell map, reduce number of times this is accessed to speed-up analysis
        if loopCount==1|~exist('Comb','var')
            Comb(:,:,1) = zeros([size(squeeze(inputImagesThres(:,:,i)))]);
            Comb(:,:,2) = zeros([size(squeeze(inputImagesThres(:,:,i)))]);
            Comb(:,:,3) = zeros([size(squeeze(inputImagesThres(:,:,i)))]);
            display(num2str([min(Comb(:)) max(Comb(:))]))
            CombTmp = Comb;
            % colormap gray
            goodImages = createObjMap(inputImagesThres(:,:,valid==1));
            if(isempty(goodImages)) goodImages = zeros(size(objMap)); end
            badImages = createObjMap(inputImagesThres(:,:,valid==0));
            if(isempty(badImages)) badImages = zeros(size(objMap)); end
            if sum(valid==3)>0
                badImages = zeros(size(objMap));
                neutralImages = createObjMap(inputImagesThres(:,:,valid==3));
            else
                neutralImages = zeros(size(objMap));
            end
            % if(isempty(badImages)) badImages = zeros(size(objMap)); end
        end
        % if mod(loopCount,20)==0|loopCount==1
        %     subplot(subplotY,subplotX,objMapPlotLoc)
        %         % set(goodFilterOverlay, 'AlphaData', goodImages);
        %         % set(badFilterOverlay, 'AlphaData', badImages);
            % set(imageOverlay, 'AlphaData', squeeze(inputImagesThres(i,:,:)));
            % goodImages2 = normalizeVector(goodImages,'normRange','zeroToOne');
            % badImages2 = normalizeVector(badImages,'normRange','zeroToOne');
        % end
        subplot(subplotY,subplotX,objMapPlotLoc)
            currentImage = squeeze(inputImages(:,:,i));
            currentImageThres = squeeze(inputImagesThres(:,:,i));
            Comb(:,:,2) = ones([size(currentImageThres)]);

            CombTmp(:,:,1) = badImages; %red
            CombTmp(:,:,2) = goodImages; %green
            CombTmp(:,:,3) = neutralImages; %blue
            % if(~isempty(neutralImages))
            %     CombTmp(:,:,1) = neutralImages; %red
            %     CombTmp(:,:,2) = neutralImages; %green
            %     CombTmp(:,:,3) = neutralImages; %blue
            % end
            switch valid(i)
                case 0
                    CombTmp(:,:,1) = CombTmp(:,:,1)-currentImageThres;
                case 1
                    CombTmp(:,:,2) = CombTmp(:,:,2)-currentImageThres;
                case 3
                    CombTmp(:,:,3) = CombTmp(:,:,3)-currentImageThres;
                    % CombTmp(:,:,2) = CombTmp(:,:,2)-currentImageThres;
                    % CombTmp(:,:,3) = CombTmp(:,:,2)-currentImageThres;
                otherwise
                    % body
            end
            CombTmp(:,:,1) = CombTmp(:,:,1)+currentImageThres;
            CombTmp(:,:,2) = CombTmp(:,:,2)+currentImageThres;
            CombTmp(:,:,3) = CombTmp(:,:,3)+currentImageThres;
            % re-color cellmap
            tmpImage = squeeze(CombTmp(:,:,3));
            neutralIdx = tmpImage>0;
            % lighten the good/bad colors
            for dimNo = 1:2
                CombTmpMain = squeeze(CombTmp(:,:,dimNo));
                alterIdx{dimNo} = CombTmpMain>0;
            end
            for dimNo = 1:2
                CombTmpMain = squeeze(CombTmp(:,:,dimNo));
                for dimNo2 = 1:size(CombTmp,3)
                    if dimNo==dimNo2;continue;end
                    CombTmp1 = squeeze(CombTmp(:,:,dimNo2));
                    CombTmp1(alterIdx{dimNo}) = CombTmpMain(alterIdx{dimNo})/4;
                    CombTmp(:,:,dimNo2) = CombTmp1;
                end
            end
            currentIdx = currentImageThres>0;
            % set not picked to gray
            for dimNo = 1:size(CombTmp,3)
                CombTmp1 = squeeze(CombTmp(:,:,dimNo));
                CombTmp1(neutralIdx) = tmpImage(neutralIdx);
                % set current to blue
                if dimNo~=3
                    CombTmp1(currentIdx) = 0;
                end
                CombTmp(:,:,dimNo) = CombTmp1;

            end

            % xCoords yCoords
            CombTmp2 = CombTmp;
            CombTmp2(yCoords(i),:,:) = 1;
            CombTmp2(:,xCoords(i),:) = 1;
            CombTmp23 = CombTmp2(:,:,1); CombTmp23(currentIdx) = NaN;
                % CombTmp23(currentIdx) = currentImageThres(currentIdx);
                CombTmp2(:,:,1) = CombTmp23;
            CombTmp23 = CombTmp2(:,:,2); CombTmp23(currentIdx) = NaN;
                % CombTmp23(currentIdx) = currentImageThres(currentIdx);
                CombTmp2(:,:,2) = CombTmp23;
            CombTmp23 = CombTmp2(:,:,3); CombTmp23(currentIdx) = NaN;
                CombTmp23(currentIdx) = currentImageThres(currentIdx);
                CombTmp2(:,:,3) = CombTmp23;

            imagesc(CombTmp2)
            % title(['signal map' 10 'green/red/gray/blue' 10 'good/bad/undecided/current'])
            title(['green/red/gray/blue' 10 'good/bad/undecided/current'])

        % if signal has peaks, plot the average signal and other info
        testpeaks = signalPeakIdx{i};
        if(~isempty(testpeaks))
            % tic
            % plot all signals and the average
            subplot(subplotY,subplotX,avgSpikeTracePlot);
                [slopeRatio] = plotPeakSignal(thisTrace,testpeaks,cellIDStr,instructionStr,minValTraces,maxValTraces,peakROI,peakOutputStat.avgSpikeTrace(i,:),peakOutputStat.slopeRatio(i),peakOutputStat.spikeCenterTrace{i},valid);
            % add in the ratio of the rise/decay slopes. Should be >>1 for calcium
            % subplot(subplotY,subplotX,filterPlotLoc)
                % title(['signal ' cellIDStr]);
            % plot the trace

            subplot(subplotY,subplotX,tracePlotLoc)
                sigDig = 100;
                thisStr = [...
                    'SNR = ' num2str(round(signalSnr(i)*sigDig)/sigDig)...
                    ' | S-ratio = ' num2str(round(peakOutputStat.slopeRatio(i)*sigDig)/sigDig)...
                    ' | # peaks = ' num2str(length(testpeaks))...
                    10 ...
                    'size (px) = ' num2str(round(inputImageSizes(i)*sigDig)/sigDig)...
                    ' | imageCorr = ' num2str(round(peakOutputStat.outputMeanImageCorrs(i)*sigDig)/sigDig)...
                    ' | Eccentricity = ' num2str(round(imgStats.Eccentricity(i)*sigDig)/sigDig)...
                    10 ...
                    'Perimeter = ' num2str(round(imgStats.Perimeter(i)*sigDig)/sigDig)...
                    ' | Solidity = ' num2str(round(imgStats.Solidity(i)*sigDig)/sigDig)...
                    ' | EquivD = ' num2str(round(imgStats.EquivDiameter(i)*sigDig)/sigDig)];
                plotSignal(thisTrace,testpeaks,'',thisStr,minValTraces,maxValTraces,options);
                if ~isempty(options.inputMovie)&options.showROITrace==1
                    hold on
                    % [tmpTrace] = applyImagesToMovie(inputImagesThres(:,:,i),options.inputMovie,'alreadyThreshold',1,'waitbarOn',0);
                    tmpTrace = ROItraces(i,:);
                    tmpTrace = squeeze(tmpTrace);
                    % nanmean(tmpTrace)
                    % nanmin(tmpTrace)
                    % nanmax(tmpTrace)
                    if abs(nanmax(tmpTrace))<abs(nanmin(tmpTrace))
                        tmpTrace = -tmpTrace;
                    end
                    % tmpTrace = tmpTrace+nanmin(tmpTrace);
                    % tmpTrace = (tmpTrace-nanmean(tmpTrace))/nanmean(tmpTrace);
                    traceRatio = nanmax(thisTrace)/nanmax(tmpTrace);
                    tmpTrace = tmpTrace*traceRatio+0.1;
                    % tmpTrace = nanmax(thisTrace)*normalizeVector(tmpTrace,'normRange','zeroToOne')+0.1;
                    plot(tmpTrace,'k');
                    legend('original','ROI');legend boxoff;
                    axis([0 length(thisTrace) minValTraces maxValTraces+0.1]);
                    hold off
                end
            % toc
        else
            subplot(subplotY,subplotX,avgSpikeTracePlot);
                plot(peakROI,thisTrace(1:length(peakROI)));
                xlabel('frames');
                ylabel('\DeltaF/F');
                ylim([minValTraces maxValTraces]);
                title(['signal peaks ' cellIDStr])
            subplot(subplotY,subplotX,tracePlotLoc)
                plot(thisTrace, 'r');
                xlabel('frames');
                % ylabel('df/f');
                axis([0 length(thisTrace) minValTraces maxValTraces]);
                thisStr = ['SNR = ' num2str(signalSnr(i)) ' | S-ratio = ' num2str(NaN) ' | # peaks = ' num2str(length(testpeaks)) ' | size (px) = ' num2str(inputImageSizes(i))];
                title([thisStr])
        end

        % get user input
        figure(mainFig);
        set(findobj(gcf,'type','axes'),'hittest','off')
        validPrevious = valid(i);
        warning('off','all')

        % show object cut movie in loop
        % if ~isempty(options.inputMovie)
        %     subplot(subplotY,subplotX,inputMoviePlotLoc)
        %         % [objCutMovie] = createObjCutMovieSignalSorter(options,testpeaks,thisTrace,inputImages,i,)
        % end
        set(gcf,'pointer','custom','PointerShapeCData',NaN([16 16]));
        if ~isempty(options.inputMovie)
            try
                % [objCutMovie] = createObjCutMovieSignalSorter(options,testpeaks,thisTrace,inputImages,i);
                objCutMovie = objCutMovieCollection{i};
                if options.movieMin<-0.025
                    objCutMovie(1,1,:) = -0.025;
                else
                    objCutMovie(1,1,:) = options.movieMin;
                end
                objCutMovie(1,2,:) = options.movieMax*0.4;

                subplot(subplotY,subplotX,inputMoviePlotLoc)
                    % set title and turn off ability to change
                    if ~isempty(objCutMovie)
                        imagesc(objCutMovie(:,:,1));drawnow
                    else

                    end
                    title(strrep(options.inputStr,'_','\_'), 'HandleVisibility' , 'off' )
                    % set(gca , 'NextPlot' , 'replacechildren');
                    axis off
                % title(['signal ' cellIDStr ' (' num2str(sum(valid==1)) ' good)']);
                frameNoMax = size(objCutMovie,3);
                frameNo = round(frameNoMax*0.40);
                % keyIn = [];
                set(gcf,'currentch','3');
                keyIn = get(gcf,'CurrentCharacter');
                % strcmp(keyIn,'3')
                options.fps = 25;

                % writerObj = VideoWriter(['cell' num2str(i) '.avi']);
                % open(writerObj);

                colormap(options.colormap);

                while strcmp(keyIn,'3')
                    keyIn = get(gcf,'CurrentCharacter');
                    if ~isempty(objCutMovie)
                        imagesc(objCutMovie(:,:,frameNo));
                        drawnow
                    else

                    end
                    % axis off;
                    % if isempty(double(keyIn))
                        % keyIn = '/';
                    % end
                    if frameNo==frameNoMax
                        frameNo = 1;
                    end
                    pause(1/options.fps);
                    frameNo = frameNo + 1;

                    % writeVideo(writerObj,getframe(mainFig));

                end

                reply = double(keyIn);
                set(gcf,'currentch','3');

                % close(writerObj);

            catch err
                display(repmat('@',1,7))
                disp(getReport(err,'extended','hyperlinks','on'));
                display(repmat('@',1,7))
                [x,y,reply]=ginput(1);
            end
        else
            [x,y,reply]=ginput(1);
        end
        % objCutMovie(:,:,1) = nanmax(objCutMovie(:));
        % playMovie(objCutMovie,'fps',15);
        warning('on','all')
        % [x,y,reply]=finput(1);
        % reply = uint8(input);
        % reply = uint8(input('enter: ','s'))
        % waitforbuttonpress
        % reply = double(get(gcf,'CurrentCharacter'));

        % 'E' make a montage of peak frames
        if isequal(reply, 109)&~isempty(options.inputMovie)
            try
                [croppedPeakImages] = viewMontage(options.inputMovie,inputImages(:,:,i),options,thisTrace,[signalPeakIdx{i}],minValMovie,maxValMovie,options.cropSizeLength);
                set(findobj(gcf,'type','axes'),'hittest','off')
                ginput(1);
            catch err
                display(repmat('@',1,7))
                disp(getReport(err,'extended','hyperlinks','on'));
                display(repmat('@',1,7))
            end
            % close(2);figure(mainFig);
        % 'C' compare signal to movie
        elseif isequal(reply, 99)&~isempty(options.inputMovie)
            signalPeakArray = [signalPeakIdx{i}];
            peakSignalAmplitude = thisTrace(signalPeakArray(:));
            [peakSignalAmplitude peakIdx] = sort(peakSignalAmplitude,'descend');
            signalPeakArray = {signalPeakArray(peakIdx)};
            try
                compareSignalToMovie(options.inputMovie, inputImages(:,:,i), thisTrace,'waitbarOn',0,'timeSeq',-10:10,'signalPeakArray',signalPeakArray,'cropSize',options.cropSizeLength);
            catch err
                display(repmat('@',1,7))
                disp(getReport(err,'extended','hyperlinks','on'));
                display(repmat('@',1,7))
                display('Error displaying movie signal')
            end
        % 'E' montage cut movie
        elseif (isequal(reply, 120)|isequal(reply, 2)|isequal(reply, 48))&~isempty(options.inputMovie)
            try
                [objCutMovie] = createObjCutMovieSignalSorter(options,testpeaks,thisTrace,inputImages,i,maxValMovie);
                % objCutMovie(:,:,1) = nanmax(objCutMovie(:));
                playMovie(objCutMovie,'fps',15);
            catch err
                display(repmat('@',1,7))
                disp(getReport(err,'extended','hyperlinks','on'));
                display(repmat('@',1,7))
            end
        else
            [valid directionOfNextChoice saveData i] = respondToUserInput(reply,i,valid,directionOfNextChoice,saveData,nImages);
        end

        % update images
        if validPrevious==valid(i)
        elseif valid(i)==1
            goodImages = goodImages+currentImageThres;goodImages(goodImages<0) = 0;
            badImages = badImages-currentImageThres;badImages(badImages<0) = 0;
            if validPrevious==3
                neutralImages = neutralImages-currentImageThres;
            end
        elseif valid(i)==0
            goodImages = goodImages-currentImageThres;goodImages(goodImages<0) = 0;
            badImages = badImages+currentImageThres;badImages(badImages<0) = 0;
            if validPrevious==3
                neutralImages = neutralImages-currentImageThres;
            end
        end
        % loop if user gets to either end
        i=i+directionOfNextChoice;
        if i<=0 i = nImages; end
        if i>nImages i = 1; end
        pause(0.001);
        figure(mainFig);

        % already checked that tmp folder exists, then save
        try
            if exist(tmpDir,'file')
                save([tmpDir filesep 'tmpDecisions_' sessionID '.mat'],'valid');
            end
        catch

        end

        loopCount = loopCount+1;
    end
    % warning on
    warning('on','all')
    warning('query','all')
end

function [objCutMovie] = createObjCutMovieSignalSorter(options,testpeaks,thisTrace,inputImages,i,cropSizeLength,maxValMovie)
    preOffset = 10;
    postOffset = 10;
    timeVector = [-preOffset:postOffset]';
    nPoints = size(options.inputMovie,3);
    maxSignalsToShow = options.maxSignalsToShow-1;
    % testpeaks= unique(testpeaks);
    % bias toward high amplitude signals
    peakSignalAmplitude = thisTrace(testpeaks(:));
    % peakSignalAmplitude
    [peakSignalAmplitude peakIdx] = sort(peakSignalAmplitude,'descend');
    % peakSignalAmplitude
    testpeaks = testpeaks(peakIdx);

    if length(testpeaks)==1
        testpeaks(end+1:end+2) = testpeaks;
    end

    testpeaks((testpeaks-preOffset)<1) = [];
    testpeaks((testpeaks+postOffset)>length(thisTrace)) = [];

    if length(testpeaks)>maxSignalsToShow
        % choose a random subset
        framesToAlign = testpeaks(1:maxSignalsToShow);
    else
        framesToAlign = testpeaks;
    end

    % thisTrace(framesToAlign(:))

    %remove points outside valid range
    % framesToAlign(find((framesToAlign-preOffset)<1)) = [];
    % framesToAlign(find((framesToAlign>(nPoints-postOffset)))) = [];
    peakIdxs = bsxfun(@plus,timeVector,framesToAlign(:)');
    nAlignPts = length(framesToAlign(:));
    % remove frame alignment outside range
    peakIdxs(find((peakIdxs<1))) = [];
    peakIdxs(find((peakIdxs>nPoints))) = [];
    % get movie cut around cell
    objCutMovie = getObjCutMovie(options.inputMovie(:,:,peakIdxs),inputImages(:,:,i),'createMontage',0,'extendedCrosshairs',2,'crossHairVal',maxValMovie*options.crossHairPercent,'outlines',1,'waitbarOn',0,'cropSize',cropSizeLength);
    objCutMovie = vertcat(objCutMovie{:});
    % objCutMovieTmp = objCutMovie{1};
    dimSize = [size(objCutMovie,1) size(objCutMovie,2) length(timeVector)];
    diffDiffMatrix = NaN(dimSize);
    diffDiffMatrix2 = diffDiffMatrix;
    % diffDiffMatrix(:,:,1) = nanmax(inputMovie(:));
    % diffDiffMatrix(:,:,1) = options.movieMean;
    % diffDiffMatrix(:,:,1) = options.movieMax*0.75;
    % diffDiffMatrix(:,:,round(size(diffDiffMatrix,3)/2)) = options.movieMax;
    nStimMovies = 2;
    for iii = 1:nStimMovies;objCutMovie = cat(3, diffDiffMatrix, objCutMovie);end
    dimDiff = maxSignalsToShow + nStimMovies - round(size(objCutMovie,3)/length(timeVector));
    % display('===')
    % dimDiff
    for diffNo = 1:dimDiff
        objCutMovie = cat(3, objCutMovie, diffDiffMatrix2);
    end

    % croppedPeakImages = compareSignalToMovie(options.inputMovie, inputImages(:,:,i), thisTrace,'crosshairs',0,'getOnlyPeakImages',1,'waitbarOn',0,'extendedCrosshairs',0,'crossHairVal',maxValMovie*options.crossHairPercent,'outlines',0,'signalPeakArray',{testpeaks},'cropSize',cropSizeLength);

    % [thresholdedImages boundaryIndices] = thresholdImages(croppedPeakImages(:,:,1),'binary',1,'getBoundaryIndex',1,'threshold',options.threshold/2);
    % for imageNo = 1:size(objCutMovie,3)
    %     tmpImg = objCutMovie(:,:,imageNo);
    %     % tmpImg(boundaryIndices{1}) = tmpImg(boundaryIndices{1})+maxValMovie*options.crossHairPercent;
    %     tmpImg(boundaryIndices{1}) = NaN;
    %     objCutMovie(:,:,imageNo) = tmpImg;
    % end

    % playMovie(objCutMovie);
    nAlignPts = nAlignPts+nStimMovies+dimDiff;
    [objCutMovie] = createStimCutMovieMontage(objCutMovie,nAlignPts,timeVector,'squareMontage',1,'addStimMovie',0);
end
function [croppedPeakImages2] = viewMontage(inputMovie,inputImage,options,thisTrace,signalPeakArray,minValMovie,maxValMovie,cropSizeLength)

    if isempty(signalPeakArray)
        imagesc(inputImage);
        colormap(options.colormap);
        axis off;
        croppedPeakImages2 = inputImage;
        return
    end
    % signalPeakArray
    maxSignalsToShow = options.maxSignalsToShow-1;
    peakSignalAmplitude = thisTrace(signalPeakArray(:));
    % peakSignalAmplitude
    [peakSignalAmplitude peakIdx] = sort(peakSignalAmplitude,'descend');
    % peakSignalAmplitude
    signalPeakArray = signalPeakArray(peakIdx);
    if length(signalPeakArray)>maxSignalsToShow
        % choose a random subset
        signalPeakArray = signalPeakArray(1:maxSignalsToShow);
    end

    signalPeakArray((signalPeakArray-31)<1) = [];
    signalPeakArray((signalPeakArray+31)>length(thisTrace)) = [];

    signalPeakArray = {signalPeakArray};
    % signalPeakArray
    croppedPeakImages = compareSignalToMovie(inputMovie, inputImage, thisTrace,'getOnlyPeakImages',1,'waitbarOn',0,'extendedCrosshairs',2,'crossHairVal',maxValMovie*options.crossHairPercent,'outlines',1,'signalPeakArray',signalPeakArray,'cropSize',cropSizeLength);
    % display cropped images
    % figure(2);
    % for i=1:size(croppedPeakImages,3)
    %     kurtosisM(1,i) = kurtosis(sum(squeeze(croppedPeakImages(i,:,:)),1));
    %     kurtosisM(2,i) = kurtosis(sum(squeeze(croppedPeakImages(i,:,:)),2));
    %     % [kurtosisX kurtosisY]
    % end
    % kurtosisM'
    meanTransientImageTmp = nanmean(croppedPeakImages(2:end-1,2:end-1,2:end),3);
    meanTransientImageTmp = padarray(meanTransientImageTmp,[1 1],maxValMovie*0.3);
    meanTransientImage = zeros([size(meanTransientImageTmp,1) size(meanTransientImageTmp,2) 1]);
    meanTransientImage(:,:,1) = meanTransientImageTmp;
    % size(meanTransientImage)
    % size(croppedPeakImages)
    croppedPeakImages = cat(3,meanTransientImage,croppedPeakImages);

    croppedPeakImages222 = compareSignalToMovie(inputMovie, inputImage, thisTrace,'getOnlyPeakImages',1,'waitbarOn',0,'extendedCrosshairs',2,'crossHairVal',maxValMovie*options.crossHairPercent,'outlines',0,'signalPeakArray',signalPeakArray,'cropSize',cropSizeLength,'crosshairs',0);
    [thresholdedImages boundaryIndices] = thresholdImages(croppedPeakImages222(:,:,1),'binary',1,'getBoundaryIndex',1,'threshold',options.threshold*0.75);
    for imageNo = 1:size(croppedPeakImages,3)
        tmpImg = croppedPeakImages(:,:,imageNo);
        % tmpImg(boundaryIndices{1}) = tmpImg(boundaryIndices{1})+maxValMovie*options.crossHairPercent;
        tmpImg(boundaryIndices{1}) = NaN;
        croppedPeakImages(:,:,imageNo) = tmpImg;
    end
    % croppedPeakImages = cat(3,nanmean(croppedPeakImages(:,:,2:end),3),croppedPeakImages);
    if size(croppedPeakImages,3)<maxSignalsToShow
        dimDiff = maxSignalsToShow-size(croppedPeakImages,3);
        croppedPeakImagesTmp = NaN([size(croppedPeakImages,1) size(croppedPeakImages,2) dimDiff]);
        croppedPeakImages = cat(3,croppedPeakImages,croppedPeakImagesTmp);
    end
    croppedPeakImages2(:,:,:,1) = croppedPeakImages;
    warning off
    montage(permute(croppedPeakImages2(:,:,:,1),[1 2 4 3]))
    croppedPeakImages2 = getimage;
    % change zeros to ones, fixes range of image display
    croppedPeakImages2(croppedPeakImages2==0)=NaN;
    % croppedPeakImages2(1,1) = minValMovie;
    croppedPeakImages2(1,1) = 0;
    croppedPeakImages2(1,2) = maxValMovie*0.4;
    imagesc(croppedPeakImages2);
    colormap(options.colormap);
    axis off;
    % title('frames at signal peaks, press any key to exit');
    % ginput(1);
    % close(2);figure(mainFig);
    % clear croppedPeakImages2
    warning on
end

function [slopeRatio] = plotPeakSignal(thisTrace,testpeaks,cellIDStr,instructionStr,minValTraces,maxValTraces,peakROI,avgSpikeTrace,slopeRatio,spikeCenterTrace,valid)
    % display plots of the signal around peaks in the signal

    % [peakOutputStat] = computePeakStatistics(thisTrace,'waitbarOn',0);
    % avgSpikeTrace = peakOutputStat.avgSpikeTrace;
    % slopeRatio = peakOutputStat.slopeRatio;
    % spikeCenterTrace = peakOutputStat.spikeCenterTrace{1};

    peakSignalAmplitude = thisTrace(testpeaks(:));
    [peakSignalAmplitude peakIdx] = sort(spikeCenterTrace(:,round(end/2)+1),'descend');
    spikeCenterTrace = spikeCenterTrace(peakIdx,:);
    if size(spikeCenterTrace,1)>20
        spikeCenterTrace = spikeCenterTrace(1:20,:);
    end

    plot(repmat(peakROI, [size(spikeCenterTrace,1) 1])', spikeCenterTrace','Color',[4 4 4]/8)
    hold on;
    plot(peakROI, avgSpikeTrace,'k', 'LineWidth',3);box off;
    plot(peakROI, nanmean(spikeCenterTrace),'Color',[1 0 0 1.0], 'LineWidth',2);box off;
    % add in zero line
    xval = 0;
    x=[xval,xval];
    y=[minValTraces maxValTraces];
    plot(x,y,'r'); box off;

    hold off;
    % title(['signal transients ' cellIDStr])
    % title(['transients'])
    title(['signal ' cellIDStr 10 '(' num2str(sum(valid==1)) ' good)']);
    xlabel('frames');
    % ylabel('df/f');
    ylabel('\DeltaF/F');
    ylim([minValTraces maxValTraces]);
end

function plotSignal(thisTrace,testpeaks,cellIDStr,instructionStr,minValTraces,maxValTraces,options)
    % length(testpeaks)
    peakSignalAmplitude = thisTrace(testpeaks(:));
    [peakSignalAmplitude peakIdx] = sort(peakSignalAmplitude,'descend');
    testpeaks = testpeaks(peakIdx);
    if length(testpeaks)>options.maxSignalsToShow
        testpeaks = testpeaks(1:options.maxSignalsToShow);
    end

    % plots a signal along with test peaks
    plot(thisTrace, 'r');
    hold on;
    scatter(testpeaks, thisTrace(testpeaks), 'LineWidth',0.5,'MarkerFaceColor',[0 0 0], 'MarkerEdgeColor',[0 0 0])
    title([cellIDStr instructionStr])
    xlabel('frames');
    % ylabel('df/f');
    axis([0 length(thisTrace) minValTraces maxValTraces]);
    box off;
    hold off;
end

function [valid directionOfNextChoice saveData i] = respondToUserInput(reply,i,valid,directionOfNextChoice,saveData,nFilters)
    % decide what to do based on input (not a switch due to multiple comparisons)
    if isequal(reply, 3)|isequal(reply, 110)|isequal(reply, 31)
        % n key or right click
        directionOfNextChoice=1;
        % display('invalid IC');
        % set(mainFig,'Color',[0.8 0 0]);
        valid(i) = 0;
    elseif isequal(reply, 28)
        % go back, left
        directionOfNextChoice=-1;
    elseif isequal(reply, 29)
        % go forward, right
        directionOfNextChoice=1;
    elseif isequal(reply, 102)
        % user clicked 'f' for finished, exit loop
        movieDecision = questdlg('Are you sure you want to exit?', ...
            'Finish sorting', ...
            'yes','no','yes');
        if strcmp(movieDecision,'yes')
            saveData=1;
        end
        % i=nFilters+1;
    elseif isequal(reply, 103)
        % if user clicks 'g' for goto, ask for which IC they want to see
        icChange = inputdlg('enter signal #');
        if ~isempty(icChange)
            icChange = str2num(icChange{1});
            if icChange>nFilters|icChange<1
                % do nothing, invalid command
            else
                i = icChange;
                directionOfNextChoice = 0;
            end
        end
    elseif isequal(reply, 115)
        movieDecision = questdlg('Are you sure you want to exit?', ...
            'Finish sorting', ...
            'yes','no','yes');
        if strcmp(movieDecision,'yes')
            % 's' if user wants to get ride of the rest of the ICs
            display(['classifying the following signals as bad: ' num2str(i) ':' num2str(nFilters)])
            valid(i:nFilters) = 0;
            saveData=1;
        end
    elseif isequal(reply, 121)|isequal(reply, 1)|isequal(reply, 30)
        % y key or left click
        directionOfNextChoice=1;
        % display('valid IC');
        % set(mainFig,'Color',[0 0.8 0]);
        valid(i) = 1;
    else
        % forward=1;
        % valid(i) = 1;
    end
end