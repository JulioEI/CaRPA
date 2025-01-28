function myData2 = getCellInfo2(inputMovie,inputImages,inputSignals,signalPeakIdx)

%Computes features for each cell

%inputMovie = options.inputMovie;
minValMovie = min(inputMovie(:));
maxValMovie = max(inputMovie(:));
cropSizeLength = 10;

 i = 2;
 
inputImage = inputImages(i,:,:);
inputSignal = inputSignals(i,:,:);
signalPeakId = signalPeakIdx{i};

%--Get and process images--%

%Get cropped images at peak value
[croppedPeakImages] = viewMontageCustom(inputMovie,inputImage,inputSignal,signalPeakId,minValMovie,maxValMovie,cropSizeLength,true);

%Get cell image, a filtered and a tresholded version
cellImg = croppedPeakImages(:,:,1);
PSF = fspecial('gaussian',ceil(size(croppedPeakImages,1)/5),5);
cellImgF = filter2(PSF,cellImg);
cellMask = cellImgF>nanmean(cellImgF(:))+1*nanstd(cellImgF(:));

%--Get cell propierties--%

%Cell smoothness
cellInfo.sScore(i) = sum(sum(abs(del2(cellImg))))/7.36; %7.36 is mean noise value

%Distance from center
cellMap = squeeze(inputImage);
[x,y] = find(cellMap==max(cellMap(:)));
cellMapCenter = size(cellMap)/2;
cellInfo.centroidDist(i) = sqrt((x-cellMapCenter(1))^2 + (y-cellMapCenter(2))^2)/sqrt((size(cellMap,1)-cellMapCenter(1))^2 + (size(cellMap,2)-cellMapCenter(2))^2);

%Global SNR
function [globalSNR] = getGlobalSNR(signalPeakId,inputSignal)
timeSeq = -3:3;
peakIdx = bsxfun(@plus,timeSeq',signalPeakId);
peakIdx = unique(peakIdx(:));
peakIdx(peakIdx>length(inputSignal))=[];
peakIdx(peakIdx<=0)=[]; 
noiseSignal = inputSignal;
noiseSignal(peakIdx) = [];
globalSNR = nanmean(inputSignal(signalPeakId)/nanstd(noiseSignal));

% cellInfo.globalSNR(i) = nanmean(inputSignal(signalPeakId)/nanstd(noiseSignal));
% 
% %RMS
% cellInfo.RMS(i) = mean(abs(inputSignal(peakIdx)));
% 
% %peakSNR
% spikeROI = -20:20;
% extractMatrix = bsxfun(@plus,signalPeakId',spikeROI);
% extractMatrix(extractMatrix<=0)=1;
% extractMatrix(extractMatrix>=size(inputSignal,2))=size(inputSignal,2);
% spikeCenterTrace = reshape(inputSignal(extractMatrix),size(extractMatrix));
% cellInfo.peakSNR(i) = mean(std(spikeCenterTrace,1))/nanstd(noiseSignal);
        
function [croppedPeakImages] = viewMontageCustom(inputMovie,inputImage,thisTrace,signalPeakArray,minValMovie,maxValMovie,cropSizeLength,doSort)

    if isempty(signalPeakArray)
        imagesc(inputImage);
        colormap(customColormap([]));
        axis off;
        croppedPeakImages2 = inputImage;
        return
    end
    % signalPeakArray
    maxSignalsToShow = 1e16;
    if doSort
        peakSignalAmplitude = thisTrace(signalPeakArray(:));
    % peakSignalAmplitude
        [peakSignalAmplitude, peakIdx] = sort(peakSignalAmplitude,'descend');
        % peakSignalAmplitude
        signalPeakArray = signalPeakArray(peakIdx);
    end
    if length(signalPeakArray)>maxSignalsToShow
        % choose a random subset
        signalPeakArray = signalPeakArray(1:maxSignalsToShow);
    end

%     signalPeakArray((signalPeakArray-31)<1) = [];
%     signalPeakArray((signalPeakArray+31)>length(thisTrace)) = [];

    signalPeakArray = {signalPeakArray};
    % signalPeakArray
    croppedPeakImages = compareSignalToMovieCustom(inputMovie, inputImage, thisTrace,'getOnlyPeakImages',1,'waitbarOn',0,'extendedCrosshairs',0,'signalPeakArray',signalPeakArray,'cropSize',cropSizeLength,'crosshairs',0);
function [croppedPeakImages] = compareSignalToMovieCustom(inputMovie, inputFilters, inputSignal, varargin)
	% shows a cropped version of inputMovie for each inputFilters and aligns it to inputSignal peaks to make sure detection is working
	% biafra ahanonu
	% started: 2013.11.04 [18:40:45]
	% inputs
		% inputMovie - matrix dims are [X Y t] - where t = number of time points
		% inputFilters - matrix dims are [n X Y] - where n = number of filters, NOTE THE DIFFERENCE
		% inputSignal - matrix dims are [n t] - where n = number of signals, t = number of time points
	% outputs
		% none, this is a display function
	% changelog
		% 2014.01.18 [12:24:29] fully implemented, cut out from controllerAnalysis, need to improve handling at beginning of movie, but that's a playMovie function issue
	% TODO
		%

	%========================
	% old way of saving, only temporary until full switch
	options.oldSave = 0;
	% size in pixels to show signal image
	options.cropSize = 20;
	% frames before/after to show
	options.timeSeq = -10:10;
	% waitbar
	options.waitbarOn = 1;
	% whether to just get the peak images and ignore showing the movie
	options.getOnlyPeakImages = 0;
	% 1 = plus shaped crosshairs, 0 = dot
	options.extendedCrosshairs = 1;
	%
	options.crosshairs = 1;
	%
	options.signalPeakArray = [];
	% set to 1 if input images should be normalized
	options.normalizeMovieImages = 0;
	% get options
	options = getOptions(options,varargin);
	% unpack options into current workspace
	% fn=fieldnames(options);
	% for i=1:length(fn)
	%     eval([fn{i} '=options.' fn{i} ';']);
	% end
	%========================

	if isempty(options.signalPeakArray)
		[signalPeaks, signalPeakArray] = computeSignalPeaks(inputSignal, 'makePlots', 0,'makeSummaryPlots',0,'waitbarOn',options.waitbarOn);
	else
		signalPeakArray = options.signalPeakArray;
	end

	% get the centroids and other info for movie
	[xCoords yCoords] = findCentroid(inputFilters,'waitbarOn',options.waitbarOn);
	cropSize = options.cropSize;
	nSignals = size(inputFilters,1);
	nPoints = size(inputMovie,3);
	movieDims = size(inputMovie);
	timeSeq = options.timeSeq;

	% inputMovie(inputMovie>1.3) = NaN;
	% inputMovie(inputMovie<0.8) = NaN;

	% loop over all signals and visualize their peaks side-by-side with movie
	exitSignal = 0;
	for signalNo=1:nSignals
		peakLocations = signalPeakArray{signalNo};

		peakIdxs = bsxfun(@plus,timeSeq',peakLocations);
		peakIdxs(find(peakIdxs<1)) = 1;
		peakIdxs(find(peakIdxs>nPoints)) = 1;
		% get region to crop
		warning off;
		xLow = xCoords(signalNo) - cropSize;
		xHigh = xCoords(signalNo) + cropSize;
		yLow = yCoords(signalNo) - cropSize;
		yHigh = yCoords(signalNo) + cropSize;
		% check that not outside movie dimensions
		xMin = 0;
		xMax = movieDims(2);
		yMin = 0;
		yMax = movieDims(1);

		% adjust for the difference in centroid location if movie is cropped
		xDiff = 0;
		yDiff = 0;
		if xLow<=xMin; xDiff = xLow-xMin; xLow = xMin+1; end
		if xHigh>=xMax; xDiff = xHigh-xMax; xHigh = xMax-1; end
		if yLow<=yMin; yDiff = yLow-yMin; yLow = yMin+1; end
		if yHigh>=yMax; yDiff = yHigh-yMax; yHigh = yMax-1; end

		% need to add a way to adjust the cropped movie target point if near the boundary

		% get the cropped movie at peaks
		% yLow
		% yHigh
		% xLow
		% xHigh
		% peakLocations
		croppedPeakImages = inputMovie(yLow:yHigh,xLow:xHigh,peakLocations);
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        globalCoord = [xCoords(signalNo),yCoords(signalNo)];
        %croppedPeakImages = shiftEvents(croppedPeakImages,inputMovie,peakLocations,cropSize,movieDims,globalCoord);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
        firstImg = squeeze(inputFilters(signalNo,yLow:yHigh,xLow:xHigh));
		firstImg = normalizeVector(firstImg,'normRange','zeroToOne');
		maxVec = nanmax(croppedPeakImages(:));
		minVec = nanmin(croppedPeakImages(:));
		firstImg = firstImg*maxVec;
		% firstImg = (firstImg-minVec)./(maxVec-minVec);
		% if options.normalizeMovieImages==1
		% end
		croppedPeakImages(:,:,end+1) = firstImg;
		% move inputImage to the front
		croppedPeakImages = circshift(croppedPeakImages,[0 0 1]);
		% croppedPeakImagesTmp = croppedPeakImages(:,:,end);
		% croppedPeakImagesTmp(:,:,end+1:end+(length(croppedPeakImages)-1)) = croppedPeakImages(:,:,1:(end-1));
		% croppedPeakImages = croppedPeakImagesTmp;
		for frameNo=1:size(croppedPeakImages,3)
			cropImg = squeeze(croppedPeakImages(:,:,frameNo));
			if options.normalizeMovieImages==1
				cropImg = normalizeVector(cropImg,'normRange','zeroToOne');
			end
			croppedPeakImages(:,:,frameNo) = cropImg;
		end
		% croppedPeakImages = normalizeMovie(croppedPeakImages,'normalizationType','meanDivision');
		cDims = size(croppedPeakImages);
		if options.getOnlyPeakImages==0
			% get cropped version of the movie
			croppedMovie = inputMovie(yLow:yHigh,xLow:xHigh,peakIdxs);
			cDims = size(croppedMovie);
			exitSignal = playMovie(inputMovie(:,:,peakIdxs),'extraMovie',croppedMovie,...
				'extraLinePlot',inputSignal(signalNo,peakIdxs),...
				'windowLength',30,...
				'colormapColor','jet',...
				'extraTitleText',['signal #' num2str(signalNo) '/' num2str(nSignals) '    peaks: ' num2str(length(peakLocations))],...
				'primaryPoint',[xCoords(signalNo) yCoords(signalNo)],...
				'secondaryPoint',crossHairLocation);
				% 'recordMovie','test.avi',...
		end
		warning on;
		if exitSignal==1
			break;
		end
    end

