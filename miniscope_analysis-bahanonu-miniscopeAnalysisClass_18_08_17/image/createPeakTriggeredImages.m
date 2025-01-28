function [outputImages outputMeanImageCorrs] = createPeakTriggeredImages(inputMovie, inputImages, inputSignals, varargin)
	% Gets event triggered average image from an input movie based on cell images located in input image and trace matrix
	% biafra ahanonu
	% started: 2015.09.28, abstracted from behaviorAnalysis
	% inputs
		% inputMovie - [x y frames]
		% inputImages - [x y nSignals]
		% inputSignals - [nSignals frames]
	% outputs
		%

	% changelog
		% 2016.01.02 [21:22:13] - added comments and refactored so that accepts inputImages as [x y nSignals] instead of [nSignals x y]
	% TODO
		%

	%========================
	% size in pixels to crop around the movie
	options.cropSize = 10;
	% hierarchy name in hdf5 where movie is
	options.inputDatasetName = '/1';
	% save time if already computed peaks
	% options.signalPeaks = [];
	options.signalPeaksArray = [];
	% show waitbar or not
	options.waitbarOn = 1;
	%
	options.normalizeOutput = 1;
	% get options
	options = getOptions(options,varargin);
	% display(options)
	% unpack options into current workspace
	% fn=fieldnames(options);
	% for i=1:length(fn)
	% 	eval([fn{i} '=options.' fn{i} ';']);
	% end
	%========================

	try
		% load movie
		if ischar(inputMovie)
			inputMovie = loadMovieList(inputMovie,'inputDatasetName',options.inputDatasetName);
		else

		end

		if isempty(options.signalPeaksArray)
			[~, signalPeaksArray] = computeSignalPeaks(inputSignals,'waitbarOn',options.waitbarOn,'makeSummaryPlots',0);
		else
			signalPeaksArray = options.signalPeaksArray;
		end

		[xCoords yCoords] = findCentroid(inputImages,'waitbarOn',options.waitbarOn);
		xCoords = round(xCoords);
		yCoords = round(yCoords);

		nSignals = size(inputSignals,1)
		movieDims = size(inputMovie);
		outputImages = NaN(size(inputImages));

		outputMeanImageCorrs = [];

		reverseStr = '';
		% loop over all signals
		for signalNo = 1:nSignals
			% get frames when cell files into [x y frames] matrix
			% signalNo
			signalImages = inputMovie(:,:,signalPeaksArray{signalNo});

			xLow = xCoords(signalNo) - options.cropSize;
			xHigh = xCoords(signalNo) + options.cropSize;
			yLow = yCoords(signalNo) - options.cropSize;
			yHigh = yCoords(signalNo) + options.cropSize;
			% check that not outside movie dimensions
			xMin = 1 ;
			xMax = movieDims(2);
			yMin = 1 ;
			yMax = movieDims(1);

			% adjust for the difference in centroid location if movie is cropped
			xDiff = 0;
			yDiff = 0;
			if xLow<xMin xDiff = xLow-xMin; xLow = xMin; end
			if xHigh>xMax xDiff = xHigh-xMax; xHigh = xMax; end
			if yLow<yMin yDiff = yLow-yMin; yLow = yMin; end
			if yHigh>yMax yDiff = yHigh-yMax; yHigh = yMax; end

			% [yLow yHigh xLow xHigh]
			% crop to a region of x-y pixels around the cell's centroid
			signalImagesCrop = signalImages(yLow:yHigh,xLow:xHigh,:);

			corrVals = [];
			signalImagesCropTmp = signalImagesCrop;
			signalImagesCropTmp(isnan(signalImagesCropTmp)) = 0;
			inputImageCrop = squeeze(inputImages(yLow:yHigh,xLow:xHigh,signalNo));
			for peakImageNo = 1:size(signalImagesCropTmp,3)
				corrVals(peakImageNo) = corr2(inputImageCrop,squeeze(signalImagesCropTmp(:,:,peakImageNo)));
			end
			% corrVals
			outputMeanImageCorrs(signalNo) = nanmean(corrVals);

			% take the average of all frames in matrix
			signalImagesCropSingle = squeeze(nanmean(signalImagesCrop,3));

			signalImage = NaN(size(squeeze(inputImages(:,:,1))));
			% size(signalImage)
			signalImage(yLow:yHigh,xLow:xHigh) = signalImagesCropSingle;

			% imagesc(signalImage);
			% store resulting image in
			% size(outputImages)
			outputImages(:,:,signalNo) = signalImage;


			reverseStr = cmdWaitbar(signalNo,nSignals,reverseStr,'inputStr','getting movie images','waitbarOn',1,'displayEvery',50);
		end

		if options.normalizeOutput==1
			outputImages = normalizeMovie(outputImages,'normalizationType','zeroToOne');
		end

		%
	catch err
		outputImages = [];
		outputMeanImageCorrs = [];
		display(repmat('@',1,7))
		disp(getReport(err,'extended','hyperlinks','on'));
		display(repmat('@',1,7))
	end
end
