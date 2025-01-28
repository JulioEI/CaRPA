function [success] = unitTestDownsampleMovie(inputMovie,varargin)
	% Tests downsampling of a 3D matrix using Matlab imresize, interp1, and imageJ then automatically compares the results
	% building off of unit_matlab_imageJ [biafra ahanonu]
	% biafra ahanonu
	% started: 2014.07.14
	% inputs
		% inputMovie - 3D matrix [x y t], char string with path to movie file, or cell array of strings each pointing to a movie file to be loaded (they will all be concatenated together)
	% outputs
		%

	% changelog
		% 2016.02.xx - Added additional comments and tests, made more user friendly.
	% TODO
		% Refactor so that each downsample vector is stored in a named fieldname within a structure and make plotting functions dynamic, would allow easier addition of new downsample methods.

	%========================
	% number of frames to use for max projection calculation
	options.nFramesMaxProject = 500;
	% downsample factor of 4
	options.downsampleFactor = 4;
	% inputMovie HDF5 datasetname
	options.inputDatasetName = '/1';
	% 'nearest', 'next', 'previous', 'linear','spline','pchip', or 'cubic'
	options.interp1Method = 'linear';
	% number of peaks to display in final plot
	options.nPeaksDisplay = 4;
	% size of crop area to keep in movie after user selects point
	options.cropSize = 10;
	% Filter order number for decimate()
	options.decimateFilterOrder = 8;
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
		success = 0;

		% check that Miji exists, if not, have user enter information
		if exist('Miji.m','file')==2
			display(['Miji located in: ' which('Miji.m')]);
			% Miji is loaded, continue
		else
			pathToMiji = uigetdir('\.','Enter path to Miji.m in Fiji (e.g. \Fiji.app\scripts)');

			addpath(pathToMiji);
			% privateLoadBatchFxnsPath = 'private\privateLoadBatchFxns.m';
			% fid = fopen(privateLoadBatchFxnsPath,'at')
			% % fprintf(fid, '\npathtoMiji = ''%s'';\n', pathToMiji);
			% fclose(fid);
		end

		% load the movie
		inputMovieClass = class(inputMovie);
		if strcmp(inputMovieClass,'char')
		    inputMovie = loadMovieList(inputMovie,'inputDatasetName',options.inputDatasetName);
		end
		inputMovie = single(inputMovie);

		% get the max projection of the first user specified number of frames
		if size(inputMovie,3)<options.nFramesMaxProject
			options.nFramesMaxProject = size(inputMovie,3);
		end
		inputMovieMaxProject = max(dfofMovie(inputMovie(:,:,1:options.nFramesMaxProject)),[],3);

		% ask the user to select a pixel for use in downsampling comparison
		[~,~] = openFigure(1, '');
		imagesc(inputMovieMaxProject);
		title('Select a pixel');
		userCoords = round(ginput(1));
		% userCoords
		rowCoord = userCoords(1);
		colCoord = userCoords(2);

		% crop the movie to save computation time and dfof
		cropSize = options.cropSize;
		subfxnCropInputMovie();
		inputMovie = dfofMovie(inputMovie);

		startTime = tic;
		%% =====================
		% grab the 1D vector corresponding to that point in space
		lineInputMovieRaw = squeeze(inputMovie(rowCoord,colCoord,:));
		subsampleIdx = floor(linspace(1,length(lineInputMovieRaw),length(lineInputMovieRaw)/options.downsampleFactor));
		lineInputMovieRawSubsampled = lineInputMovieRaw(subsampleIdx);

		%% =====================
		% get matlab downsample
		subfxnImresizeDownsample()

		%% =====================
		% get interp1 based downsampling
		subfxnInterp1DecimateDownsample()

		%% =====================
		subfxnImageJDownsample()

		%% =====================
		% get peaks and setup plot legends
		[signalPeaks, signalPeakIdx] = computeSignalPeaks(double(lineInputMovieImresizeBilinear'),'makePlots', 0,'makeSummaryPlots',0,'waitbarOn',1);
		signalPeakIdx = signalPeakIdx{1};
		% signalPeakIdx
		% display(signalPeakIdx{1})
		[xPlot yPlot] = getSubplotDimensions(2*options.nPeaksDisplay+1);
		legendStr = {...
		['raw | ' num2str(nanvar(lineInputMovieRaw(:)))],...
		['imresize | bilinear | ' num2str(nanvar(lineInputMovieImresizeBilinear(:)))],...
		['imresize | bicubic | '  num2str(nanvar(lineInputMovieImresizeBicubic(:)))],...
		['interp1 | ' options.interp1Method ' | ' num2str(nanvar(lineInputMovieInterp1Down(:)))],...
		['decimate | Chebyshev | ' num2str(nanvar(lineInputMovieDecimateDown(:)))],...
		['imagej | avg | '  num2str(nanvar(lineInputMovieImageJDownAvg(:)))],...
		['imagej | no avg | '  num2str(nanvar(lineInputMovieImageJDownNoAvg(:)))]};

		%% =====================
		% normalize all plots with dfof
		% lineInputMovieRaw = normalizeVector(lineInputMovieRaw,'normRange','dfof');
		% lineInputMovieImresizeBilinear = normalizeVector(lineInputMovieImresizeBilinear,'normRange','dfof');
		% lineInputMovieImresizeBicubic = normalizeVector(lineInputMovieImresizeBicubic,'normRange','dfof');
		% lineInputMovieInterp1Down = normalizeVector(lineInputMovieInterp1Down,'normRange','dfof');
		% lineInputMovieDecimateDown = normalizeVector(lineInputMovieDecimateDown,'normRange','dfof');
		% lineInputMovieImageJDownAvg = normalizeVector(lineInputMovieImageJDownAvg,'normRange','dfof');
		% lineInputMovieImageJDownNoAvg = normalizeVector(lineInputMovieImageJDownNoAvg,'normRange','dfof');

		%% =====================
		% plot data

		% plot whole trial on single plot
		[~,~] = openFigure(2, '');
			clf
			hold off;
			inputVector{1} = lineInputMovieRaw;
			inputVector{end+1} = lineInputMovieImresizeBilinear;
			inputVector{end+1} = lineInputMovieImresizeBicubic;
			inputVector{end+1} = lineInputMovieInterp1Down;
			inputVector{end+1} = lineInputMovieDecimateDown;
			inputVector{end+1} = lineInputMovieImageJDownAvg;
			inputVector{end+1} = lineInputMovieImageJDownNoAvg;
			input1X = linspace(1,length(inputVector{end}),length(inputVector{1}));
			subfxnPlotVectors(inputVector,input1X,legendStr);

		[~,~] = openFigure(3, '');
			clf
			% signalPeakIdx = signalPeakIdx(randperm(length(signalPeakIdx)));
		% plot each downsampled trace cropped to peaks
		for subplotNo = 0:2:2*options.nPeaksDisplay
			% subplotNo
			if subplotNo==0
				subplot(xPlot,yPlot,1)
				subfxnPlotVectors(inputVector,input1X,legendStr);
			else
				try
					cutRange = {[-20:20]',[-3:3]'};
					cutRangeRaw = {[(-20*options.downsampleFactor):(20*options.downsampleFactor)]',[(-3*options.downsampleFactor):(3*options.downsampleFactor)]'};
					nameStr = {'zoom out', 'close-up'};
					for plotType = 1:2
						hold off;inputVector = {};
						subplot(xPlot,yPlot,subplotNo+(plotType-1))

						% get raw vector
						signalIdx = round(options.downsampleFactor*(signalPeakIdx(subplotNo/2)-0.5));
						signalIdx = bsxfun(@plus,signalIdx,cutRangeRaw{plotType});
						inputVector{1} = lineInputMovieRaw(signalIdx);
						input1X = linspace(1,length(cutRange{plotType}),length(inputVector{1}));

						% get downsampled vectors
						signalIdx = signalPeakIdx(subplotNo/2);
						signalIdx = bsxfun(@plus,signalIdx,cutRange{plotType});
						inputVector{end+1} = lineInputMovieImresizeBilinear(signalIdx);
						inputVector{end+1} = lineInputMovieImresizeBicubic(signalIdx);
						inputVector{end+1} = lineInputMovieInterp1Down(signalIdx);
						inputVector{end+1} = lineInputMovieDecimateDown(signalIdx);
						inputVector{end+1} = lineInputMovieImageJDownAvg(signalIdx);
						inputVector{end+1} = lineInputMovieImageJDownNoAvg(signalIdx);

						% !
						subfxnPlotVectors(inputVector,input1X,'');
						title(['peak #' num2str(subplotNo/2) ' ' nameStr{plotType}])
						xlabel('frames');ylabel('dfof');
					end
				catch
				end
			end
		end
		suptitle('Downsample methods comparison | raw line is sampled at native rate')

		success = 1;

	catch err
		success = 0;
		display(repmat('@',1,7))
		disp(getReport(err,'extended','hyperlinks','on'));
		display(repmat('@',1,7))
	end
	function subfxnCropInputMovie()
		movieDims = size(inputMovie);
		xLow = rowCoord - cropSize;
		xHigh = rowCoord + cropSize;
		yLow = colCoord - cropSize;
		yHigh = colCoord + cropSize;
		% check that not outside movie dimensions
		xMin = 1;
		xMax = movieDims(2);
		yMin = 1;
		yMax = movieDims(1);

		% adjust for the difference in centroid location if movie is cropped
		xDiff = 0;
		yDiff = 0;
		if xLow<xMin xDiff = xLow-xMin; xLow = xMin; end
		if xHigh>xMax xDiff = xHigh-xMax; xHigh = xMax; end
		if yLow<yMin yDiff = yLow-yMin; yLow = yMin; end
		if yHigh>yMax yDiff = yHigh-yMax; yHigh = yMax; end

		% crop movie and readjust coords
		inputMovie = inputMovie(yLow:yHigh,xLow:xHigh,:);
		cDims = size(inputMovie);
		crossHairLocation = [round(cDims(2)/2+xDiff/2) round(cDims(1)/2+yDiff/2)];
		rowCoord = crossHairLocation(1);
		colCoord = crossHairLocation(2);
	end
	function subfxnImresizeDownsample()
		startTime = tic;
		tic
		inputMovieMatlabDown = downsampleMovie(inputMovie,'downsampleFactor',options.downsampleFactor,'downsampleDimension','time','downsampleType','bilinear');
		toc
		timeInputMovieImresizeBilinear = toc(startTime);
		lineInputMovieImresizeBilinear = squeeze(inputMovieMatlabDown(rowCoord,colCoord,:));

		startTime = tic;
		tic
		inputMovieMatlabDown = downsampleMovie(inputMovie,'downsampleFactor',options.downsampleFactor,'downsampleDimension','time','downsampleType','bicubic');
		lineInputMovieImresizeBicubic = squeeze(inputMovieMatlabDown(rowCoord,colCoord,:));
		toc
		timeInputMovieImresizeBicubic = toc(startTime);
		% lineInputMovieImresizeDown
	end
	function subfxnInterp1DecimateDownsample()
		display('interp1/decimate downsampling')
		inputMovieInterp1Down = zeros(size(inputMovie));
		inputMovieDecimateDown = zeros(size(inputMovie));
		filterOrder = options.decimateFilterOrder;
		for rowNo = 1:size(inputMovie,1)
			for colNo = 1:size(inputMovie,2)
				lineInputMovieRawTmp = double(squeeze(inputMovie(rowNo,colNo,:)));
				lineInputMovieInterp1Down = interp1(1:length(lineInputMovieRaw),lineInputMovieRaw,linspace(1,length(lineInputMovieRaw),length(lineInputMovieRaw)/options.downsampleFactor),options.interp1Method);
				inputMovieInterp1Down(rowNo,colNo,1:length(lineInputMovieInterp1Down)) = lineInputMovieInterp1Down;

				lineInputMovieDecimateDown = decimate(lineInputMovieRawTmp,options.downsampleFactor,filterOrder);
				inputMovieDecimateDown(rowNo,colNo,1:length(lineInputMovieDecimateDown)) = lineInputMovieDecimateDown;
			end
		end
		inputMovieInterp1Down = inputMovieInterp1Down(:,:,1:length(lineInputMovieInterp1Down));
		lineInputMovieInterp1Down = squeeze(inputMovieInterp1Down(rowCoord,colCoord,:));
		inputMovieDecimateDown = inputMovieDecimateDown(:,:,1:length(lineInputMovieDecimateDown));
		lineInputMovieDecimateDown = squeeze(inputMovieDecimateDown(rowCoord,colCoord,:));
	end
	function subfxnImageJDownsample()
		% get imagej downsampled
		Miji;
		MIJ.createImage('imagej', inputMovie, true);
		% run imagej downsampling with averaging
		MIJ.run('Scale...', ['x=1.0 y=1.0 z=' num2str(1/options.downsampleFactor) ' width=' num2str(size(inputMovie,1)) ' height=' num2str(size(inputMovie,2)) ' depth=' num2str(floor(size(inputMovie,3)/options.downsampleFactor)) ' interpolation=Bilinear average process create title=imagej_down']);
		inputMovieImageJDownAvg = MIJ.getCurrentImage;
		MIJ.run('Close');
		% run imagej downsampling without averaging
		MIJ.run('Scale...', ['x=1.0 y=1.0 z=' num2str(1/options.downsampleFactor) ' width=' num2str(size(inputMovie,1)) ' height=' num2str(size(inputMovie,2)) ' depth=' num2str(floor(size(inputMovie,3)/options.downsampleFactor)) ' interpolation=Bilinear process create title=imagej_down']);
		inputMovieImageJDownNoAvg = MIJ.getCurrentImage;
		MIJ.run('Close');MIJ.run('Close');
		MIJ.exit

		lineInputMovieImageJDownAvg = squeeze(inputMovieImageJDownAvg(rowCoord,colCoord,:));
		lineInputMovieImageJDownNoAvg = squeeze(inputMovieImageJDownNoAvg(rowCoord,colCoord,:));
	end
end
function subfxnPlotVectors(inputVector,input1X,legendStr)
	plot(input1X,inputVector{1},'Color',[0.5 0.5 0.5],'LineWidth',5);
	hold on;
	plot(inputVector{2},'r','LineWidth',3);
	plot(inputVector{3},'k','LineWidth',1);
	plot(inputVector{4},'g','LineWidth',1);
	plot(inputVector{5},'y','LineWidth',1);
	plot(inputVector{6},'b','LineWidth',1);
	plot(inputVector{7},'c','LineWidth',1);
	box off;axis tight;
	yL = get(gca,'YLim');
	xL = get(gca,'XLim');
	line([xL(2)/2 xL(2)/2],yL,'Color','r','LineWidth',1);
	if ~isempty(legendStr)
		legend(legendStr);
	end
end
function subfxnCompareHistograms()
	% % send over matlab downsampled and look at histogram of diff from imagej
	% MIJ.createImage('matlab_down', testMovieMatlabDown, true);
	% MIJ.run('Image Calculator...','image1=matlab_down operation=Subtract image2=imagej_down create 32-bit stack');
	% MIJ.run('Histogram', 'bins=256 use x_min=0 x_max=0 y_max=Auto stack');
	% uiwait(msgbox('press OK to end downsample test','Success','modal'));
	% MIJ.run('Close');MIJ.run('Close');MIJ.run('Close');MIJ.run('Close');MIJ.run('Close');
	% MIJ.exit
end