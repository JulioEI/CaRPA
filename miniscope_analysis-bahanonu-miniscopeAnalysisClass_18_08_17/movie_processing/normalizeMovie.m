function [inputMovie] = normalizeMovie(inputMovie, varargin)
	% takes an input movie and applies a particular normalization (e.g. lowpass divisive).
	% biafra ahanonu
	% started: 2013.11.09 [09:25:48]
	% inputs
		% inputMovie = [x y frames] 3D matrix
	% outputs
		% inputMovie = [x y frames] 3D matrix normalized

	% changelog
		% 2014.02.17 added in mean subtraction/division to function.
	% TODO
		%

	% input is an image, convert to movie
	if length(size(inputMovie))==2
		inputMovieDims = 2;
		inputMovieTmp(:,:,1) = inputMovie;
		inputMovie = inputMovieTmp;
	else
		inputMovieDims = 3;
	end
	%========================
	% fft,bandpassDivisive,lowpassFFTDivisive,imfilterSmooth,imfilter,meanSubtraction,meanDivision,negativeRemoval
	options.normalizationType = 'meanDivision';
	% for fft
	options.secondaryNormalizationType = [];
	% maximum frame to normalize
	options.maxFrame = size(inputMovie,3);
	% use parallel registration (using matlab pool)
	options.parallel = 1;
	% ===
	% options for fft
	% for bandpass, low freq to pass
	options.freqLow = 10;
	% for bandpass, high freq to pass
	options.freqHigh = 50;
	% imageJ normalization options
	options.imagejFFTLarge = 10000;
	options.imagejFFTSmall = 80;
	% highpass, lowpass, bandpass
	options.bandpassType = 'highpass';
	% binary or gaussian
	options.bandpassMask = 'gaussian';
	% show the frequency spectrum and images
	% 0 = no, 1 = yes
	options.showImages = 0;
	% ===
	% fspecial, 'disk' option
	option.imfilterType = 'disk';
	% how to deal with boundaries, see http://www.mathworks.com/help/images/ref/imfilter.html
	options.boundaryType = 'circular';
	% 'disk' option: pixel radius to blur
	options.blurRadius = 35;
	% 'gaussian' option
	options.sizeBlur = 80;
	options.sigmaBlur = 3;
	% ===
	% cmd line waitbar on?
	options.waitbarOn = 1;
	% get options
	options = getOptions(options,varargin);
	% unpack options into current workspace
	fn=fieldnames(options);
	for i=1:length(fn)
	    eval([fn{i} '=options.' fn{i} ';']);
	end
	if strcmp(normalizationType,'bandpassDivisive')
		normalizationType = 'fft';
	end
	if strcmp(normalizationType,'lowpassFFTDivisive')
		normalizationType = 'fft';
		options.secondaryNormalizationType = 'lowpassFFTDivisive';
		options.bandpassType = 'lowpass';
	end
	%========================
	switch normalizationType
		case 'imagejFFT'
			% opens imagej
			mijiCheck()
			% MUST ADD \Fiji.app\scripts
			% open imagej instance
			% Miji(false);
			Miji;
			startTime = tic;
			% pass matrix to imagej
			MIJ.createImage('result', inputMovie, true);
			% settings taken from original imagej implementation
			% bpstr= ' filter_large=10000 filter_small=80 suppress=None tolerance=5 process';
			bpstr= [' filter_large=' num2str(options.imagejFFTLarge) ' filter_small=' num2str(options.imagejFFTSmall) ' suppress=None tolerance=5 process'];
			MIJ.run('Bandpass Filter...',bpstr);
			% grab the image from imagej
			inputMovieFFT = MIJ.getCurrentImage;
			% close imagej instance
			MIJ.run('Close');
			MIJ.exit;
			toc(startTime);
			% divide lowpass from image
			inputMovie = bsxfun(@rdivide,single(inputMovie),single(inputMovieFFT));
		case 'imagejFFT_test'
			mijiCheck()
			reverseStr = '';
			% inputImage = squeeze(inputMovie(:,:,1));
			options.runfftTest=0;
			lowFreqList = [10 50 80];
			highFreqList = [10000 5000 120];
			fontSize = 30;
			userFreqList = inputdlg({'filter_small','filter_large','font size (pt)'},'bandpass parameters',[1 100],{num2str(lowFreqList),num2str(highFreqList),num2str(fontSize)});
			lowFreqList = str2num(userFreqList{1});
			highFreqList = str2num(userFreqList{2});
			fontSize = str2num(userFreqList{3});

			[lowFreqList,highFreqList] = meshgrid(lowFreqList, highFreqList);
			lowFreqList = lowFreqList(:);
			highFreqList = highFreqList(:);
			nFreqs = length(lowFreqList);
			% pairs = [p(:) q(:)];

			% inputImageTest = zeros([size(inputImage,1) size(inputImage,2) nFreqs]);
			inputImageTest = {};
			inputMovieDuplicate = {};
			inputMovieDivide = {};
			% size(inputImageTest)
			% Miji(false);
			Miji;
			startTime = tic;
			for freqNo = 1:nFreqs
				lowFreq = lowFreqList(freqNo);
				highFreq = highFreqList(freqNo);
				% pass matrix to imagej
				MIJ.createImage('result', inputMovie, true);
				% settings taken from original imagej implementation
				bpstr= [' filter_large=' num2str(highFreq) ' filter_small=' num2str(lowFreq) ' suppress=None tolerance=5 process'];
				MIJ.run('Bandpass Filter...',bpstr);
				% grab the image from imagej
				inputMovieFFT = MIJ.getCurrentImage;
				% inputImageTest(:,:,freqNo) = inputMovieFFT;
				inputImageTest{freqNo} = inputMovieFFT;
				inputMovieDuplicate{freqNo} = inputMovie;
				% divide lowpass from image
				inputMovieDivide{freqNo} = bsxfun(@rdivide,single(inputMovie),single(inputMovieFFT));

				moptions.identifyingText{freqNo} = [num2str(lowFreq) ' | ' num2str(highFreq)];
				reverseStr = cmdWaitbar(freqNo,nFreqs,reverseStr,'inputStr','normalizing movie','displayEvery',5);
				MIJ.run('Close');
			end
			MIJ.exit;
			% moptions.identifyingText = strsplit(num2str(freqList),' ');
			moptions.singleRowMontage = 1;
			moptions.fontSize = fontSize;
			[inputMovieDuplicate] = createMontageMovie(inputMovieDuplicate,'options',moptions);
			moptions.identifyingText = [];
			[inputImageTest] = createMontageMovie(inputImageTest,'options',moptions);
			[inputMovieDivide] = createMontageMovie(inputMovieDivide,'options',moptions);
			inputMovie = permute(cat(2,inputMovieDuplicate,inputImageTest,inputMovieDivide),[2 1 3]);
			% close imagej instance
			% [inputImageTest] = addText(inputImageTest,freqList,42);
			% inputImageTestArray(:,:,:,1) = inputImageTest;
			% figure(10)
			% montage(permute(inputImageTestArray(:,:,:,1),[1 2 4 3]))
			% inputMovie = inputImageTestArray;
		case 'fft'
			display('Running Matlab FFT')
			bandpassMatrix = zeros(size(inputMovie));
			% get options
			ioptions.showImages = options.showImages;
			ioptions.lowFreq = options.freqLow;
			ioptions.highFreq = options.freqHigh;
			ioptions.bandpassType = options.bandpassType;
			ioptions.bandpassMask = options.bandpassMask;
			ioptions.padImage = 1;
			% convert movie to correct class output by fft
			outputClass = class(fftImage(squeeze(inputMovie(:,:,1)),'options',ioptions));
			inputMovie = cast(inputMovie,outputClass);
			% pre-calculate filter to save time
			testImage = squeeze(inputMovie(:,:,1));
			[cutoffFilter] = createCutoffFilter(testImage,ioptions.bandpassMask,ioptions.lowFreq,ioptions.highFreq,ioptions.bandpassType,ioptions.padImage);
			ioptions.cutoffFilter = cutoffFilter;
			figure;imagesc(cutoffFilter);title([num2str(ioptions.lowFreq) ' | ' num2str(ioptions.highFreq)])
			% ========================
			manageParallelWorkers('parallel',options.parallel);
			% ========================
			%Get dimension information about 3D movie matrix
			[inputMovieX inputMovieY inputMovieZ] = size(inputMovie);
			reshapeValue = size(inputMovie);
			%Convert array to cell array, allows slicing (not contiguous memory block)
			inputMovie = squeeze(mat2cell(inputMovie,inputMovieX,inputMovieY,ones(1,inputMovieZ)));

			reverseStr = '';
			% parfor_progress(options.maxFrame);
			% dispstat('','init');
			parfor frame=1:options.maxFrame
				% thisFrame = squeeze(inputMovie(:,:,frame));
				% if isempty(options.secondaryNormalizationType)
				% 	inputMovie(:,:,frame) = fftImage(thisFrame,'options',ioptions);
				% else
				% 	tmpFrame = fftImage(thisFrame,'options',ioptions);
				% 	inputMovie(:,:,frame) = thisFrame./tmpFrame;
				% end
				% reverseStr = cmdWaitbar(frame,options.maxFrame,reverseStr,'inputStr','normalizing movie','waitbarOn',options.waitbarOn,'displayEvery',5);
				thisFrame = squeeze(inputMovie{frame});
				if isempty(options.secondaryNormalizationType)
					inputMovie{frame} = fftImage(thisFrame,'options',ioptions);
				else
					tmpFrame = fftImage(thisFrame,'options',ioptions);
					inputMovie{frame} = thisFrame./tmpFrame;
				end
				% if mod(frame,round(options.maxFrame/10))==0
				% 	fprintf ('%2.0f-',frame/options.maxFrame*100);drawnow
				% end
				% percent = parfor_progress;
				% dispstat(num2str(percent),'keepthis');
				% parfor_progress
				% drawnow
				% bandpassMatrix(:,:,frame) = fftImage(thisFrame,'options',ioptions);
				% bandpassMatrix(:,:,frame) = imcomplement(bandpassMatrix(:,:,frame));
				% = bsxfun(@ldivide,squeeze(movie20hz(:,:,1)),filteredFrame
			end
			% parfor_progress(0);
			inputMovie = cat(3,inputMovie{:});

			% options.secondaryNormalizationType
			% if isempty(options.secondaryNormalizationType)
			% 	inputMovie = bsxfun(@ldivide,inputMovie,inputMovieTmp);
			% end
			% inputMovie = bandpassMatrix;
		case 'matlabFFT_test'
			reverseStr = '';
			% inputImage = squeeze(inputMovie(:,:,1));
			% options.runfftTest=0;
			% lowFreqList = [10 50 80];
			% highFreqList = [50 100 500];
			lowFreqList = [1 3];
			highFreqList = [7 14];
			fontSize = 30;
			userFreqList = inputdlg({'lowFreq','highFreq','font size (pt)'},'bandpass parameters',[1 100],{num2str(lowFreqList),num2str(highFreqList),num2str(fontSize)});
			lowFreqList = str2num(userFreqList{1});
			highFreqList = str2num(userFreqList{2});
			fontSize = str2num(userFreqList{3});

			[lowFreqList,highFreqList] = meshgrid(lowFreqList, highFreqList);
			lowFreqList = lowFreqList(:);
			highFreqList = highFreqList(:);
			nFreqs = length(lowFreqList);
			% pairs = [p(:) q(:)];

			% inputImageTest = zeros([size(inputImage,1) size(inputImage,2) nFreqs]);
			inputImageTest = {};
			inputMovieDuplicate = {};
			inputMovieDivide = {};
			% size(inputImageTest)


			startTime = tic;
			for freqNo = 1:nFreqs
				display([num2str(freqNo) '/' num2str(nFreqs)])
				lowFreq = lowFreqList(freqNo);
				highFreq = highFreqList(freqNo);
				inputMovieFFT = inputMovie;
				% set FFT options
				ioptions.showImages = options.showImages;
				ioptions.lowFreq = lowFreq;
				ioptions.highFreq = highFreq;
				ioptions.bandpassType = options.bandpassType;
				ioptions.bandpassMask = options.bandpassMask;
				ioptions.padImage=1;
				% covert to correct output class
				outputClass = class(fftImage(squeeze(inputMovie(:,:,1)),'options',ioptions));
				inputMovieFFT = cast(inputMovieFFT,outputClass);
				% ============
				% % pre-calculate filter to save time
				testImage = squeeze(inputMovieFFT(:,:,1));
				[cutoffFilter] = createCutoffFilter(testImage,ioptions.bandpassMask,ioptions.lowFreq,ioptions.highFreq,ioptions.bandpassType,ioptions.padImage);
				ioptions.cutoffFilter = cutoffFilter;
				figure;imagesc(cutoffFilter);title([num2str(lowFreq) ' | ' num2str(highFreq)])
				% FFT each image
				reverseStr = '';
				options.secondaryNormalizationType
				nanmean(inputMovieFFT(:))
				for frame=1:options.maxFrame
					thisFrame = squeeze(inputMovieFFT(:,:,frame));
					if isempty(options.secondaryNormalizationType)
						inputMovieFFT(:,:,frame) = fftImage(thisFrame,'options',ioptions);
					else
						% tmpFrame = fftImage(thisFrame,'options',ioptions);
						% inputMovieFFT(:,:,frame) = thisFrame./tmpFrame;
						inputMovieFFT(:,:,frame) = fftImage(thisFrame,'options',ioptions);
					end
					% bandpassMatrix(:,:,frame) = fftImage(thisFrame,'options',ioptions);
					% bandpassMatrix(:,:,frame) = imcomplement(bandpassMatrix(:,:,frame));
					reverseStr = cmdWaitbar(frame,options.maxFrame,reverseStr,'inputStr','normalizing movie','waitbarOn',options.waitbarOn,'displayEvery',5);
					% = bsxfun(@ldivide,squeeze(movie20hz(:,:,1)),filteredFrame
				end
				% inputImageTest(:,:,freqNo) = inputMovieFFT;
				inputImageTest{freqNo} = inputMovieFFT;
				inputMovieDuplicate{freqNo} = cast(inputMovie,outputClass);
				% divide lowpass from image
				inputMovieDivide{freqNo} = bsxfun(@rdivide,inputMovieDuplicate{freqNo},inputImageTest{freqNo});
				% if isempty(options.secondaryNormalizationType)
				% else
				% 	inputMovieDivide{freqNo} = inputImageTest{freqNo};
				% end
				% nanmean(inputMovieDuplicate{freqNo}(:))
				% nanmean(inputImageTest{freqNo}(:))
				% inputMovieDivide{freqNo} = bsxfun(@minus,inputMovieDuplicate{freqNo},inputImageTest{freqNo});
				inputMovieDfof{freqNo} = dfofMovie(inputMovieDivide{freqNo});

				moptions.identifyingText{freqNo} = [num2str(lowFreq) ' | ' num2str(highFreq)];
				% reverseStr = cmdWaitbar(freqNo,nFreqs,reverseStr,'inputStr','normalizing movie','displayEvery',5);
			end
			% moptions.identifyingText = strsplit(num2str(freqList),' ');
			moptions.singleRowMontage = 1;
			moptions.fontSize = fontSize;
			[inputMovieDuplicate] = createMontageMovie(inputMovieDuplicate,'options',moptions);
			moptions.identifyingText = [];
			[inputImageTest] = createMontageMovie(inputImageTest,'options',moptions);
			[inputMovieDivide] = createMontageMovie(inputMovieDivide,'options',moptions);
			[inputMovieDfof] = createMontageMovie(inputMovieDfof,'options',moptions);
			inputMovie = permute(cat(2,inputMovieDuplicate,inputImageTest,inputMovieDivide,inputMovieDfof),[2 1 3]);
		case 'imfilterSmooth'
			% create filter
			switch option.imfilterType
				case 'disk'
					movieFilter = fspecial('disk', options.blurRadius);
				case 'gaussian'
					movieFilter = fspecial('gaussian', [options.sizeBlur options.sizeBlur], options.sigmaBlur);
				otherwise
					return
			end
			nFrames = size(inputMovie,3);
			inputMovieFiltered = zeros(size(inputMovie));
			reverseStr = '';
			for frame=1:nFrames
			    thisFrame = squeeze(inputMovie(:,:,frame));
			    thisFrame(find(isnan(thisFrame)))=nanmean(thisFrame(:));
			    inputMovieFiltered(:,:,frame) = imfilter(thisFrame, movieFilter,options.boundaryType);
			    reverseStr = cmdWaitbar(frame,nFrames,reverseStr,'inputStr','normalizing movie','waitbarOn',options.waitbarOn,'displayEvery',5);
			end
			% divide each frame by the filtered movie to remove 'background'
			inputMovie = inputMovieFiltered;
		case 'imfilter'
			% create filter
			switch option.imfilterType
				case 'disk'
					movieFilter = fspecial('disk', options.blurRadius);
				case 'gaussian'
					movieFilter = fspecial('gaussian', [options.sizeBlur options.sizeBlur], options.sigmaBlur);
				otherwise
					return
			end
			nFrames = size(inputMovie,3);
			inputMovieFiltered = zeros(size(inputMovie));
			reverseStr = '';
			for frame=1:nFrames
			    thisFrame = squeeze(inputMovie(:,:,frame));
			    thisFrame(find(isnan(thisFrame)))=nanmean(thisFrame(:));
			    inputMovieFiltered(:,:,frame) = imfilter(thisFrame, movieFilter,options.boundaryType);
			    reverseStr = cmdWaitbar(frame,nFrames,reverseStr,'inputStr','normalizing movie','waitbarOn',options.waitbarOn,'displayEvery',5);
			end
			% divide each frame by the filtered movie to remove 'background'
			inputMovie = bsxfun(@ldivide,inputMovieFiltered,inputMovie);
		case 'meanSubtraction'
			inputMean = nanmean(nanmean(inputMovie,1),2);
			inputMean = cast(inputMean,class(inputMovie));
			inputMovie = bsxfun(@minus,inputMovie,inputMean);
		case 'meanDivision'
			inputMean = nanmean(nanmean(inputMovie,1),2);
			inputMean = cast(inputMean,class(inputMovie));
			inputMovie = bsxfun(@rdivide,inputMovie,inputMean);
			% inputMean = nansum(nansum(inputMovie,1),2);
			% inputMean = cast(inputMean,class(inputMovie))
		case 'negativeRemoval'
			inputMin = abs(nanmin(inputMovie(:)))
			inputMin = cast(inputMin,class(inputMovie));
			inputMovie = bsxfun(@plus,inputMovie,inputMin);
		case 'zeroToOne'
			nFrames = size(inputMovie,3);
			reverseStr = '';
			for frame=1:nFrames
			    thisFrame = squeeze(inputMovie(:,:,frame));
				maxVec = nanmax(thisFrame(:));
				minVec = nanmin(thisFrame(:));
				meanVec = nanmean(thisFrame(:));
				inputMovie(:,:,frame) = (thisFrame-minVec)./(maxVec-minVec);
			    reverseStr = cmdWaitbar(frame,nFrames,reverseStr,'inputStr','normalizing movie','waitbarOn',options.waitbarOn,'displayEvery',5);
			end
		otherwise
			inputMovie = NaN
			return;
	end

	if inputMovieDims==2
		inputMovie = squeeze(inputMovie(:,:,1));
	end
end
function mijiCheck()
	if exist('Miji.m','file')==2
		display(['Miji located in: ' which('Miji.m')]);
		% Miji is loaded, continue
	else
		pathToMiji = inputdlg('Enter path to Miji.m in Fiji (e.g. \Fiji.app\scripts):',...
		             'Miji path', [1 100]);
		pathToMiji = pathToMiji{1};
		privateLoadBatchFxnsPath = 'private\privateLoadBatchFxns.m';
		fid = fopen(privateLoadBatchFxnsPath,'at')
		fprintf(fid, '\npathtoMiji = ''%s'';\n', pathToMiji);
		fclose(fid);
	end
end
function [movieTmp] = addText(movieTmp,inputText,fontSize)
	nFrames = size(movieTmp,3);
	maxVal = nanmax(movieTmp(:));
	minVal = nanmin(movieTmp(:));
	reverseStr = '';
	for frameNo = 1:nFrames
		movieTmp(:,:,frameNo) = squeeze(nanmean(...
			insertText(movieTmp(:,:,frameNo),[0 0],num2str(inputText(frameNo)),...
			'BoxColor',[maxVal maxVal maxVal],...
			'TextColor',[minVal minVal minVal],...
			'AnchorPoint','LeftTop',...
			'FontSize',fontSize,...
			'BoxOpacity',1)...
		,3));
		reverseStr = cmdWaitbar(frameNo,nFrames,reverseStr,'inputStr','adding text to movie','waitbarOn',1,'displayEvery',10);
	end
	% maxVal = nanmax(movieTmp(:))
	% movieTmp(movieTmp==maxVal) = 1;
	% 'BoxColor','white'
end
function [cutoffFilter] = createCutoffFilter(testImage,bandpassMask,lowFreq,highFreq,bandpassType,padImage)
	cutoffFilter = [];
	if padImage==1
		padSize = round(1.0*mean(size(testImage)));
		testImage = padarray(testImage,[padSize padSize],'symmetric');
	end
	testImageFFT = fft2(testImage);
	testImageFFT = fftshift(testImageFFT);
	[imFFTX imFFTY] = size(testImageFFT);

	switch bandpassMask
		case 'gaussian'
			% implemented using fspecial
			if lowFreq==0
				highpassFilter = ones([imFFTX imFFTY]);
			else
				highpassFilter = 1-normalizeVector(fspecial('gaussian', [imFFTX imFFTY],lowFreq),'normRange','zeroToOne');
			end
			lowpassFilter = normalizeVector(fspecial('gaussian', [imFFTX imFFTY],highFreq),'normRange','zeroToOne');
			switch bandpassType
				case 'highpass'
					cutoffFilter = highpassFilter;
				case 'lowpass'
					cutoffFilter = lowpassFilter;
				case 'bandpass'
					cutoffFilter = highpassFilter.*lowpassFilter;
				otherwise
					% do nothing
			end
		case 'binary'
			% create binary mask, tried with fspecial but this is easier
			[ffty fftx] = size(testImageFFT);
			cx = round(fftx/2);
			cy = round(ffty/2);
			[x,y] = meshgrid(-(cx-1):(fftx-cx),-(cy-1):(ffty-cy));
			highpassFilter = ((x.^2+y.^2)>lowFreq^2);
			lowpassFilter = ((x.^2+y.^2)<highFreq^2);
			switch bandpassType
				case 'highpass'
					cutoffFilter = highpassFilter;
				case 'lowpass'
					cutoffFilter = lowpassFilter;
				case 'bandpass'
					cutoffFilter = highpassFilter.*lowpassFilter;
				otherwise
					% do nothing
			end
		otherwise
			display('invalid option given')
			filtered_image = inputImage;
			return
	end
end