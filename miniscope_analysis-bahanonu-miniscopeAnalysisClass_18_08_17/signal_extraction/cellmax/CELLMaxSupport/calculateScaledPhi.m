function scaledPhi = calculateScaledPhi(cellImages, dsScaledPhi, DFOF, CELLMaxOptions,varargin)
	% example function with outline for necessary components
	% Lacey Kitsch
	% updated by: biafra ahanonu
	% started: xxxx.xx.xx
	% inputs
	    %
	% outputs
	    %

	% changelog
	    % 2017.04.xx - added support for changing iteration number and ensure proper calc of mfilt - biafra
	    % 2017.04.xx - change phi estimate upsampling function to imresize bilinear
	% TODO
	    %

	%========================
	% number of iterations to run to get new phi estimate
	options.nIterations = 20;
	% bilinear or bicubic
	options.downsampleType = 'bilinear';

	options = getOptions(options, varargin);
	%========================

	nFrames=size(DFOF,3);
	nFramesDS=size(dsScaledPhi,2);
	nCells=size(cellImages,3);

	cellImages=reshape(cellImages, [numel(cellImages(:,:,1)) size(cellImages,3)]);
	cellImages=cellImages';

	if nFramesDS<nFrames
		nFramesInt = round(nFrames/nFramesDS);
		fprintf('Ratio of downsampled to upsampled is %d\n',nFramesInt);
		% linearly interpolate new phi estimate
	    phiEst = imresize(dsScaledPhi,[size(dsScaledPhi,1) nFrames],options.downsampleType);

	    % old way using filter
	    % fm=mfilt.firsrc(nFramesInt,1);
	    % phiEst=zeros(nCells,nFrames);
	    % for cInd=1:nCells
	    %     phiEst(cInd,:)=filter(fm,dsScaledPhi(cInd,:));
	    % end
	    phiEst(phiEst<0)=0;
	    phiEst=bsxfun(@rdivide, phiEst, sum(phiEst,1));
	else
	    phiEst=dsScaledPhi;
	end

	for iterNum = 1:options.nIterations
	    [~,phiEst,~,scaledPhi] = EM_genFilt_oneIteration(DFOF, cellImages, phiEst, 'options', CELLMaxOptions);
	end
	scaledPhi=double(scaledPhi);

	% nCells=size(cellImages,3);
	% cellImages=reshape(cellImages, [numel(cellImages(:,:,1)) size(cellImages,3)]);
	% imgsInverse=(cellImages'*cellImages)\cellImages';
	%
	% fInc=options.noiseSigma/5;
	% thresh=options.numSigmasThresh*options.noiseSigma;
	%
	% nFrames=size(DFOF,3);
	% frInc=200;
	% scaledPhi=zeros(nCells,nFrames);
	% for fr=1:frInc:nFrames
	%     frameInds=fr:(fr+199);
	%     frameInds(frameInds>nFrames)=[];
	%
	%     thisCounts=DFOF(:,:,frameInds);
	%     if max(thisCounts(:))>1
	%         thisCounts=thisCounts-1;
	%     end
	%     thisCounts(thisCounts<=thresh)=0;
	%     thisCounts=floor(thisCounts/fInc);
	%
	%     thisCounts=reshape(thisCounts, [numel(DFOF(:,:,1)) length(frameInds)]);
	%     totCounts=sum(thisCounts,1);
	%
	%     scaledPhi(:,frameInds)=imgsInverse*thisCounts;
	%
	% end
	% scaledPhi=scaledPhi*fInc;
end