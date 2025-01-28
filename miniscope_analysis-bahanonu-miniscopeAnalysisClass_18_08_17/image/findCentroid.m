function [xCoords yCoords] = findCentroid(inputMatrix,varargin)
	% finds the x,y centroid coordinates of each 2D in the 3D input matrix
	% biafra ahanonu
	% started: 2013.10.31 [19:39:33]
	% adapted from SpikeE code
	% inputs
		%
	% outputs
		%
	% changelog
		% 2016.01.02 [21:22:13] - added comments and refactored so that accepts inputImages as [x y nSignals] instead of [nSignals x y]
		% 2016.08.06 - some changes to speed up algorithm by using thresholded rather than weighted sum of image.
		% 2017.01.14 [20:06:04] - support switched from [nSignals x y] to [x y nSignals]
	% TODO
		%

	%========================
	options.waitbarOn = 1;
	options.thresholdValue = 0.3;
	% threshold for images
	options.imageThreshold = 0.5;
	options.runImageThreshold = 1;
	% get options
	options = getOptions(options,varargin);
	% display(options)
	% unpack options into current workspace
	fn=fieldnames(options);
	for i=1:length(fn)
		eval([fn{i} '=options.' fn{i} ';']);
	end
	%========================

	inputDims = size(inputMatrix);
	inputDimsLen = length(inputDims);
	if inputDimsLen==3
	    nImages = size(inputMatrix,3);
	elseif inputDimsLen==2
	    nImages = 1;
	    tmpImage = inputMatrix; clear inputMatrix;
	    inputMatrix(:,:,1) = tmpImage;
	    options.waitbarOn = 0;
	else
	    return
	end

	if options.runImageThreshold==1
		inputMatrixThreshold = thresholdImages(inputMatrix,'waitbarOn',options.waitbarOn,'threshold',options.imageThreshold);
	else
		inputMatrixThreshold = inputMatrix;
	end

	reverseStr = '';
	for imageNum=1:nImages
		% threshold image
		thisImage = squeeze(inputMatrixThreshold(:,:,imageNum));
		% get the sum of the image
		imagesum = sum(thisImage(:));
		% get coordinates
		[i,j] = find(thisImage > options.thresholdValue);
		centroidPos = [mean(i(:)) mean(j(:))];
		if length(centroidPos)==1
			% size(i)
			% size(j)
			% centroidPos
			centroidPos = [centroidPos 1];
		end
		xCoords(imageNum) = round(centroidPos(2));
		yCoords(imageNum) = round(centroidPos(1));

		% clf;imagesc(thisImage);colorbar; hold on;
		% plot(xCoords(imageNum),yCoords(imageNum),'r+')
		% pause

		% xTmp = repmat(1:size(thisImage,2), size(thisImage,1), 1);
		% yTmp = repmat((1:size(thisImage,1))', 1,size(thisImage,2));
		% xCoords(imageNum) = sum(sum(thisImage.*xTmp))/imagesum;
		% yCoords(imageNum) = sum(sum(thisImage.*yTmp))/imagesum;
		% use median instead of mean?

		if (mod(imageNum,20)==0|imageNum==nImages)&options.waitbarOn==1
		    reverseStr = cmdWaitbar(imageNum,nImages,reverseStr,'inputStr','finding centroids');
		end
	end
end