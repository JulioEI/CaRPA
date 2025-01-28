function [sqSizeX, sqSizeY, nSqHorz, nSqVert, sqIncX, sqIncY] = findOptimalSquareSize(maxSqSize,desOverlap,imgSize)
	% calculates the chunk size based on movie statistics

	% nSquares=1:round(max(size(imgs,1),size(imgs,2))/(maxSqSize/3));
	% imgsSizeY = size(imgs,1);
	% imgsSizeX = size(imgs,2);
	imgsSizeY = imgSize(1);
	imgsSizeX = imgSize(2);

	nSquares=1:round(max(imgsSizeY,imgsSizeX)/(maxSqSize/3));
	optSqSizesX=zeros(size(nSquares));
	optSqSizesY=zeros(size(nSquares));
	% imgSize=size(imgs(:,:,1));
	imgSize = [imgsSizeY imgsSizeX];

	for nSqInd=1:length(nSquares)
	    nSq=nSquares(nSqInd);
	    optSqSizesX(nSqInd)=(imgSize(2)+(nSq-1)*desOverlap)/nSq;
	    optSqSizesY(nSqInd)=(imgSize(1)+(nSq-1)*desOverlap)/nSq;
	end

	sqSizeX=round(optSqSizesX(find(optSqSizesX<=(maxSqSize+1), 1, 'first')));
	sqSizeY=round(optSqSizesY(find(optSqSizesY<=(maxSqSize+1), 1, 'first')));
	nSqHorz=nSquares(find(optSqSizesX<=(maxSqSize+1), 1, 'first'));
	nSqVert=nSquares(find(optSqSizesY<=(maxSqSize+1), 1, 'first'));

	sqIncX=sqSizeX-desOverlap;
	sqIncY=sqSizeY-desOverlap;
end