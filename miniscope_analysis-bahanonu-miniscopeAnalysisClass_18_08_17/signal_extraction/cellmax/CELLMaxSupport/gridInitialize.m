function [cellFitParams, options] = gridInitialize(imgs,varargin)

	options.gridSpacing=6;
	options.gridWidth=8;
	options=getOptions(options,varargin);

	% xCoords=-(centroidSpacing/2):centroidSpacing:size(imgs,2)+(centroidSpacing/2);
	% yCoords=-(centroidSpacing/2):centroidSpacing:size(imgs,1)+(centroidSpacing/2);
	xCoords=0:options.gridSpacing:size(imgs,2);
	yCoords=0:options.gridSpacing:size(imgs,1);

	[xCentroids, yCentroids]=meshgrid(xCoords,yCoords);

	evenRows=1:2:size(xCentroids,1);

	xCentroids(evenRows,:)=xCentroids(evenRows,:)+options.gridSpacing/2;

	nGaussians=numel(xCentroids);
	cellFitParams=[xCentroids(:), yCentroids(:), options.gridWidth*ones(nGaussians,2), zeros(nGaussians,1)];
end