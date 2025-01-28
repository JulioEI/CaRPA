function [initImages, initProbs] = initImagesFromICs(icImgs, icTraces, varargin)

options=getDefaultCELLMaxOptions();
options=getOptions(options,varargin);

% get sizes of all ICs
icSizes=zeros(size(icImgs,3),1);
for k=1:size(icImgs,3)
    icSizes(k)=sum(sum(icImgs(:,:,k)>0.3*max(max(icImgs(:,:,k)))));
end

% if using elimination heuristic, calculate max and set min size
if options.useInitialElimHeuristic
    numICsToUseForSizeCalc=ceil(size(icImgs,3)/5);
    maxSize=max(4*prctile(icSizes(1:numICsToUseForSizeCalc),50),100);
    minSize=8;
else
    maxSize=1000000;
    minSize=0;
end

% truncate ICs, eliminate negative values, and smooth slightly
for cInd=1:size(icImgs,3)
    thisImg=icImgs(:,:,cInd);
    maxVal=max(thisImg(:));
    thisSize=icSizes(cInd);
    if thisSize<maxSize && thisSize>minSize
        thisImg(thisImg<0.3*maxVal)=0.001*maxVal;
        icImgs(:,:,cInd)=thisImg;
    else
        icImgs(:,:,cInd)=nan;
    end
end 

initImages = permute(icImgs, [3 1 2]);
initImages = reshape(initImages, [size(initImages,1), numel(icImgs(:,:,1))]);
icTraces(isnan(initImages(:,1)),:)=[];
initImages(isnan(initImages(:,1)),:)=[];
initImages=bsxfun(@rdivide, initImages, sum(initImages,2));

initProbs=icTraces;
initProbs=bsxfun(@rdivide, initProbs, sum(initProbs,1));