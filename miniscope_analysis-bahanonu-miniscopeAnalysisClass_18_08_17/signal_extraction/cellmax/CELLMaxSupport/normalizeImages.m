function normImgs = normalizeImages(origImgs,varargin)

options.threshFrac=0;
options=getOptions(options,varargin);

normImgs=zeros(size(origImgs));

for imgInd=1:size(origImgs,3)
    thisImg=origImgs(:,:,imgInd);
    maxVal=max(thisImg(:));
    if maxVal>0
        normImgs(:,:,imgInd)=thisImg/maxVal;
    end
    if options.threshFrac>0
        thisNormImg=normImgs(:,:,imgInd);
        thisNormImg(thisImg<options.threshFrac*maxVal)=0;
        normImgs(:,:,imgInd)=thisNormImg;
    end
end
