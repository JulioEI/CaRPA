function imageCentroids = getImageCentroids(cImages, imgSize, varargin)

options.thresh=0.2;
options=getOptions(options,varargin);

% if cImages is 2-dimensional, we assume the first dimension is image
% index, and the rows are the images flattened
if numel(size(cImages))==2
    if prod(imgSize)==size(cImages,1)
        cImages=cImages';
    end
    cImages=reshape(cImages, [size(cImages,1) imgSize]);
    cImages=permute(cImages, [2 3 1]);
end

nImages=size(cImages,3);
imageCentroids=zeros(nImages,2);
for cInd=1:nImages
    
    thisImage=cImages(:,:,cInd);
    
    [maxVal,maxLoc]=max(thisImage(:));

    thisImage(thisImage<=options.thresh*maxVal)=0;
    thisImage(thisImage>0)=1;
    
    thisImagelabel=bwlabel(logical(thisImage),4);
    bwIC=thisImagelabel==thisImagelabel(maxLoc);
    %imageSize=sum(bwIC(:));
    
    info=regionprops(bwIC,'Centroid');
    imageCentroids(cInd,:)=info.Centroid;
end