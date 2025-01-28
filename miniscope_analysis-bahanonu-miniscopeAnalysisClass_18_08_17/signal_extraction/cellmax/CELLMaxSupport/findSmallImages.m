function smallCellInds = findSmallImages(cellImages, imgSize, sizeThresh)

% if cellImages is 2-dimensional, we assume the first dimension is image
% index, and the rows are the images flattened
if numel(size(cellImages))==2
    if prod(imgSize)==size(cellImages,1)
        cellImages=cellImages';
    end
    cellImages=reshape(cellImages, [size(cellImages,1) imgSize]);
    cellImages=permute(cellImages, [2 3 1]);
end

nImages=size(cellImages,3);
imageAreas=zeros(nImages,1);
for cInd=1:nImages
    
    thisImage=cellImages(:,:,cInd);
    
    [maxVal,maxLoc]=max(thisImage(:));

    thisImage(thisImage<=0.2*maxVal)=0;
    thisImage(thisImage>0)=1;
    
    thisImagelabel=bwlabel(logical(thisImage),4);
    bwImg=thisImagelabel==thisImagelabel(maxLoc);
    %imageSize=sum(bwIC(:));
    
    info=regionprops(bwImg,'Centroid', 'Area');
    imageAreas(cInd)=info.Area;
end

smallCellInds=find(imageAreas<sizeThresh);
