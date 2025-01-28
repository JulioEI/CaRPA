function images = removeDiscontig(images)

nImages=size(images,3);
fractionMaxValThresh=0.07;

reverseStr = '';
for k=1:nImages
    thisImage=images(:,:,k);
    [maxVal,maxLoc]=max(thisImage(:));
    thisImageBin=images(:,:,k);
    thisImageBin(thisImageBin<=(fractionMaxValThresh*maxVal))=0;
    thisImageBin(thisImageBin>0)=1;
    thisImageBin=logical(thisImageBin);
    imageLabels=bwlabel(thisImageBin,8);
    thisImage(imageLabels~=imageLabels(maxLoc))=0;
    images(:,:,k)=thisImage;

    if mod(k,25)==0
    	reverseStr = cmdWaitbar(k,nImages,reverseStr,'inputStr','removing discontiguous regions','waitbarOn',1,'displayEvery',25);
    end
end