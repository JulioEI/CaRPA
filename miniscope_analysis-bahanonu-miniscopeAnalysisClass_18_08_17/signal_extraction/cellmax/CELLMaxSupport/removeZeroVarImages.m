function [cellImages, goodInds] = removeZeroVarImages(cellImages)

	% 2017.03.27 - display progress, etc. - biafra

	nImages=size(cellImages,3);
	imgsToDelete=zeros(nImages,1);
	reverseStr = '';
	for k=1:nImages
		if mod(k,25)==0
		    reverseStr = cmdWaitbar(k,nImages,reverseStr,'inputStr','remove zero variance images','waitbarOn',1,'displayEvery',25);
		end

	    thisImage=cellImages(:,:,k);
	    thisImage(thisImage==0)=max(thisImage(:));
	    if max(thisImage(:))==min(thisImage(:))
	        imgsToDelete(k)=1;
	    end
	end
	imgsToDelete=logical(imgsToDelete);
	cellImages(:,:,imgsToDelete)=[];
	goodInds=1:nImages;
	goodInds(imgsToDelete)=[];
end