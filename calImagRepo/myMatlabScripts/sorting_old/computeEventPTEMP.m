probLog = [];
fail = 0;
for i = 1:size(inputSignals,1)
 movieDims = size(options.inputMovie);
    %Get event images
    [croppedPeakImages] = viewMontageCustom(options.inputMovie,inputImages(i,:,:),inputSignals(i,:,:),[signalPeakIdx{i}],minValMovie,maxValMovie,options.cropSizeLength,true);
    
    windowSize =10;
    %weightVectR = exp((1:windowSize)-1/(windowSize)-1)/(1.5*max(exp((1:windowSize)-1/(windowSize)-1)));
    weightDecay = 2;
    weightVectR = ((1:windowSize).^weightDecay)/max(1.2*((1:windowSize).^weightDecay));
    weightVect = [weightVectR,1,flip(weightVectR)];
    
    windowFrames = bsxfun(@plus,signalPeakIdx{i}',-windowSize:windowSize);
    windowFrames(windowFrames<=0)=1;
    windowFrames(windowFrames>=size(options.inputMovie,3))=size(options.inputMovie,3);
        
    allFramesLine = viewMontageCustom(options.inputMovie,inputImages(i,:,:),inputSignals(i,:,:),reshape(windowFrames,[1,numel(windowFrames)]),minValMovie,maxValMovie,options.cropSizeLength,false);

    allFrames = reshape(allFramesLine(:,:,2:end),[size(allFramesLine,1),size(allFramesLine,2),length(signalPeakIdx{i}),length(weightVect)]);
     
    eventsTimeAvg = ones(size(croppedPeakImages(:,:,2:end)));
    for k = 1:size(allFrames,3)
        averageMatTemp = repmat(reshape(weightVect,[1,1,numel(weightVect)]),[size(allFrames,1),size(allFrames,2),1]).*squeeze(allFrames(:,:,k,:));
        averageMat = max(averageMatTemp(:))*(sum(averageMatTemp,3)/max(max(sum(averageMatTemp,3))));
        eventsTimeAvg(:,:,k) = averageMat;
    end
    
    PSF = fspecial('gaussian',ceil(size(croppedPeakImages,1)/5),5);
    %Get cell mask, center and diameters
    
    cellImg = croppedPeakImages(:,:,1);
    cellImgF = filter2(PSF,cellImg);
    cellMask = cellImgF>nanmean(cellImgF(:))+1*nanstd(cellImgF(:));   

        
    %Compute mask means of events
    signalImagesCrop = eventsTimeAvg;%croppedPeakImages(:,:,2:end);

    %signalImagesCrop = cat(3,signalImagesCrop,nanmean(croppedPeakImages(:,:,2:end),3));
    meanValVect = zeros([1,size(signalImagesCrop,3)]);
    meanOutVect = zeros([1,size(signalImagesCrop,3)]); 
    for j = 1:size(signalImagesCrop,3)
        tempEvent = signalImagesCrop(:,:,j);
        meanValVect(j) = nanmean(tempEvent(cellMask));
        meanOutVect(j) = nanmean(tempEvent(~cellMask));
    end
    [meanValVect,sortingIdx] = sort(meanValVect,'descend');
    signalImagesCrop = signalImagesCrop(:,:,sortingIdx);

    signalImagesCropF = zeros(size(signalImagesCrop));
    for j = 1:size(signalImagesCrop,3)
        signalImagesCropF(:,:,j) = filter2(PSF,signalImagesCrop(:,:,j));
    end
     
    
    tresholdedImg = zeros(size(signalImagesCrop));
    tresholdedOutside = zeros(size(signalImagesCrop));
    overlap = zeros([1,size(signalImagesCropF,3)]);
    for k = 1:size(signalImagesCropF,3)
        tempImg = signalImagesCropF(:,:,k);
        tresholdedImg(:,:,k) = meanValVect(k)*((tempImg>nanmean(tempImg(:))+1*nanstd(tempImg(:)))).*(cellImgF>nanmean(cellImgF(:))+1*nanstd(cellImgF(:)));
        tresholdedOutside(:,:,k) = (tempImg>nanmean(tempImg(:))+1*nanstd(tempImg(:))) .* (1-tresholdedImg(:,:,k));
        missed = sum(sum(0~=(cellImgF>nanmean(cellImgF(:))+1*nanstd(cellImgF(:))))) - sum(sum(0~=tresholdedImg(:,:,k)));
        overlap(k) = sum(sum(0~=tresholdedImg(:,:,k)))/(missed+sum(sum(0~=tresholdedOutside(:,:,k))));
        %overlap(k) = 1./(1 + exp(-30*((overlap(k)-0.5)-0)));
        probCell = overlap(k);%*(meanValVect(k)/max(meanValVect));
%         if probCell< 0; probCell = 0; end
        probLog = [probLog,probCell];
    end

end