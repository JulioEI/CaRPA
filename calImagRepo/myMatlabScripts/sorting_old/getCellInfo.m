function [cellInfo] = getCellInfo(inputSignals,testpeaksArray,cellImages,inputMovie, varargin)
%Gathers various parameters about traces, shapes and events for clustering
%and classification.
%
%	cellInfo = getCellInfo(inputSignals,testpeaksArray,cellImages,inputMovie);
%
%You can specify a set of features by specifying its name plus get, e.g:
%
%   cellInfo = getCellInfo(inputSignals,testpeaksArray,cellImages,inputMovie,'getexpRatio','getglobalSNR');
%
%You can show the related plots by including 'showPlots':
%
%   cellInfo = getCellInfo(inputSignals,testpeaksArray,cellImages,inputMovie,'showPlots','getexpRatio','getglobalSNR');
%
%FEATURE INFO:
%
%Traces:
%   expRatio: Ratio of exponential fit of the curves expRatio:
%       tau_post/tau_post_pre   where tau_post  comes from the fitted exponential curve to the trace f(t)= a*e^(-t/tau_post)
%       This gives an idea of the speed of the dynamics of the trace. tau_post_post/tau_post_pre >1  is desired
%   expPeak: Average of trace peak over all the events.
%   globalSNR: SNR peaks vs signal:
%   peakSNR: SNR peak variability:
%   RMS: rms of signal
%
%Shapes:
%   cellEcc: Eccentricity of thresholded cell:
%       Ratio of the distance between the foci of the ellipse and its major axis length.
% 	cellSol: Solidity of thresholded cell:
% 		Proportion of the pixels in the convex hull that are also in the region.
%   pScore: Relation between perimeter and perimeter of the fitted ellipse:
% 		pScore= abs(1-realPerimeter/ellipseFitedPerimter)
% 		This gives an idea of how smooth is the contour of the cell. Ideally pScore=0
% 	cScore: Similitude to circle:
%       Computed as the std of the distances of the perimeter points to the centroid. Is 0 if the shape is a sphere.
%   sScore: Smoothnes of the shape, computed as the sum of the absolute
%   values of the laplacian of the shape.
%   cellPer: Perimeter of the cell
%
%   centroidDist: Distance of the cell centroid from the center of the cellMap.
%
% Events:
% 		imMovCorr:
% 			Correlation between movie and image 
%       getinterCorr:
%           Correlation between movie events.
%       eventEcc
%       eventSol
%       eventP
%       eventC
%       eventPer

%Initialize options:
params.showPlots = false;
params.getexpRatio = true;
params.getexpPeak = true;
params.getglobalSNR = true;
params.getpeakSNR = true;
params.getcellEcc = true;
params.getcellSol = true;
params.getpScore = true;
params.getcScore = true;
params.getsScore = true;
params.getcellPer = true;
params.getRMS = true;
params.getimMovCorr = true;
params.getinterSNR = true;
params.getinterCorr = true;
params.getcentroidDist = true;
params.geteventEcc = true;
params.geteventSol = true;
params.geteventP = true;
params.geteventC = true;
params.geteventS = true;
params.geteventPer = true;
params.geteventDist = true;

paramsFields = fieldnames(params);

%If the user selects certain feat., only compute those.
paramFeatStart = 2; %When the features start in the params
slectFeat = true;
for i = 1:numel(varargin)
    %paramI = ~cellfun(@isempty,strfind(paramsFields,varargin{i}));
    paramI = strcmp(paramsFields,varargin{i});
    
    if slectFeat && sum(paramI(paramFeatStart:end))
        for k = paramFeatStart:size(paramsFields)
            params.(paramsFields{k}) = false;
        end
        slectFeat = false;
    end
    try
        params.(paramsFields{paramI}) = true;
    catch
        error('Feature requested not found')
    end
end

%Initializes struct for desired features
nCells = size(inputSignals,1);
for k = paramFeatStart:size(paramsFields)
    if params.(paramsFields{k})
        paramsName = paramsFields{k};
        cellInfo.(paramsName(4:end)) =  NaN([nCells,1]);
    end
end

%Loops over all cells to extract the selected features.
for i=1:nCells
    %disp([num2str(100*i/nCells),'%'])
    
    %TRACES
    
    if params.getglobalSNR || params.getRMS || params.getpeakSNR
        timeSeq = -3:3;
        spikeROI = -20:20;
		if ~isempty(testpeaksArray{i})
			[globalSNR,RMS,peakSNR] = snrRmsFun(inputSignals(i,:),testpeaksArray{i}, timeSeq,spikeROI , params);
			if params.getglobalSNR
				cellInfo.globalSNR(i) = globalSNR;
			end
			if params.getRMS
				cellInfo.RMS(i) = RMS;
			end
			if params.getpeakSNR
				cellInfo.peakSNR(i) = peakSNR;
			end
		end
    end
    if params.getexpPeak
		if ~isempty(testpeaksArray{i})
			cellInfo.expPeak(i) = nanmean(inputSignals(i,testpeaksArray{i}));
			if params.showPlots
			   disp(['expPeak: ',num2str(cellInfo.expPeak(i))])
			end
		end
    end    
    if params.getexpRatio
        spikeROI = -40:40;
        slopeFrameWindow = 10;
		if ~isempty(testpeaksArray{i})
			[expRatio] = traceFun(inputSignals(i,:),testpeaksArray{i},spikeROI,slopeFrameWindow,params);
			if params.getexpRatio
				cellInfo.expRatio(i) = expRatio;
			end
		end
    end
    %SHAPES
    if  params.getsScore
        cellInfo.sScore(i) = sum(sum(abs(del2(cellImages(:,:,i)))));
        if params.showPlots
            disp(['Cell smooth score: ',num2str(cellInfo.sScore(i))])
        end
    end
    if params.getcellEcc || params.getcellSol || params.getpScore || params.getcScore || params.getcellPer || params.getcentroidDist 
    	cellMap = cellImages(:,:,i);
        [cellEcc,cellSol,pScore,cScore,cellPer,centroidDist] = shapeFun(cellMap, mean(cellMap(:))+3*std(cellMap(:)),params);
        if params.getcellEcc
            cellInfo.cellEcc(i) = cellEcc;
        end  
        if params.getcellSol
            cellInfo.cellSol(i) = cellSol;
        end         
        if params.getpScore
            cellInfo.pScore(i) = pScore;
        end         
        if params.getcScore
            cellInfo.cScore(i) = cScore;
        end 
        if params.getcellPer
            cellInfo.cellPer(i) = cellPer;
        end 
        
        if params.getcentroidDist
            cellInfo.centroidDist(i) = centroidDist;
        end
        
    end
    
    %EVENTS 
    if params.getimMovCorr || params.getinterCorr|| params.getinterSNR || params.geteventEcc || params.geteventSol || params.geteventP || params.geteventC || params.geteventS || params.geteventPer || params.geteventDist
        cropSize = 10;
		if ~isempty(testpeaksArray{i})
			[imMovCorr,interCorr,interSNR,eventEcc,eventSol,eventP,eventC,eventPer,eventDist,eventS] = corrFun(cellImages(:,:,i),inputMovie(:,:,testpeaksArray{i}),size(inputMovie),cropSize,params);
			if params.getimMovCorr
				cellInfo.imMovCorr(i) = imMovCorr;
			end
			if params.getinterCorr
				cellInfo.interCorr(i) = interCorr;
			end
			if params.getinterSNR
				cellInfo.interSNR(i) = interSNR;
			end
			if params.geteventEcc
				cellInfo.eventEcc(i) = eventEcc;
			end
			if params.geteventSol
				cellInfo.eventSol(i) = eventSol;
			end        
			if params.geteventP
				cellInfo.eventP(i) = eventP;
			end
			if params.geteventC
				cellInfo.eventC(i) = eventC;
			end
			if params.geteventPer
				cellInfo.eventPer(i) = eventPer;
			end
			if params.geteventDist
				cellInfo.eventDist(i) = eventDist;
            end
            if params.geteventS
                cellInfo.eventS(i) = eventS;
            end
		end			
    end
    
end
end

function [expRatio] = traceFun(inputSignal,testPeak,spikeROI,slopeFrameWindow,params)
    expRatio = NaN;

    extractMatrix = bsxfun(@plus,testPeak',spikeROI);
    extractMatrix(extractMatrix<=0)=1;
    extractMatrix(extractMatrix>=size(inputSignal,2))=size(inputSignal,2);
    spikeCenterTrace = reshape(inputSignal(extractMatrix),size(extractMatrix));
    
    if params.getexpRatio
        prePeakIdx = find(spikeROI==-(slopeFrameWindow)):find(spikeROI==0);
        postPeaklIdx = find(spikeROI==0):find(spikeROI==slopeFrameWindow);

        meanSignalPre =  mean(spikeCenterTrace(:,prePeakIdx),1);
        meanSignalPost =  mean(spikeCenterTrace(:,postPeaklIdx),1);

        f = {NaN,NaN};
        ind = 0;
        for meanSignal = [flip(meanSignalPre)',meanSignalPost']
            ind = ind + 1;
    %         [~,locs] = findpeaks(diff(meanSignal));
    %         if isempty(locs) %meanSignal has no minima
    %             [~,locs] = max(diff(meanSignal));
    %             locs = flip(locs);
    %         end
            trimMeanSignal = meanSignal;%meanSignal(1:locs(1));
            x = 1:size(trimMeanSignal);
            try
                f{ind} = fit(x',trimMeanSignal,'exp1');
            catch
                try
                    coeff = polyfit(x',log(trimMeanSignal),1);
                    f{ind} = struct('a' , 10^coeff(2), 'b', coeff(1));
                catch
                    warning('Could not fit the trace')
                    f{ind} = struct('a' , NaN, 'b', NaN);
                end
            end
        end

        fittedPre = f{1}.a*exp(f{1}.b*(1:0.01:size(meanSignalPre,2)));
        fittedPost = f{2}.a*exp(f{2}.b*(1:0.01:size(meanSignalPost,2)));
        
        expRatio = f{2}.b/f{1}.b;
    end
    if params.showPlots
        %plot([prePeakIdx,postPeaklIdx],[meanSignalPre,meanSignalPost]);
        myStd = [std(spikeCenterTrace(:,prePeakIdx),1),std(spikeCenterTrace(:,postPeaklIdx),1)];
        shadedErrorBar([prePeakIdx,postPeaklIdx],[meanSignalPre,meanSignalPost],myStd,{'b','linewidth',2})
        hold on;
        if params.getexpRatio
            h = plot(linspace(prePeakIdx(1),prePeakIdx(end),size(fittedPre,2)),flip(fittedPre),'r','linewidth',2);
            plot(linspace(postPeaklIdx(1),postPeaklIdx(end),size(fittedPost,2)),fittedPost,'r','linewidth',2);
            disp(['expRatio: ',num2str(expRatio)])
            legend(h,'Fitted')
        end

        plot([prePeakIdx(end),prePeakIdx(end)],[0,max(meanSignalPre)],'b--','linewidth',2)
        hold off;
        title('traceInfo')
        axis square
        grid on
        pause()
    end
end

function [cellEcc,cellSol,pScore,cScore,cellPer,centroidDist] = shapeFun(cellMap,threshold,params)
    [cellEcc,cellSol,pScore,cScore,cellPer,centroidDist] = deal(NaN);
    ellipse = [NaN;NaN];
    
    cellMapT = cellMap > threshold;

    if params.getpScore || params.getcScore || params.getcentroidDist
        
        [x,y] = find(cellMapT);
        
        if params.getcentroidDist
            centroid = [mean(x),mean(y)];
            cellMapCenter = size(cellMap)/2;
            centroidDist = sqrt((centroid(1)-cellMapCenter(1))^2 + (centroid(2)-cellMapCenter(2))^2)/sqrt((size(cellMap,1)-cellMapCenter(1))^2 + (size(cellMap,2)-cellMapCenter(2))^2);
            if params.showPlots
                disp(['centroidDist: ' , num2str(centroidDist)])
            end
        end
    
        [x,y] = find(bwperim(cellMapT));
        
        if params.getpScore
            try
                warning off
                [z, a, b, alpha] = fitellipse([x,y]');
                warning on
                perimeter = regionprops(cellMapT,'perimeter');
                realPerimeter = max([perimeter.Perimeter]);
                fittedPerimeter = ellipsePerimeter([z(1),z(2), a, b, alpha]);
                pScore = abs(1-(realPerimeter/fittedPerimeter));
                if params.showPlots
                    Q = [cos(alpha), -sin(alpha) 
                         sin(alpha), cos(alpha)];
                    ellipseTemp = Q * [a * cos(0:0.01:2*pi); b * sin(0:0.01:2*pi)];
                    ellipse = [z(1) + ellipseTemp(1,:);z(2) + ellipseTemp(2,:)];
                    disp([' '])
                    disp(['Perimeter score: ' , num2str(pScore)])
                end
            catch
                disp('Could not fit ellipse')
            end                
        end
        if params.getcScore
            euclDist = zeros([1,size(x)]);
            for k = 1:size(x)
                euclDist(k) = sqrt((x(k)-centroid(1))^2+(y(k)-centroid(2))^2);
            end
            cScore = std(euclDist); %This gives how similar the shape is to a circle.
            if params.showPlots
                disp([' '])
                disp(['Circle score: ',num2str(cScore)])
            end
        end
    end
    
    if params.getcellEcc || params.getcellSol || params.getcellPer
       regProps = regionprops('table',cellMapT,'perimeter','Eccentricity','Solidity');
       [~,index] = max([regProps.Perimeter]);
       if ~isempty(index);
        cellEcc = regProps.Eccentricity(index);
        cellSol = regProps.Solidity(index);
        cellPer = regProps.Perimeter(index);
       end
       if params.showPlots
        disp('')
       	disp(regProps)  
        disp('')
       end
    end
  
    if params.showPlots
        imagesc(cellMapT')
        title('Shapes')
        hold on
        %plot(ellipse(1,:),ellipse(2,:),'r','linewidth',2)
        axis square
        grid on
        hold off
        pause()
    end
end

function [globalSNR,RMS,peakSNR] = snrRmsFun(inputSignal,testPeak,timeSeq, spikeROI, params)
    [globalSNR,RMS,peakSNR] = deal(NaN);
    loopSignal = inputSignal;
    peaks = testPeak;
    peakIdx = bsxfun(@plus,timeSeq',peaks);
    peakIdx = unique(peakIdx(:));
    % remove peaks outside range of signal
    peakIdx(peakIdx>length(loopSignal))=[];
    peakIdx(peakIdx<=0)=[]; 
    % remove signal then add back in noise based on signal statistics
    noiseSignal = loopSignal;
    % noiseSignal(peakIdx) = NaN;
    noiseSignal(peakIdx) = [];
    
    if params.getpeakSNR
        
        extractMatrix = bsxfun(@plus,testPeak',spikeROI);
        extractMatrix(extractMatrix<=0)=1;
        extractMatrix(extractMatrix>=size(inputSignal,2))=size(inputSignal,2);
        spikeCenterTrace = reshape(inputSignal(extractMatrix),size(extractMatrix));
        
        peakSNR = mean(std(spikeCenterTrace,1))/nanstd(noiseSignal);%mean(mean(spikeCenterTrace,1)./std(spikeCenterTrace,1));
        if params.showPlots
            disp(['peakSNR: ',num2str(peakSNR)])  
        end
    end
    if params.getglobalSNR
        % compute SNR
        tmpSignal = inputSignal;
        % x_snr = nanmean(loopSignal)/nanstd(noiseSignal);
        globalSNR = nanmean(tmpSignal(peaks)/nanstd(noiseSignal));
        if params.showPlots
            disp(['globalSNR: ', num2str(globalSNR)])
        end
    end
    if params.getRMS
        % remove noise from signal vector
        xtmp = zeros([1 length(loopSignal)]);
        xtmp(peakIdx) = 1;
        loopSignal(~logical(xtmp)) = NaN;
        RMS = nanmean(abs(loopSignal)); %sqrt(nanmean(loopSignal.^2)) %Before it was computed as
        if params.showPlots
            disp(['RMS: ', num2str(RMS)])
        end
    end
end

function [imMovCorr,interCorr,interSNR,eventEcc,eventSol,eventP,eventC,eventPer,eventDist,eventS] = corrFun(inputImage,signalImages,movieDims,cropSize,params)
    [imMovCorr,interCorr,interSNR,eventEcc,eventSol,eventP,eventC,eventPer,eventDist,eventS] = deal(NaN);
    threshold = mean(inputImage(:))+3*std(inputImage(:));
    inputImageT = inputImage > threshold;

    [x,y] = find(inputImageT');   
    
    xHigh = max(x);
    yHigh = max(y);
    xLow = min(x);
    yLow = min(y);

    
    inputImageCrop = inputImage(yLow:yHigh,xLow:xHigh);
    signalImagesCrop = signalImages(yLow:yHigh,xLow:xHigh,:);
    eventMean = mean(signalImagesCrop,3);
    
    if params.geteventS
        eventS = sum(sum(abs(del2(eventMean))));
        if params.showPlots
            disp(['Event smooth score: ',num2str(eventS)])
        end
    end
    
    if params.geteventEcc || params.geteventSol || params.geteventP || params.geteventC || params.geteventPer ||  params.geteventDist
        
        if params.geteventEcc;eventParams.getcellEcc = true;else eventParams.getcellEcc = false;end
        if params.geteventSol;eventParams.getcellSol = true;else eventParams.getcellSol = false;end
        if params.geteventP;eventParams.getpScore = true;else eventParams.getpScore = false;end
        if params.geteventC;eventParams.getcScore = true;else eventParams.getcScore = false;end
        if params.geteventPer;eventParams.getcellPer = true;else eventParams.getcellPer = false;end
        if params.geteventDist;eventParams.getcentroidDist = true;else eventParams.getcentroidDist = false;end
        if params.showPlots;eventParams.showPlots = true;else eventParams.showPlots = false;end
        [eventEcc,eventSol,eventP,eventC,eventPer,eventDist] = shapeFun(eventMean,nanmedian(eventMean(:)),eventParams);
    end
    
    %       eventEcc
%       eventSol
%       eventP
%       eventC
%       eventPer

    if params.getinterCorr || params.getinterSNR
        corrMat = corrcoef(reshape(signalImagesCrop,[size(signalImagesCrop,1)*size(signalImagesCrop,2),size(signalImagesCrop,3)]));
        corrMatUpper = triu(corrMat,1);
        interCorr = mean(corrMatUpper(corrMatUpper~=0));
        interSNR = 1/(1/mean(corrMatUpper(corrMatUpper~=0))-1);
    end
    if params.getimMovCorr
        cellMask = find(inputImageCrop == 0);

        inputImageCrop = -1 + 2.*(inputImageCrop - min(min(inputImageCrop)))./(max(max((inputImageCrop))) - min(min(inputImageCrop)));

        inputImageCrop(cellMask) = NaN;
        for k = 1:size(signalImagesCrop,3)
            tmp = signalImagesCrop(:,:,k);
            tmp2 = tmp(:);
            tmp2(cellMask) = NaN;
            tmp2 = -1 + 2.*(tmp2 - min(tmp2))./(max(tmp2) - min(tmp2)); %Normalize betweeen 0 and 1.
            signalImagesCrop(:,:,k) = reshape(tmp2,size(tmp));
        end
    end
    peakImageN = size(signalImagesCrop,3);
    corrVals = zeros([peakImageN,1]);
  
    if params.showPlots
        title('Correlations between events and filter')
        myAxes = zeros(1, ceil(sqrt(peakImageN+2)));
        myAxes(1) = subplot(ceil(sqrt(peakImageN+2)),ceil(sqrt(peakImageN+2)),1);
        imagesc(inputImageCrop)
    end

    for k = 1:peakImageN
        if params.getimMovCorr
            corrVals(k) = myCorr2(inputImageCrop,signalImagesCrop(:,:,k));
        end
        if params.showPlots
            myAxes(k+2) = subplot(ceil(sqrt(peakImageN+2)),ceil(sqrt(peakImageN+2)),k+2);
            %imagesc(signalImagesCrop(:,:,k))
            imagesc(signalImages(yLow:yHigh,xLow:xHigh,k))
            colormap gray
        end
    end
    if params.getimMovCorr
        imMovCorr = nanmean(corrVals);
    end
    if params.showPlots
        disp(' ')
        if params.getimMovCorr
            disp(['imMovCorr ',num2str(imMovCorr)])
        end
        if params.getinterCorr
            disp(['InterCorrelation: ',num2str(interCorr)])
        end
        if params.getinterSNR
            disp(['SNRCorrelation: ',num2str(interSNR)])
        end
        
        pause()
        delete(myAxes)
    end
    
    
%     threshold = mean(inputImage(:))+3*std(inputImage(:));
%     inputImageT = inputImage > threshold;
% 
%     [x,y] = find(inputImageT');
% 
%     centroidCoords = round(mean([x,y]));
% 
%     xLow = centroidCoords(1) - cropSize;
%     xHigh = centroidCoords(1) + cropSize;
%     yLow = centroidCoords(2) - cropSize;
%     yHigh = centroidCoords(2) + cropSize;
%     % check that not outside movie dimensions
%     xMin = 1 ;
%     xMax = movieDims(2);
%     yMin = 1 ;
%     yMax = movieDims(1);
% 
%     % adjust for the difference in centroid location if movie is cropped
%     xDiff = 0;
%     yDiff = 0;
%     if xLow<xMin xDiff = xLow-xMin; xLow = xMin; end
%     if xHigh>xMax xDiff = xHigh-xMax; xHigh = xMax; end
%     if yLow<yMin yDiff = yLow-yMin; yLow = yMin; end
%     if yHigh>yMax yDiff = yHigh-yMax; yHigh = yMax; end
% 
%     % [yLow yHigh xLow xHigh]
%     % crop to a region of x-y pixels around the cell's centroid
%     signalImagesCrop = signalImages(yLow:yHigh,xLow:xHigh,:);
% 
%     signalImagesCropTmp = signalImagesCrop;
%     signalImagesCropTmp(isnan(signalImagesCropTmp)) = 0;
%     inputImageCrop = inputImage(yLow:yHigh,xLow:xHigh);
%     peakImageN = size(signalImagesCropTmp,3);
%     corrVals = zeros([peakImageN,1]);
%     
%     if params.showPlots
%         title('Correlations between events and filter')
%         myAxes = zeros(1, ceil(sqrt(peakImageN+2)));
%         myAxes(1) = subplot(ceil(sqrt(peakImageN+2)),ceil(sqrt(peakImageN+2)),1);
%         imagesc(inputImageCrop)
%     end
% 
%     for k = 1:peakImageN
%         corrVals(k) = corr2(inputImageCrop,signalImagesCropTmp(:,:,k));
%         if params.showPlots
%             myAxes(k+2) = subplot(ceil(sqrt(peakImageN+2)),ceil(sqrt(peakImageN+2)),k+2);
%             imagesc(signalImagesCropTmp(:,:,k))
%         end
%     end
%     imMovCorr = nanmean(corrVals);
%     if params.showPlots
%         disp(' ')
%         disp(['imMovCorr ',num2str(imMovCorr)])
%         pause()
%         delete(myAxes)
%     end
end

function c2 = myCorr2(A,B)
    AMean = nanmean(A(:));
    BMean = nanmean(B(:));
    
    c2 = nansum((A(:)-AMean).*(B(:)-BMean))./sqrt(nansum((A(:)-AMean).^2).*nansum((B(:)-BMean).^2)); 
end
