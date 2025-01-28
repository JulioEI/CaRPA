function [eventTimes, cellTraceSigmas, options] = detectEventsOnPhi(cellTraces, scaledPhi, cellImages, centroids, DFOF, varargin)
    % [eventTimes, options] = detectEventsOnPhi(cellTraces, scaledPhi, cellImages, centroids, imgs, varargin)
    %
    % Detects events on the scaledProbability output of EM using either a
    % correlation based method or a simple thresholding approach.
    %
    % Inputs:
    % cellTraces: EM output dsCellTraces
    % scaledPhi: EM output dsScaledProbability
    % cellImages: EM output cellImages
    % centroids: EM output centroids
    % imgs: temporally downsampled DFOF
    %
    % Options:
    %
    %	options.framerate: Framerate of imgs in Hz.
    %       Default = 4
    %	options.burstDuration: Length of time (in seconds) that burst is
    %       required to remain above threshold.
    %       Default = 1
    %	options.refractoryPeriod: Minimum time (in seconds) between bursts.
    %   	Default = 1.5
    %   options.useImageCorrelation: Determines which method to use to identify
    %       valid events. If useImageCorrelation = 1, will use the correlation
    %       between the cell's image and the movie frame at the detected peak,
    %       works well in eliminating crosstalk from neighboring cells. When
    %       useImageCorrelation = 0, will usea  thresholding approach.
    %       Default = 1
    %   options.useImageCorrelationTrace: Whether to use scaled probability trace for peak detection, make compatible with other trace types
    %
    % Options if useImageCorrelation = 1
    %   options.nBaselineFrames: The number of frames that will be used to
    %       estimate the baseline correlation
    %       Default = 500;
    %   options.correlationWindow: The number of pixels out from the centroid
    %       of the image that will be used for the correlation.
    %       Default = 10
    %   options.prctileCorr: Prctile of baseline measure that EM-Event correlation
    %       must exceed for event to be considered valid.
    %       Default = 95
    %   options.prctileCorrThresh: Prctile of baseline measure that EM-Event
    %       thresholded correlation must exceed for event to be considered valid.
    %       Default = 95
    %   options.prctileDistance: Prctile of baseline measure that EM-Event
    %       distance must exceed for event to be considered valid.
    %       Default = 5
    %   options.minThresh: Threshold used to isolate events within the event image.
    %       Default = 0.75
    %
    % Options if useImageCorrelation = 0
    %   options.numSigmasThresh: The number of std devs above baseline used to
    %       identify valid events if useImageCorrelation=0.
    %       Default = 5
    %   options.diffThreshold: The distance above their neighbors for a peak to
    %       be considered valid.
    %       Default = 0
    %
    % Written by Lacey Kitch and Maggie Carr Larkin in 2013

    % 2017.02.18 [17:16:34] updated to make more compatible with non-Phi traces - biafra
    %--------------------------------------------------------------------------

    %Set default options
    options = getDefaultEventDetectionOptions();

    % Get user inputs
    options=getOptions(options,varargin);
    warning('off', 'signal:findpeaks:largeMinPeakHeight')

    % Initialize gassian for smoothing eventImages
    g = options.smoothFilter;

    %Preallocate outputs
    nCells = size(cellTraces,1);
    eventTimes=cell(nCells,1);
    if options.readMovieChunks==0
        size(DFOF)
        movieSize=size(DFOF);
    else
        % movieSize = DFOF;
        movieDims = loadMovieList(options.movieFilename,'getMovieDims',1,'inputDatasetName',options.movieDatasetName);
        yDim = movieDims.one;
        xDim = movieDims.two;
        zDim = movieDims.three;
        movieSize = [yDim xDim zDim];
    end

    % Perform event detection
    if options.useImageCorrelation
        cellTraceSigmas = [];
        reverseStr = '';
        for cInd=1:nCells
            if mod(cInd,25)==0
                reverseStr = cmdWaitbar(cInd,nCells,reverseStr,'inputStr','calculating traces','waitbarOn',1,'displayEvery',25);
            end

            clear burst_thresh
            if options.useImageCorrelationTrace==0
                thisPhi = scaledPhi(cInd,:);
                burst_thresh = min(0.75*max(thisPhi),0.5);
            else
                thisPhi = cellTraces(cInd,:);
                %Find all potential spiketimes
                inputSignalStd = std(thisPhi(:));
                burst_thresh = inputSignalStd*options.numSigmasThresh;
            end

            [~, spiketimes] = findpeaksFast(thisPhi,'minpeakheight',burst_thresh,...
                'minpeakdistance',round(options.refractoryPeriod*options.framerate));
            if options.displayInfo==1
                figure;plot(thisPhi);hold on;plot(spiketimes,thisPhi(spiketimes),'k.');
            end

            if ~isempty(spiketimes)
                %Identify this EM, thresholded EM, and the EM weighted centroid
                xLims=max(1,round(centroids(cInd,1)-options.correlationWindow)):min(movieSize(2), round(centroids(cInd,1)+options.correlationWindow));
                % xLims
                yLims=max(1,round(centroids(cInd,2)-options.correlationWindow)):min(movieSize(1), round(centroids(cInd,2)+options.correlationWindow));
                % yLims

                referenceImage = cellImages(yLims,xLims,cInd);
                threshReferenceImage = getBinImage(referenceImage,g);
                [threshReferenceImage,nRegions] = bwlabel(threshReferenceImage);
                if nRegions>1
                    Area = hist(threshReferenceImage(:),0:1:nRegions);
                    [~,loc] = max(Area(2:end));
                    threshReferenceImage = threshReferenceImage==loc; clear loc Area
                end
                referenceCentroid = calculateWeightedCentroids(threshReferenceImage,referenceImage);

                %Determine baseline values for this EM
                if options.useImageCorrelationTrace==0
                    baseline_times = find(thisPhi<=options.phiBaselineThresh);
                else
                    baseline_times = find(thisPhi<=burst_thresh);
                end
                baseline_ind = baseline_times(randperm(length(baseline_times),min(length(baseline_times),options.nBaselineFrames))); clear baseline_times
                baseline_corr = zeros(size(baseline_ind));
                baseline_corr_thresh = zeros(size(baseline_ind));
                baseline_centroidDistance = zeros(size(baseline_ind));

                if options.displayInfo==1
                    fprintf('baseline_ind length = %d\n',length(baseline_ind))
                end
                threshCentroidAll = [];

                for j=1:length(baseline_ind)
                    try
                        if options.readMovieChunks==0
                            thisImage=DFOF(yLims,xLims,baseline_ind(j));
                        else
                            thisImage = readHDF5Subset(options.movieFilename,[min(yLims) min(xLims) baseline_ind(j)-1],[length(yLims) length(xLims) 1],'datasetName',options.movieDatasetName,'displayInfo',0);
                            thisImage(isnan(thisImage)) = 0;
                        end
                    catch
                        continue;
                    end
                    thisImage = double(thisImage);
                    %Calcualte correlation of the image with referenceImage
                    tmp = corrcoef(referenceImage(:),thisImage(:));
                    baseline_corr(j) = tmp(2); clear tmp

                    %Determine the value
                    %Determine the thresholded image
                    % [thisImageThresh,threshCentroid] = calculateThreshImage(thisImage,referenceCentroid,g,options.minThresh); clear thisImage
                    [thisImageThresh,threshCentroid] = calculateThreshImage(thisImage,[round(rand*size(thisImage,1)) round(rand*size(thisImage,2))],g,options.minThresh); clear thisImage

                    %Calculate correlation between thresholded images
                    tmp = corrcoef(threshReferenceImage(:),double(thisImageThresh(:)));
                    baseline_corr_thresh(j) = tmp(2); clear tmp thisImageThresh

                    %Measure the distance between this images centroid and the referenceImage
                    baseline_centroidDistance(j) = sqrt(sum((referenceCentroid-threshCentroid).^2));

                    threshCentroidAll = [threshCentroidAll;threshCentroid];
                end
                threshCentroidAll = mean(threshCentroidAll,1);
                baseline_corr = prctile(baseline_corr,options.prctileCorr);
                baseline_corr_thresh = prctile(baseline_corr_thresh,options.prctileCorr);
                if options.displayInfo==1
                    fprintf('nanmean(baseline_centroidDistance) = %f\n prctile(baseline_centroidDistance,options.prctileDistance) = %f',nanmean(baseline_centroidDistance),prctile(baseline_centroidDistance,options.prctileDistance))
                end
                baseline_centroidDistance = prctile(baseline_centroidDistance,options.prctileDistance);
                % baseline_centroidDistance = nanmean(baseline_centroidDistance);
                clear baseline_ind

                %Determine the correlation between each detected burst and the cellImage
                event_corr = zeros(size(spiketimes));
                event_corr_thresh = zeros(size(spiketimes));
                event_centroidDistance = zeros(size(spiketimes));
                for j = 1:length(spiketimes)
                    try
                        if options.readMovieChunks==0
                            thisImage=DFOF(yLims,xLims,spiketimes(j));
                        else
                            thisImage = readHDF5Subset(options.movieFilename,[min(yLims) min(xLims) spiketimes(j)-1],[length(yLims) length(xLims) 1],'datasetName',options.movieDatasetName,'displayInfo',0);
                            thisImage(isnan(thisImage)) = 0;
                        end
                    catch
                        continue;
                    end
                    thisImage = double(thisImage);
                    %Calcualte correlation of the image with referenceImage
                     tmp = corrcoef(referenceImage(:),thisImage(:));
                     event_corr(j) = tmp(2); clear tmp

                    %Determine the thresholded image
                    % [thisImageThresh,threshCentroid] = calculateThreshImage(thisImage,referenceCentroid,g,options.minThresh);
                    [thisImageThresh,threshCentroid] = calculateThreshImage(thisImage,[round(rand*size(thisImage,1)) round(rand*size(thisImage,2))],g,options.minThresh);

                    %Calculate correlation between thresholded images
                    tmp = corrcoef(threshReferenceImage(:),double(thisImageThresh(:)));
                    event_corr_thresh(j) = tmp(2);

                    %Measure the distance between this images centroid and the referenceImage
                    event_centroidDistance(j) = sqrt(sum((referenceCentroid-threshCentroid).^2));

                    if options.displayInfo==1
                        event_corr1 = event_corr(j)>baseline_corr;
                        event_corr_thresh1 = event_corr_thresh(j)>baseline_corr_thresh;
                        event_centroidDistance1 = event_centroidDistance(j)<=baseline_centroidDistance;
                        caxisMin = 0;
                        caxisMax = 0.1;
                        xplot = 6;
                        yplot = 2;
                        subplot(xplot,yplot,1);
                            imagesc(threshReferenceImage);hold on
                            plot(referenceCentroid(1),referenceCentroid(2),'r+','MarkerSize',12)
                            plot(threshCentroid(1),threshCentroid(2),'g.','MarkerSize',12);
                            plot(threshCentroidAll(1),threshCentroidAll(2),'bx','MarkerSize',12);hold off;
                            legend({'ref','event','baseline'},'Location','eastoutside')
                            caxis([caxisMin caxisMax])
                            colormap gray
                            title('reference cell thresholded')
                        subplot(xplot,yplot,2);
                            imagesc(thisImageThresh);hold on;
                            plot(referenceCentroid(1),referenceCentroid(2),'r+','MarkerSize',12)
                            plot(threshCentroid(1),threshCentroid(2),'g.','MarkerSize',12);
                            plot(threshCentroidAll(1),threshCentroidAll(2),'bx','MarkerSize',12);hold off;
                            legend({'ref','event','baseline'},'Location','eastoutside')
                            caxis([caxisMin caxisMax])
                            title('event movie thresholded')
                        subplot(xplot,yplot,[3 4]);
                            plot(event_corr_thresh);
                            hline = refline([0 baseline_corr_thresh]);
                            hline.Color = 'r';
                            hold on;plot(j,event_corr_thresh(j),'k.');hold off;
                        subplot(xplot,yplot,5);
                            imagesc(referenceImage);
                            caxis([caxisMin max(referenceImage(:))])
                        subplot(xplot,yplot,6);
                            imagesc(thisImage);
                            caxis([caxisMin caxisMax])
                        subplot(xplot,yplot,[7 8]);
                            plot(event_corr);
                            hline = refline([0 baseline_corr]);
                            hline.Color = 'r';
                            hold on;plot(j,event_corr(j),'k.');hold off;
                        subplot(xplot,yplot,[9 10]);
                            plot(event_centroidDistance);
                            hline = refline([0 baseline_centroidDistance]);
                            hline.Color = 'r';
                            hold on;plot(j,event_centroidDistance(j),'k.');hold off;
                        subplot(xplot,yplot,[11 12])
                            plot(event_corr>baseline_corr&event_corr_thresh>baseline_corr_thresh&event_centroidDistance<=baseline_centroidDistance)

                        backgroundGood = [208,229,180]/255;
                        backgroundBad = [244,166,166]/255;
                        backgroundNeutral = repmat(230,[1 3])/255;
                        decisionVal = event_corr1&event_corr_thresh1&event_centroidDistance1;
                        if decisionVal==1
                            set(gcf,'Color',backgroundGood);
                        elseif decisionVal==0
                            set(gcf,'Color',backgroundBad);
                        else
                            set(gcf,'Color',backgroundNeutral);
                        end

                        titleStr = fprintf('%f | %f | %f | %f | %f | %f\n',event_corr(j),baseline_corr,event_corr_thresh(j),baseline_corr_thresh,event_centroidDistance(j),baseline_centroidDistance);
                        title([num2str(event_corr1&event_corr_thresh1&event_centroidDistance1) ' | ' num2str(event_corr1&event_corr_thresh1)])
                        drawnow
                        pause
                    end
                end

                %Valid events are more correlated than baseline and their centroids
                %are closer to the EM centroid than baseline values.
                % event_corr
                % baseline_corr
                % baseline_corr_thresh
                % event_corr_thresh
                % baseline_corr_thresh
                event_corr = event_corr>baseline_corr;
                event_corr_thresh = event_corr_thresh>baseline_corr_thresh;
                event_centroidDistance = event_centroidDistance<=baseline_centroidDistance;

                if options.displayInfo==1
                    figure;
                    % plot(event_corr & event_corr_thresh & event_centroidDistance);
                    subplot(3,1,1);plot(event_corr);subplot(3,1,2);plot(event_corr_thresh);subplot(3,1,3);plot(event_centroidDistance);drawnow;
                end
                eventTimes{cInd} = spiketimes(event_corr & event_corr_thresh & event_centroidDistance);
                clear xLims yLims referenceImage threshReferenceImage
            end
        end
    else
        cellTraceSigmas = zeros(nCells,1);
        for cInd=1:nCells

            thisPhi = scaledPhi(cInd,:);
            thisTrace = cellTraces(cInd,:);

            %Ignore really small values
            threshTrace = thisTrace(thisPhi<0.1);
            if ~isempty(threshTrace)
                cellTraceSigmas(cInd)=std(threshTrace);
                probabilityThreshold = options.numSigmasThresh*cellTraceSigmas(cInd);

                %Find peaks that are greater than the probability threshold and far enough apart
                [~, spiketimes] = findpeaks(thisPhi,'minpeakheight',probabilityThreshold,...
                    'minpeakdistance',round(options.refractoryPeriod*options.framerate),...
                    'threshold', options.diffThreshold);

                %Require that the trace remains above threshold for long enough
                eventTimes{cInd} = intersect(spiketimes, find(filtfilt(ones(1,...
                    round(options.burstDuration*options.framerate))/...
                    round(options.burstDuration*options.framerate),1,...
                    thisPhi)>(probabilityThreshold)));
            end
        end
    end
end

%--------------------------------------------------------------------------
function [pks,locs] = findpeaksFast(X,varargin)
    % FINDPEAKSFAST Find local peaks in data

    %Set default options
    minPeakDistance = 1;
    threshold = 0;
    minPeakHeight = Inf;

    for option = 1:2:length(varargin)-1
        if ischar(varargin{option})
            switch(lower(varargin{option}))
                case 'minpeakheight'
                    minPeakHeight = varargin{option+1};
                case 'minpeakdistance'
                    minPeakDistance = varargin{option+1};
                case 'threshold'
                    threshold = varargin{option+1};
                otherwise
                    error(['Option ',varargin{option},' unknown.']);
            end
        else
            error('Options must be strings, followed by the variable');
        end
    end

    pks = []; locs = [];
    if all(isnan(X)),
        return,
    end

    % Replace Inf by realmax because the diff of two Infs is not a number
    infIdx = isinf(X);
    if any(infIdx),
        X(infIdx) = sign(X(infIdx))*realmax;
    end
    infIdx = infIdx & X>0; % Keep only track of +Inf

    % Get Peaks Above Min Peak Height
    Indx = find(X > minPeakHeight);
    if(isempty(Indx))
        warning(message('signal:findpeaks:largeMinPeakHeight', 'MinPeakHeight', 'MinPeakHeight'));
        return
    end

    trend = sign(diff(X));
    idx = find(trend==0); % Find flats
    N = length(trend);
    for i=length(idx):-1:1,
        % Back-propagate trend for flats
        if trend(min(idx(i)+1,N))>=0,
            trend(idx(i)) = 1;
        else
            trend(idx(i)) = -1; % Flat peak
        end
    end

    idx  = find(diff(trend)==-2)+1;  % Get all the peaks
    locs = ismember(Indx,idx);      % Keep peaks above MinPeakHeight
    locs = Indx(locs);
    pks  = X(locs);

    %Remove Peaks Below Threshold
    idelete = [];
    for i = 1:length(pks),
        delta = min(pks(i)-X(locs(i)-1),pks(i)-X(locs(i)+1));
        if delta<threshold,
            idelete = [idelete i]; %#ok<AGROW>
        end
    end
    if ~isempty(idelete),
        locs(idelete) = [];
    end

    if any(infIdx)
        X(infIdx) = Inf;                 % Restore +Inf
        locs = union(locs,find(infIdx)); % Make sure we find peaks like [realmax Inf realmax]
    end
    pks  = X(locs);

    % Remove Peaks Separated By Less Than Min Peak Distance
    if ~isempty(pks) && minPeakDistance~=1

        % Order peaks from large to small
        [pks, idx] = sort(pks,'descend');
        locs = locs(idx);

        idelete = ones(size(locs))<0;
        for i = 1:length(locs),
            if idelete(i) == 0
                % If the peak is not in the neighborhood of a larger peak, find
                % secondary peaks to eliminate.
                idelete = idelete | (locs>=locs(i)-minPeakDistance)&(locs<=locs(i)+minPeakDistance);
                idelete(i) = 0; % Keep current peak
            end
        end
        pks(idelete) = [];
        locs(idelete) = [];
    end

    %Order Peaks

    if isempty(pks), return; end

    [locs, idx] = sort(locs);
    pks = pks(idx);
end

%--------------------------------------------------------------------------
function WeightedCentroid = calculateWeightedCentroids(bw,I)

    WeightedCentroid = [0 0];

    sum_region = sum(I(bw~=0));

    [y,x] = find(bw); %Indices

    WeightedCentroid(1) = sum(x .* double(I(bw~=0))) / sum_region; %X
    WeightedCentroid(2) = sum(y .* double(I(bw~=0))) / sum_region; %Y
end

%--------------------------------------------------------------------------
function binImage = getBinImage(cellImage,g)
    if ~isempty(g)
        binImage = conv2(cellImage,g,'same');
    else
        binImage = cellImage;
    end

    maxVal=max(binImage(:));
    binImage = binImage>=0.5*maxVal;
end

%--------------------------------------------------------------------------
function [labeledImage,weighted_centroid] = calculateThreshImage(image,peak,smoothfilt,minThresh)
    %Isolate thresholded event image

    %Smooth Image
    if ~isempty(smoothfilt)
        image = conv2(image,smoothfilt,'same');
    end

    %Normalize
    imageNorm = (image-min(image(:)))./(max(image(:))-min(image(:)));

    %Identify closest peak
    [imageThresh,nRegions] = bwlabel(imageNorm > minThresh); clear imageNorm
    closestPeak = nan(nRegions,2);
    for val = 1:nRegions
        [~,maxInd]=max(image(imageThresh==val));
        [y,x] = find(imageThresh==val);
        closestPeak(val,:) = [x(maxInd) y(maxInd)]; clear maxInd x y
    end
    clear val

    if nRegions>1
        distance_to_center = bsxfun(@minus,closestPeak,peak).^2;
        distance_to_center = sqrt(distance_to_center(:,1)+distance_to_center(:,2));
        [~,minLoc] = min(distance_to_center);
        closestPeak = closestPeak(minLoc,:);
        clear minLoc distance_to_center imageNorm

        %Isolate each burst to create a thresholded event image
        imageThresh = (image-min(image(:)))./(image(closestPeak(2),closestPeak(1))-min(image(:)));

        %Identify whether there are two peaks in image
        maxVal = 1;
        [thresh,threshLoc] = max(imageThresh(:));
        [threshLocx,threshLocy]=ind2sub(size(imageThresh),threshLoc);
        if sqrt((closestPeak(2)-threshLocx).^2 + (closestPeak(1)-threshLocy).^2) < 5
            maxVal = thresh;
            closestPeak = [threshLocy threshLocx];
        end

        if thresh>maxVal
            %Need to identify the minThresh possible to ignore this second peak
            labeledImage = bwlabel(imageThresh>=thresh);
            thresh = maxVal;
            while labeledImage(closestPeak(2),closestPeak(1))~=labeledImage(threshLocx,threshLocy) && thresh>minThresh
                thresh = thresh - 0.05;
                labeledImage = bwlabel(imageThresh>=thresh);
            end
            minThresh = max(minThresh,thresh + 0.1);
            thresh = maxVal;
        end
        labeledImage = bwlabel(imageThresh>=thresh & imageThresh<=maxVal);

        while thresh>minThresh
            %Iterate lowering the threshold and removing nonrelevant peaks
            thresh = thresh - 0.05;
            labeledImage = bwlabel(imageThresh>=thresh & imageThresh<=maxVal);
            valid = labeledImage(closestPeak(2),closestPeak(1));
            imageThresh(labeledImage~=valid & labeledImage>0) = maxVal+1;
        end
        labeledImage = labeledImage== labeledImage(closestPeak(2),closestPeak(1));
    else
        labeledImage = imageThresh;
    end

    %Measure weighted centroid
    weighted_centroid = calculateWeightedCentroids(labeledImage,image);
end