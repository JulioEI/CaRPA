function [eventTimes, cellTraceSigmas, options] = detectEventsOnPhiOld(cellTraces, scaledPhi, cellImages, centroids, imgs, varargin)

% Written by Lacey Kitch and Maggie Carr Larkin in 2013

% cellTraces is nCells x nFrames

%%%%% Input options:

% framerate of input movie
options.framerate=4;

% option to use the correlation between the cell's image and the movie
% frame at the detected peak - works well in eliminating crosstalk from
% neighboring cells
options.useImageCorrelation=0;

% length of time (in seconds) that burst is required to remain above
% threshold (on average)
options.burstDuration=1;    % duration in seconds

% minimum time (in seconds) between bursts
options.refractoryPeriod=1.5;

% the number of std devs above baseline we require (if not using image
% correlation)
options.numSigmasThresh=3;

% option to require peaks to be a minimum distance above their neighboring
% peaks (in time, in the same trace)
options.diffThreshold=0;

% If options.useImageCorrelation=1, this is the number of pixels out
% from the centroid of the image that will be used for the correlation
options.correlationWindow=10;

% If options.useImageCorrelation=1, this is the number of frames that
% will be used to estimate the baseline correlation
options.nBaselineFrames=1000;

% If options.useImageCorrelation=1, this is the number of correlation
% standard deviations that the image must be above (where the distribution
% of correlations is calculated on the baseline times)
options.numSigmasCorrThresh=4;

% get the input options
options=getOptions(options,varargin);
warning('off', 'signal:findpeaks:largeMinPeakHeight')


% perform event detection
nCells = size(cellTraces,1);
movieSize=size(imgs);
cellTraceSigmas = zeros(nCells,1);
eventTimes=cell(nCells,1);
for cInd=1:nCells

    if options.useImageCorrelation

        thisPhi = scaledPhi(cInd,:);

        burst_thresh = min(0.2*max(thisPhi),0.2);
        [~, spiketimes] = findpeaks(thisPhi,'minpeakheight',burst_thresh,...
            'minpeakdistance',round(options.refractoryPeriod*options.framerate));

        %Goal: identify the smallest amplitude spiketime

        %Work from the minimum possible spiketime, determine if the event image is
        %correlated enough, build evidence

        %Determine the baseline correlation
        xLims=max(1,round(centroids(cInd,1)-options.correlationWindow)):min(movieSize(2), round(centroids(cInd,1)+options.correlationWindow));
        yLims=max(1,round(centroids(cInd,2)-options.correlationWindow)):min(movieSize(1), round(centroids(cInd,2)+options.correlationWindow));

        referenceImage = cellImages(yLims,xLims,cInd);

        baseline_times = find(thisPhi<=0.1);
        baseline_ind = baseline_times(randperm(length(baseline_times),min(length(baseline_times),options.nBaselineFrames)));
        eventImages=imgs(yLims,xLims,baseline_ind);
        baseline_corr = zeros(size(baseline_ind));
        for j=1:length(baseline_ind)
            thisImage = eventImages(:,:,j);
            baseline_corr(j) = corr(referenceImage(:),thisImage(:),'rows','pairwise');
        end
        baseline_corr = std(baseline_corr);

        %Determine the correlation between each detected burst and the cellImage
        eventImages = imgs(yLims,xLims,spiketimes);
        event_corr = zeros(size(spiketimes));
        for j = 1:length(spiketimes)
            thisImage = eventImages(:,:,j);
            event_corr(j) = corr(referenceImage(:),thisImage(:),'rows','pairwise');
        end

        %Valid events are more correlated than baseline
        event_corr = event_corr./baseline_corr;

        eventTimes{cInd} = spiketimes(event_corr>options.numSigmasCorrThresh);

    else
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