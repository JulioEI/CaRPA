function [eventTimes, peakTimes, eventAmps, eventBin, cellTraceSigmas, options] = detectEvents(cellTraces, varargin)

    % Written by Lacey Kitch in 2013
    % Based on code written by Laurie Burns in 2011

    % input an options structure that overwrites default fields with
    %    detectEvents(cellTraces, 'options', options)

    % cellTraces is nCells x nFrames
    % options.numSigmasThresh is the number of std devs above baseline we require
    % options.calcSigma is a toggle, if 1, we fit the std of each trace
    % options.noiseSigma is the noise std dev of the movie, we use this if options.calcSigma=0

    options.numSigmasThresh=3;
    options.calcSigma=1;
    options.noiseSigma=0.005;
    options.doSubMedian=0;
    options.medFiltSize=200;
    options.doMovAvg=0;
    options.framerate=5;    % in Hz
    options.movAvgFiltSize=0.25; % in seconds
    options.movAvgReqSize=0.5; % in seconds
    options.reportMidpoint=0;
    options.minTimeBtEvents=0.25; % in seconds
    options.offsetFrames=0; % offset to account for ca2+ rise
    options.displayPeaks=0;
    options.window=800;
    options=getOptions(options,varargin);

    options.window=min(options.window, size(cellTraces,2)-1);
    options.movAvgFiltSize=ceil(options.movAvgFiltSize*options.framerate);
    options.movAvgReqSize=ceil(options.movAvgReqSize*options.framerate);
    options.minTimeBtEvents=ceil(options.minTimeBtEvents*options.framerate);
    nCells = size(cellTraces,1);
    cellTraceSigmas = zeros(nCells,1);


    % find trace std dev, or fit it
    if options.calcSigma

        %fInc=0.0005;
        %fVals=(-0.05:fInc:0.05)';

        % commented code does this:
        % Get SD by fitting gaussian to the histogrammed data (because of heavy pos
        % tail from bursting)
        % current code takes the std dev in sliding windows along the trace and
        % finds the minimum as an estimate of noise std dev
        reverseStr = '';
        for cellnum = 1 : nCells
            if mod(cellnum,25)==0
                reverseStr = cmdWaitbar(cellnum,nCells,reverseStr,'inputStr','detecting events','waitbarOn',1,'displayEvery',25);
            end
            %xdata = linspace(-0.3,0.3,1000)';         %%%% might want to take a look at these values
    %         if options.doSubMedian
    %             ydata = hist(cellTraces(cellnum,:) - ...
    %                 ordfilt2(cellTraces(cellnum,:)',ordFiltOrder,ones(ordFiltDomain,1),'symmetric')',fVals)';
    %         else
    %             ydata=hist(cellTraces(cellnum,:)-mean(cellTraces(cellnum,:)), fVals)';
    %         end

            if options.doSubMedian
                cellTraces(cellnum,:) = cellTraces(cellnum,:) - ordfilt2(cellTraces(cellnum,:)',round(options.medFiltSize/2),ones(options.medFiltSize,1),'symmetric');
            end

    %         options = fitoptions('method','nonlinearleastsquares',...
    %             'startpoint',[100 0 0.005]);
    %         try
    %             fitres = fit(fVals(2:end-1),ydata(2:end-1),'gauss1',options);
    %             cellTraceSigmas(cellnum,1) = (fitres.c1/sqrt(2));
    %         catch %#ok<CTCH>
    %             disp('fit failed...')
    %             cellTraceSigmas(cellnum,1) = std(cellTraces(cellnum,:));
    %         end

            trace=cellTraces(cellnum,:);
            windowStarts=1:20:(length(trace)-options.window);
            stdDevs=zeros(length(windowStarts),1);
            for wInd=1:length(windowStarts)
                windowStart=windowStarts(wInd);
                stdDevs(wInd)=std(trace(windowStart:min(length(trace),windowStart+options.window-1)));
            end
            cellTraceSigmas(cellnum,1) = min(stdDevs);
        end
    else
        cellTraceSigmas=options.noiseSigma*ones(nCells,1);
    end


    % Run the peak finding
    eventTimes = cell(nCells,1);
    peakTimes = cell(nCells,1);
    offsetpeaks = cell(nCells,1);
    riseTimes = cell(nCells,1);
    thresh = zeros(nCells,1);


    % get the traces shifted by the filtered trace
    if options.doSubMedian
        paddedCellTraces=padarray(cellTraces, [0 ceil(options.medFiltSize/2)]);
        medSubTraces=paddedCellTraces-medfilt1(paddedCellTraces',options.medFiltSize)';
        medSubTraces(:,1:ceil(options.medFiltSize/2))=[];
        medSubTraces(:,end-ceil(options.medFiltSize/2)+1:end)=[];
    else
        medSubTraces=cellTraces;
    end

    reverseStr = '';
    nCellTraces = size(cellTraces,1);
    for c=1:nCellTraces

        if mod(c,25)==0
            reverseStr = cmdWaitbar(c,nCellTraces,reverseStr,'inputStr','detecting events','waitbarOn',1,'displayEvery',25);
        end

        offsetpeaks{c} = zeros(1,0);
        % set the threshold to a multiple of the SD of the trace for that cell
        thresh(c) = options.numSigmasThresh*cellTraceSigmas(c);

        if options.doMovAvg
            inputsignal = filtfilt(ones(1,options.movAvgFiltSize)/options.movAvgFiltSize,1,medSubTraces(c,:));
        else
            inputsignal=medSubTraces(c,:);
        end

        % find peaks that satisfy the minimum height, the minimum number of
        % frames between peaks, and the required moving average size above
        % thresh
        if max(inputsignal)>=thresh(c)
            [~,testpeaks] = findpeaks(inputsignal,'minpeakheight',thresh(c));
            [~,testpeaks2] = findpeaks(inputsignal,'minpeakdistance',options.minTimeBtEvents);
            testpeaks = intersect(testpeaks,testpeaks2);
            clear testpeaks2
            testpeaks = intersect(testpeaks,find(filtfilt(ones(1,options.movAvgReqSize)/options.movAvgReqSize,1,...
                medSubTraces(c,:))>thresh(c)));
        else
            testpeaks=[];
        end

        % if there are none, move on to next cell
        if isempty(testpeaks)
            eventTimes{c} = zeros(1,0);
            riseTimes{c} = zeros(1,0);
            offsetpeaks{c} = zeros(1,0);
            continue
        end

        % find the trough, and thus the rise time and rise amplitude
        [theseRiseTimes,theseRiseHeights] = findRecentTroughs(inputsignal,inputsignal,testpeaks);
        okpeaks = theseRiseHeights > thresh(c); % get the peaks with increase greater than thresh

        peakTimes{c} = testpeaks(okpeaks);
        if options.reportMidpoint
            eventTimes{c} = round(1/2*(testpeaks(okpeaks) + theseRiseTimes(okpeaks)));
            eventTimes{c}(eventTimes{c}<1)=1;
        else
            eventTimes{c} = testpeaks(okpeaks);
        end
        riseTimes{c} = theseRiseHeights(okpeaks);

    end

    % what check is this doing?
    for cellNum=1:size(cellTraces,1)
        thesePeaks=offsetpeaks{cellNum};
        theseRises=riseTimes{cellNum};
        thesePeaks=thesePeaks-options.offsetFrames;
        riseTimes{cellNum}=theseRises(thesePeaks>0);
        offsetpeaks{cellNum}=thesePeaks(thesePeaks>0);
    end
    eventAmps=offsetpeaks;

    if options.displayPeaks
        figure;
        for cInd=1:size(cellTraces,1)
            plot(cellTraces(cInd,:)); hold on; plot(eventTimes{cInd}, cellTraces(cInd,eventTimes{cInd}), 'r.')
            xlim([0 size(cellTraces,2)])
            pause(1)
            hold off
        end
    end

    eventBin=zeros(size(cellTraces));
    for cInd=1:size(cellTraces,1)
        eventBin(cInd,eventTimes{cInd})=1;
    end
end