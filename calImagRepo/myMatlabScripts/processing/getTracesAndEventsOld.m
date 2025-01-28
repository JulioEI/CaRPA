function [output] = getTracesAndEvents(extractionFile, decisions, calcium, missing, behavior)

if ~isempty(behavior)
    if ~iscell(behavior)
        behavior = {behavior};
    end
    output.position = [];
    output.velocity = [];
    for behaviorFile = behavior
        try
            [position,velocity,score] = getMouseTrajectory(behaviorFile{1});
        catch
            score = 1;
        end
        if score > .1
            disp('First tracker with grave errors, trying python...')
            disp(behaviorFile{1})
            [positionPy,velocityPy,scorePy] = pyTrackerInterface(behaviorFile{1});
            if scorePy > .1
                disp('First and second trackers with grave errors, taking best one...')
            end
            if scorePy <= score
                position = positionPy;
                velocity = velocityPy;
            end
        end
        output.position = [output.position;position];
        output.velocity = [output.velocity;velocity];
    end
end
 
if ~isempty(missing) %No missing frames specified
    if ~iscell(missing)
        missing = {missing};
    end
    fullMissing = [];
    for missingFramesC = missing
        missingFrames = missingFramesC{1};
        if isstr(missingFrames) %Missing frames specified as txt file
            missingFrames = readFramesFromLog(missingFrames);
            if any(missingFrames<0)
               warning(['Negative frames in ',missingFramesC{1},' aborting file'])
               output = [];
               return
            end
        end
        fullMissing = [fullMissing,missingFrames];
    end
end

if isstr(calcium)
    hinf = hdf5info(calcium);
    calcium = hdf5read(hinf.GroupHierarchy.Datasets);
end

try
    cellImages = extractionFile.cellImages;
    method = 'em';
catch
    try
        cellImages = extractionFile.IcaFilters;
        method = 'ica';
    catch
        error('ExtractionFile not understood')
    end
end

if isempty(decisions)
    decisions = logical(ones([1,size(cellImages,3)]));
else
    decisions(decisions == 3) = 0; %remove 3 (undefined);
    decisions = logical(decisions);
end

%COMPUTE UPSCALED TRACES%
disp('Computing upscaled traces...')
robustMovie = calcium;
robustMovie(isnan(calcium)) = min(calcium(:));
switch method
    case 'em'
        %%%upscale em%%
        [upScaledPhi, upFilteredTraces] = recalcPhiAndDetectEvents(robustMovie,cellImages(:,:,decisions),extractionFile.dsCellTraces(decisions,:),extractionFile.CELLMaxoptions,'runEventDetection',0);
        %Interpolate missing frames
        if ~isempty(fullMissing)
            fprintf('Interpolating %d missing frames: \n',length(fullMissing));fprintf('%d ',fullMissing);fprintf('\n')
            upScaledPhi = interpolMissingFrames(upScaledPhi,fullMissing);
            upFilteredTraces = interpolMissingFrames(upFilteredTraces,fullMissing);
        end
        upScaledPhi = upScaledPhi*10^5; %revise
        output.rawProb = upScaledPhi';
        output.rawTraces = upFilteredTraces';
        output.zTraces = zscore(upFilteredTraces');
        output.zProb = zscore(upScaledPhi');
        tracesToGetEvents = output.rawProb; %From which trace value are the events going to be extracted
    case 'ica'
        %%%upscale ica%%
        flTraces = extractFTfromMovie(robustMovie,cellImages(:,:,decisions));
        if ~isempty(fullMissing)
            fprintf('Interpolating %d missing frames: \n',length(fullMissing));fprintf('%d ',fullMissing);fprintf('\n')
            flTraces = interpolMissingFrames(flTraces,fullMissing);
        end
        output.rawTraces =  flTraces';
        output.zTraces = zscore(flTraces');
        tracesToGetEvents = output.rawTraces; %From which trace value are the events going to be extracted
    otherwise
        error('Method not supported')
end

%FIX CALCIUM/BEHAVIOR DIFFERENCES
if size(tracesToGetEvents,1) ~= size(output.position,1)
    disp('Fixing calcium/behavior differences')
    manual = 1;    
    switch method
        case 'em'
            [ ~,~,behavIdx,caIdx] = fixBehaviorCalciumLength(output.position,output.rawProb,manual);
            output.rawProb = output.rawProb(caIdx,:);
            output.rawTraces = output.rawTraces(caIdx,:);
            output.zTraces = output.zTraces(caIdx,:);
            output.zProb = output.zProb(caIdx,:);
        case 'ica'
            [ ~,~,behavIdx,caIdx] = fixBehaviorCalciumLength(output.position,output.rawTraces,manual);
            output.rawTraces = output.rawTraces(caIdx,:);
            output.zTraces = output.zTraces(caIdx,:);
        otherwise
            error('Method not supported')
    end
    tracesToGetEvents = tracesToGetEvents(caIdx,:);
    output.position = output.position(behavIdx,:);
    output.velocity = output.velocity(behavIdx,:);
end
disp('Computing events...')
%COMPUTE EVENTS BY TREHSOLDING%
tresholdEvents = zeros(size(tracesToGetEvents));
robustInputSignals = tracesToGetEvents';
robustInputSignals(isnan(robustInputSignals)) = 0;
[~, signalPeakIdx] = computeSignalPeaks(robustInputSignals,'numStdsForThresh',2);
for k = 1:size(tresholdEvents,2)
    tresholdEvents(signalPeakIdx{k},k) = 1;
end
output.tresholdEvents = tresholdEvents;

%COMPUTE EVENTS WITH OASIS%
hardTresh = .1;
spikeDeconvTrace = zeros(size(tracesToGetEvents));
spikeDeconv = zeros(size(tracesToGetEvents));
for k = 1:size(spikeDeconvTrace,2)
    [c, s] = deconvolveCa(tracesToGetEvents(:,k),'exp2','thresholded');
    spikeDeconvTrace(:,k) = c;
    sTreshIdx = find(s > hardTresh);
    sTresh = s(sTreshIdx);
    spikeDeconv(sTreshIdx,k) = round(sTresh/hardTresh);
end
output.spikeDeconvTrace = spikeDeconvTrace;
output.spikeDeconv = spikeDeconv;

%COMPUTE EVENTS WITH MLSPIKES
load MLspikesParameters
spikeML = zeros(size(tracesToGetEvents));
for k = 1:size(spikeML,2)
    spikeML(:,k) = double(myBackward_driftstate(tracesToGetEvents(:,k),par));
end
output.spikeML = spikeML;

%RETURN EACH CELL PEAK%
disp('Computing locations of peaks...')
cellAnatomicLocat = zeros([size(cellImages,3),2]);
for k = 1:size(cellImages,3)
    [xtemp,ytemp] = find(cellImages(:,:,k)==max(max(cellImages(:,:,k))),1);
    cellAnatomicLocat(k,:) = [xtemp,ytemp];
end
output.cellAnatomicLocat = cellAnatomicLocat;
end

