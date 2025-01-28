function [output] = getTracesAndEvents(extractionFile, decisions, calcium, logs, behavior, options)
%%
if nargin < 6
    options = struct();
end

disp(repmat(' ',[1,100]));
disp(repmat(' ***',[1,20]));
disp(' (1) DeepLabCut');
disp(' (2) CARPA');
disp(' (3) NO Behavior required');
%  wbt = 1;
wbt = input(' Which Behaviours do you want ? (choose one of the above) : ');
disp([' @ Option ',num2str(wbt),' is selected']);
disp(repmat(' ',[1,100]));
TR = true;
while TR
    if wbt==1
        What_Behav_traj = 'DeepLabCut';
        TR = false;
    elseif wbt==2
        What_Behav_traj = 'CARPA';
        TR = false;
    elseif wbt==3
        What_Behav_traj = 'NOBEHAV';
        TR = false;        
    else
        wbt = input(' Which Behaviours do you want ? (choose one of the above) : ');
        TR = true;
    end
end
disp(repmat(' ***',[1,20]));

% disp([' @@@@@@ for detecting behavioural positions, Tracker of ', What_Behav_traj, ' is being used @@@@@@ ']);
% disp([' @@@@@@ If needed you can go "getTracesAndEvents" (Line 7-8) and change it Manually'])

logs_true_length = 0;
for k = 1:length(logs)
    if iscell(logs{k})
        if length(logs{k})>1 %len == 1 indicates logs is up for deletion
            logs_true_length = logs_true_length + 1;
        end
    else
        logs_true_length = logs_true_length + 1;
    end
end
if logs_true_length ~= length(behavior)
    error('ERROR: DIFFERENT NUM OF LOGS THAN BEHAVIORAL FILES')
end

%First we need to separate the extraction file into parts to interpolate
%correctly:
missingFrames = cell([1,length(logs)]);
totalFrames = zeros([1,length(logs)]);
delete_ca_blocks = [];
for k = 1:length(logs)
    if iscell(logs{k})
        if length(logs{k}) == 1
            delete_ca_blocks = [delete_ca_blocks,k];
        end
        all_missing =  [];
        all_total = 0;
        for k_merged = 1:length(logs{k})
            [this_missing, this_total] = readFramesFromLog(logs{k}{k_merged});
            all_missing = [all_missing,this_missing];
            all_total = all_total + this_total;
        end
        missingFrames{k} = all_missing;
        totalFrames(k) = all_total;
    else
        [missingFrames{k}, totalFrames(k)] = readFramesFromLog(logs{k});
    end
    
    if any(missingFrames{k}<0)
       missingFrames{k}(missingFrames{k}<0) = [];
       warning(['Negative frames in ',logs{k},' aborting file'])
    end
end
%%
%Set fix mode
fixMode = 'shift';
askFlag = 1;
%%
%Load calcium file
hinf = hdf5info(calcium);
fullCalcium = hdf5read(hinf.GroupHierarchy.Datasets);
%%
%Get extraction method and cellImages.
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
%%
%Create variables
switch method
    case 'em'
        output.rawProb = cell([1,length(logs_true_length)]);
        output.rawTraces = cell([1,length(logs_true_length)]);
    case 'ica'
        output.rawTraces = cell([1,length(logs_true_length)]);
end
output.position = cell([1,length(logs_true_length)]);
output.velocity = cell([1,length(logs_true_length)]);
output.behavIdx = cell([1,length(logs_true_length)]);
output.caIdx = cell([1,length(logs_true_length)]);
%%
%Load decisions
if isempty(decisions)
    decisions = logical(ones([1,size(cellImages,3)]));
else
    decisions(decisions == 3) = 0; %remove 3 (undefined);
    decisions = logical(decisions);
end
%%

%% Iterate trough the frame blocks of each separate file
blockFrames = [0,cumsum(totalFrames)];
if length(fullCalcium) < blockFrames(end)
    warning('Calcium frames in logs exceeds actual frames')
    blockFrames = blockFrames(blockFrames < length(fullCalcium));
    blockFrames = [blockFrames,length(fullCalcium)];
end
beh_block = 0;
for ca_block = 1:length(totalFrames)
    if any(delete_ca_blocks == ca_block)
        continue
    end
    beh_block = beh_block + 1;

    %%  Load calcium block
    calcium = fullCalcium(:,:,(1+blockFrames(ca_block)):(blockFrames(ca_block+1)));
    
    %% Find behavior for block    
    if strcmp(What_Behav_traj  , 'CARPA')
        disp(['Computing behavior for block ',num2str(beh_block),'\',num2str(length(behavior))])
        if iscell(behavior{beh_block})
            position = [];
            velocity = [];
            for k_beh = 1:length(behavior{beh_block})
                [this_position,this_velocity] = carpaUtilities.compute_trajectory(behavior{beh_block}{k_beh},options);
                position = [position,this_position];
                velocity = [velocity,this_velocity];
            end
        else
            [position,velocity] = carpaUtilities.compute_trajectory(behavior{beh_block},options);
        end
    elseif strcmp(What_Behav_traj  , 'DeepLabCut')
        Animal_name_ind = findstr(behavior{ca_block} , '\Mouse-');
        Animal_name = behavior{ca_block}(Animal_name_ind(1)+7:Animal_name_ind(1)+10);
        date_ind = findstr(behavior{ca_block} , [Animal_name,'-']);
        date =  behavior{ca_block}(date_ind(end)+5:date_ind(end)+12);
        disp(['Loading DeepLabCut behavior for block ',num2str(beh_block),'\',num2str(length(behavior))])
        if ~isempty(findstr(behavior{ca_block} , 'GlobalRemapping')) && ca_block==1
            SESSION = 'GR_A';
        elseif ~isempty(findstr(behavior{ca_block} , 'GlobalRemapping')) && ca_block==2
            SESSION = 'GR_B';
        elseif ~isempty(findstr(behavior{ca_block} , 'NOLTraining'))
            SESSION = 'NOLTraining';
        elseif ~isempty(findstr(behavior{ca_block} , 'NOLTest'))
            SESSION = 'NOLTest';
        end
        bName = ['DeepLabCut_Position_Mouse_', Animal_name,'_',date,'_',SESSION,'_Head.mat'];
        zz = strfind(behavior{ca_block},'\');
        fullpath = behavior{ca_block}(1:zz(end));
        Deeplabcut_File = findfiles(bName , fullpath);
        load(char(Deeplabcut_File)); 
        position = Position; velocity = [];
    elseif strcmp(What_Behav_traj  , 'NOBEHAV')
        position = [];
        velocity = [];
    end
    
    %% Compute upscaled traces
    disp(['Computing upscaled traces for block ',num2str(ca_block),'\',num2str(length(totalFrames))])
    robustMovie = calcium;
    robustMovie(isnan(calcium)) = min(calcium(:));
    switch method
        case 'em'
            %%%upscale em%%
            dsRatio = size(fullCalcium,3)/size(extractionFile.dsCellTraces,2);
            downsampledTraces = extractionFile.dsCellTraces(decisions,floor(1+blockFrames(ca_block)/dsRatio):floor(blockFrames(ca_block+1)/dsRatio)); %Downsampled block
            [upScaledPhi, upFilteredTraces] = recalcPhiAndDetectEvents(robustMovie,cellImages(:,:,decisions),downsampledTraces,extractionFile.CELLMaxoptions,'runEventDetection',0);
            %Interpolate missing frames
            if ~isempty(missingFrames{ca_block})
                fprintf('Interpolating %d missing frames: \n',length(missingFrames{ca_block}));fprintf('%d ',missingFrames{ca_block});fprintf('\n')
                upScaledPhi = interpolMissingFrames(upScaledPhi,missingFrames{ca_block});
                upFilteredTraces = interpolMissingFrames(upFilteredTraces,missingFrames{ca_block});
            end
            upScaledPhi = upScaledPhi*10^5; %revise
            rawProb = upScaledPhi';
            rawTraces = upFilteredTraces';
        case 'ica'
            %%%upscale ica%%
            flTraces = extractFTfromMovie(robustMovie,cellImages(:,:,decisions));
            if ~isempty(missingFrames{ca_block})
                fprintf('Interpolating %d missing frames: \n',length(fullMissing));fprintf('%d ',fullMissing);fprintf('\n')
                flTraces = interpolMissingFrames(flTraces,missingFrames{ca_block});
            end
            rawTraces = flTraces';
    end
    %% For some cases
    if strcmp(What_Behav_traj  , 'NOBEHAV')
        position = zeros(size(rawTraces,1) , 2);
    end
    %%
    %% Fix calcium/behavior differences
    out_behavIdx = 1:size(position,1);
    out_caIdx = 1:size(rawTraces,1);
    if size(rawTraces,1) ~= size(position,1)
        mismath = num2str(abs(size(rawTraces,1)-size(position,1)));
        if strcmp(What_Behav_traj  , 'DeepLabCut')
            disp([' @@@ ',date , ' ' , SESSION , ' behaviour and calcium mismatch (',mismath,' frames)']);
        end
    try
        disp(['Fixing calcium/behavior differences for block ',num2str(ca_block),'\',num2str(length(totalFrames))])
        
        if ~strcmp(fixMode,'manual') && ~strcmp(fixMode,'shift')
            if askFlag
                askFlag = 0;
                choice = questdlg('Use same method next session?', 'Help!','Yes','No','Yes, all sessions','No and stop asking','No'); 
                switch choice
                    case 'Yes'
                        askFlag = 1;
                    case 'Yes, all sessions'
                        askFlag = 0;
                    case 'No'
                        fixMode = 'manual';
                        askFlag = 1;
                    case 'No and stop asking'
                        fixMode = 'manual';
                        askFlag = 0;
                end                           
            end
        end
        
        [ ~,~,behavIdx,caIdx,fixMode] = fixBehaviorCalciumLength(position,rawTraces,fixMode);
        switch method
            case 'em'
                if length(caIdx) > size(rawProb,1)
                    rawProb = fixVecInterpolate(rawProb,length(caIdx));
                else
                    rawProb = rawProb(caIdx,:);
                end
        end
        
        if length(caIdx) > size(rawTraces,1)
            rawTraces = fixVecInterpolate(rawTraces,length(caIdx));
        else
            rawTraces = rawTraces(caIdx,:);            
        end
        
        if length(behavIdx) > size(position,1)
            position = fixVecInterpolate(position,length(behavIdx));
            velocity = fixVecInterpolate(velocity,length(behavIdx));
        else
            position = position(behavIdx,:);
            velocity = velocity(behavIdx,:);
        end
        out_caIdx = nan([1,length(out_caIdx)]);
        out_caIdx(caIdx) = caIdx;
        out_behavIdx = nan([1,length(out_behavIdx)]);
        out_behavIdx(behavIdx) = behavIdx;
    catch        
        if length(rawTraces)> length(position)
            rawTraces = rawTraces(1:length(position),:);
        else
            position = position(1:length(rawTraces),:);
        end
    end
    end
    
    %% Save variables
    output.position{beh_block} = position;
    output.velocity{beh_block} = velocity;
    output.rawTraces{beh_block} = rawTraces;
    output.caIdx{beh_block} = out_caIdx;
    output.behavIdx{beh_block} = out_behavIdx;
    switch method
        case 'em'
            output.rawProb{beh_block} = rawProb;
    end
end
%%
%% Concat variables
output.rawTraces = cat(1,output.rawTraces{:}); 
output.position = cat(1,output.position{:});
output.velocity = cat(1,output.velocity{:});
switch method
    case 'em'
        output.rawProb = cat(1,output.rawProb{:}); 
end

%%
%Get variable to compute events
switch method
    case 'ica' %
        tracesToGetEvents = output.rawTraces;
    case 'em' % cellmax
        tracesToGetEvents = output.rawProb;
end

%% Interpolating values to replace NaNs
IndNans=find(isnan(tracesToGetEvents(:)));
if ~isempty(find(isnan(tracesToGetEvents(:))))
    for k = 1:size(tracesToGetEvents,2)% Number of cells
        IndNans=find(isnan(tracesToGetEvents(:,k)));
        if isempty(IndNans); continue; end
        %% When the first element is NaN
        if ismember(1,IndNans)
            TT = true; i=1;
            while TT
                if ismember(i+1,IndNans); i=i+1; else; TT = false; end
            end
            tracesToGetEvents(1:i , k) = tracesToGetEvents(i+1,k);
        end
        %% When the last element is NaN
        if ismember(size(tracesToGetEvents,1),IndNans)
            TT = true;
            i=size(tracesToGetEvents,1);
            while TT
                if ismember(i-1,IndNans); i=i-1; else; TT = false; end
            end
            tracesToGetEvents(end:-1:i , k) = tracesToGetEvents(i-1,k);
        end
        %% When other elements are NaN
        IndNans=find(isnan(tracesToGetEvents(:,k)));
        if isempty(IndNans); continue; end
        IndListJumps=diff(IndNans');
        Change_points = find([IndListJumps inf]>1);
        Length_of_Sequences = diff([0 Change_points]); % length of the sequences
        EndPoints = cumsum(Length_of_Sequences); % endpoints of the sequences
        StartPoint = 1;
        for kk =1:length(EndPoints)
            A1 = tracesToGetEvents(IndNans(StartPoint)-1,k);
            A2 = tracesToGetEvents(IndNans(EndPoints(kk))+1,k);
            STEP = (A2-A1)/(Length_of_Sequences(kk)+1);
            if A1==A2
                A = ones(1,2+Length_of_Sequences(kk))*A1;
            else
                A = A1:STEP:A2;
            end
            tracesToGetEvents(IndNans(StartPoint):IndNans(EndPoints(kk)),k) = A(2:end-1);
            StartPoint = EndPoints(kk)+1;
        end
    end
end

%% Normalizing the tracesToGetEvents by 99th non-zero percentile
tracesToGetEvents_nonzero = tracesToGetEvents;
tracesToGetEvents_nonzero(tracesToGetEvents==0)=nan;
for i=1:size(tracesToGetEvents,2) % Number of cells
    tracesToGetEvents(:,i) = tracesToGetEvents(:,i)/prctile(tracesToGetEvents_nonzero(:,i),99);
end

disp('Computing events...')
%% COMPUTE EVENTS WITH (1) TresholdEvents
[~, signalPeakIdx] = computeSignalPeaks(tracesToGetEvents','numStdsForThresh',2.3);
tresholdEvents = zeros(size(tracesToGetEvents));
for k = 1:size(tresholdEvents,2)
    tresholdEvents(signalPeakIdx{k},k) = 1;
end
% output.tresholdEvents = tresholdEvents;

%% COMPUTE EVENTS WITH (2)OASIS and (3)MLSPIKES
% hardTresh = .1;
load MLspikesParameters
spikeDeconvTrace = nan(size(tracesToGetEvents));
spikeDeconv = nan(size(tracesToGetEvents));
spikeML = zeros(size(tracesToGetEvents));

for k = 1:size(spikeDeconvTrace,2)
    try
        [c, s , O] = deconvolveCa(tracesToGetEvents(:,k),'exp2','thresholded');
    catch        
        warning(['Error processing traces with oasis in cell ', num2str(k)]);
        continue
    end
    hardTresh = 5*O.sn;
    spikeDeconvTrace(:,k) = c;
    sTreshIdx = find(s > hardTresh);
    sTresh = s(sTreshIdx);
    spikeDeconv(sTreshIdx,k) = round(sTresh/hardTresh);
    
    par.drift.parameter = 0.1*O.sn;
    par.F0 = [0 , prctile(tracesToGetEvents(:,k),90)];    
    spikeML(:,k) = double(myBackward_driftstate(tracesToGetEvents(:,k),par));
end

output.spikeDeconvTrace = spikeDeconvTrace;
output.spikeDeconv = spikeDeconv;
output.spikeML = spikeML;

%% RETURN EACH CELL PEAK%
disp('Computing locations of peaks...')
cellAnatomicLocat = zeros([size(cellImages,3),2]);
for k = 1:size(cellImages,3)
    [xtemp,ytemp] = find(cellImages(:,:,k)==max(max(cellImages(:,:,k))),1);
    cellAnatomicLocat(k,:) = [xtemp,ytemp];
end
output.cellAnatomicLocat = cellAnatomicLocat;

end

