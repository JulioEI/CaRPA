classdef PCAnalysis < handle

    properties
        %Folder structures
        rootFolder = '';  
        dayStruct = '';
        mouse = '';
        folderName = '';
        
        %Analysis parameters
        dtCamera = 0.05; %period of miniscope (seconds/frame)
        dT = 0.5; %Integration parameter (seconds)
        binN = [20,1]; %Number of bins in X Y
        minVel = [4, 0]; %Frames under this speed will be discarded (cm/s)
        pxPerCm = [6.5, 6.5]; %How many px a centimeter is (px/cm)
        traceField = 'rawProb'; %Default field to extract traces
        spikeField = 'spikeDeconv'; %Default field to extract spikes
        filterPF = false; %Remove cells with non significant place fields
    end
    
    properties (Constant)
        %File regexp
        tracesEventsRegexp = {'TracesAndEvents'};
        analysisFileRegexp = {'emAnalysis_?\d*.mat$','pcaicaAnalysis_?\d*.mat$'};
        decisionsRegexp = {'Sorted_?\d*.mat$','decisions_?\d*.mat$'}
        
        %Time parsers
        dateParser = 'yyyymmdd'; %Reads and prints time with this format
        timeParser = 'HHMMSS';
    end
    
    methods
        
        function obj = PCAnalysis(folder)
            
            if nargin < 1
                obj.rootFolder = uigetdir('','Select the animal to analyze: ');
            else
                obj.rootFolder = folder;
            end
            
            nameSep = strsplit(obj.rootFolder,filesep);nameSep = nameSep(~cellfun(@isempty,nameSep));
            obj.folderName = nameSep{end};
            regexpOut = regexp(obj.folderName, 'ouse\D*(\d+)','tokens');
            try
                obj.mouse = regexpOut{1}{1};
            end
            obj.getAnalyzableSessions;
        end

        function sessions = chooseSessions(obj,sessionVar)
            %Enter sessions with cell e.g.{1,{1,2},3}means the
            %first and second sessions of day 1 and all the sessions of day
            %3 or use the gui
            sessions = [];
            if nargin < 2 || isempty(sessionVar)
                %If no sessions are specified bring up the gui
                dayList = cellfun(@(x) datestr(x,obj.dateParser),{obj.dayStruct.date},'UniformOutput',0);
                s = listdlg('PromptString','Select day(s):',...
                    'SelectionMode','multiple',...
                    'ListString',dayList);
                for si = s
                    if length(obj.dayStruct(si).sessions) == 1
                        sessions = [sessions,obj.dayStruct(si).sessions];
                    elseif length(obj.dayStruct(si).sessions) > 1
                        timeList = {};
                        for fileN = {obj.dayStruct(si).sessions(:).tracesEventsFileName}
                            regexpOut = regexp(fileN{1}, '(\d+&\d+)','tokens');
                            if isempty(regexpOut)
                                regexpOut = regexp(fileN{1}, '(\d+)','tokens');     
                            end
                            timeList = [timeList,regexpOut{end}{1}];
                        end

                        s2 = listdlg('PromptString',['Select session(s) for day: ',dayList(si)],...
                            'SelectionMode','multiple',...
                            'ListString',timeList);
                        for s2i = s2
                            sessions = [sessions,obj.dayStruct(si).sessions(s2i)];
                        end
                    else
                        warning('The selected day has no avalaible sessions. Run fileSolution to fix')
                    end
                end
            elseif isstruct(sessionVar)
                sessions = sessionVar;
            else
                if isscalar(sessionVar);sessionVar = {sessionVar};end
                for k = 1:length(sessionVar)
                    if iscell(sessionVar{k})
                        sessIdx = [sessionVar{k}{:}]; %Take selected sessions
                        for i = sessIdx
                            sessions = [sessions,obj.dayStruct(sessionVar{k-1}).sessions(i)];
                        end
                    else
                        if k ~= length(sessionVar)
                            if iscell(sessionVar{k+1})
                                continue;
                            end
                        end
                        sessIdx = 1:length(obj.dayStruct(sessionVar{k}).sessions);
                        for i = sessIdx
                            sessions = [sessions,obj.dayStruct(sessionVar{k}).sessions(i)];
                        end
                    end
                end
            end
        end
        
        function turnoverAlignmentMat = quantifyTurnoverAlignment(obj,varargin)
            %OPTIONAL INPUTS: filterPF, sessions, rField, binN, minVel, pxPerCm, dtCamera, dT, 
            
            %Get parameters
            param = obj.getDefaultParams({'binN','minVel','pxPerCm','dtCamera','dT','filterPF'},varargin);
                        
            %Get the sessions
            sessions = obj.chooseSessions(PCAnalysis.parseInput({varargin,'sessions',[]}));
            
            %Remove cells with non significant place fields
            if length(param.filterPF) > 1
                PCdecisions = param.filterPF;
            elseif param.filterPF
                [pos,r] = obj.getXR(varargin{:});
                [~,~,PCdecisions] = PCAnalysis.isPF(r,pos,param,varargin{:});                
            else
                PCdecisions = [];
            end
            
            %Get alignment for sessions
            globalIDs = PCAnalysis.alignSessions(sessions,PCdecisions);
            
            
            turnoverAlignmentMat = zeros(length(sessions));
            for i = 1:length(sessions)
                for j = 1:length(sessions)
                    %Find the cells in current session that align to the reference session
                    ijGlobalID = find(globalIDs(:,i) & globalIDs(:,j));
                    %Ratio of cells that align to reference session i
                    turnoverAlignmentMat(i,j) = length(ijGlobalID)/sum(globalIDs(:,i)~=0);
                end
            end
            
        end
        
        function autoCorrMat = PFautoCorr(obj,varargin)
            %OPTIONAL INPUTS: filterPF, sessions, rField, binN, minVel, pxPerCm, dtCamera, dT, 
            
            %Get parameters
            param = obj.getDefaultParams({'binN','minVel','pxPerCm','dtCamera','dT','filterPF'},varargin);
            
            %Get the sessions
            sessions = obj.chooseSessions(PCAnalysis.parseInput({varargin,'sessions',[]}));
            
            %Get traces and positions
            [pos,r] = obj.getXR(varargin{:});
            
            %Remove cells with non significant place fields
            if length(param.filterPF) > 1
                PCdecisions = param.filterPF;
            elseif param.filterPF
                [r,pos,PCdecisions] = PCAnalysis.isPF(r,pos,param,varargin);               
            else
                PCdecisions = cellfun(@(x) ones([1,size(x,2)]),r,'UniformOutput',false);
            end
            
            %Get alignment for sessions
            globalIDs = PCAnalysis.alignSessions(sessions,PCdecisions);
            
            for k = 1:length(sessions)
                
                halfT = round(length(r{k})/2);
                %Compute PF for first half of the session
                PF.firstHalf = PCAnalysis.buildPF(r{k}(1:halfT,:), pos{k}(1:halfT,:), param);
                %Compute PF for second half of the session
                PF.secondHalf = PCAnalysis.buildPF(r{k}(halfT+1:end,:), pos{k}(halfT+1:end,:), param); 
                       
                %Compute autocorrelation between first and second halfs
                autoCorr{k} = zeros([1,size(r{k},2)]);
                for cellK = 1:size(r{k},2)
                    autoCorr{k}(cellK) = corr2(PF.firstHalf(:,:,cellK),PF.secondHalf(:,:,cellK));
                end
            end           
            
            %Get matrix for all sessions
            PCidx = cellfun(@find, PCdecisions, 'UniformOutput', 0);
            autoCorrMat = zeros(size(globalIDs));
            for cellK = 1:size(globalIDs,1)
                for k = 1:length(sessions)
                    cellInThisSess = globalIDs(cellK,k);
                    if cellInThisSess
                        autoCorrMat(cellK,k) = autoCorr{k}(PCidx{k} == cellInThisSess);
                    end
                end
            end
            
            %Show results
            figure;imagesc(autoCorrMat);colorbar;title('AutoCorr')
            alphaMap = (autoCorrMat == 0);
            white = cat(3, ones(size(autoCorrMat)), ones(size(autoCorrMat)), ones(size(autoCorrMat))); 
            hold on; h = imagesc(white); hold off; set(h, 'AlphaData', alphaMap);
            
        end

        function [xCorrMat, meanCorrelation,steCorrelation] = PFxCorr(obj,varargin)
            %OPTIONAL INPUTS: sessions, rField, binN, minVel, pxPerCm, dtCamera, dT, 

            %Get parameters
            param = obj.getDefaultParams({'binN','minVel','pxPerCm','dtCamera','dT','filterPF'},varargin);
            
            %Get the sessions
            sessions = obj.chooseSessions(PCAnalysis.parseInput({varargin,'sessions',[]}));
            
            %Get traces and positions
            [pos,r] = obj.getXR(varargin{:});
            
            %Remove cells with non significant place fields
            if param.filterPF
                [r,pos,PCdecisions] = PCAnalysis.isPF(r,pos,param,varargin);
            else
                PCdecisions = cellfun(@(x) ones([1,size(x,2)]),r,'UniformOutput',false);
            end
            
            %Get alignment for sessions
            globalIDs = PCAnalysis.alignSessions(sessions,PCdecisions);

            %Build all PF
            allPF = arrayfun(@(i) PCAnalysis.buildPF(r{i}, pos{i}, param),1:length(sessions),'UniformOutput',0);
            
            %Compute PF correlation
            PCidx = cellfun(@find, PCdecisions, 'UniformOutput', 0);
            xCorrMat = zeros([length(sessions),length(sessions),size(globalIDs,1)]);
            for i = 1:length(sessions)
                for j = 1:length(sessions)
                    %Find the cells in current session that align to the reference session
                    ijGlobalID = find(globalIDs(:,i) & globalIDs(:,j));
                    %Compute xcorrelation between reference
                    for k = 1:length(ijGlobalID)
                        cellPFi = allPF{i}(:,:,PCidx{i} == globalIDs(ijGlobalID(k),i));
                        cellPFj = allPF{j}(:,:,PCidx{j} == globalIDs(ijGlobalID(k),j));
                        xCorrMat(i,j,ijGlobalID(k)) = corr2(cellPFi,cellPFj);
                    end
                end
            end
            
            %Show results
            thisCorrMat = squeeze(xCorrMat(1,:,:))';
            
            figure;imagesc(thisCorrMat);colorbar;title('xCorr')
            alphaMap = (thisCorrMat == 0);
            white = cat(3, ones(size(thisCorrMat)), ones(size(thisCorrMat)), ones(size(thisCorrMat))); 
            hold on; h = imagesc(white);hold off; set(h, 'AlphaData', alphaMap);
            ylim([1,max(arrayfun(@(x) find(thisCorrMat(:,x),1,'last'),1:size(thisCorrMat,2)))])
            
            %Show mean results
            meanCorrelation = zeros(length(sessions));
            steCorrelation = zeros(length(sessions));
            for i = 1:length(sessions)
                for j = 1:length(sessions)
                    ijAllCells = squeeze(xCorrMat(i,j,:));
                    meanCorrelation(i,j) = nanmean(ijAllCells(ijAllCells~=0));
                    steCorrelation(i,j) = nanstd(ijAllCells(ijAllCells~=0))/sqrt(length(ijAllCells(ijAllCells~=0)));
                end
            end
            figure;pcolor(meanCorrelation);colorbar;axis square
            
            figure;
            myColormap = lines;
            hold on
            h = shadedErrorBar(2:length(meanCorrelation),diag(meanCorrelation,1),diag(steCorrelation,1),'lineprops',{'-o','color',myColormap(1,:),'linewidth',2},'patchSaturation',0.2);h.mainLine.DisplayName = ' with the previous day';               
            h = shadedErrorBar(1:length(meanCorrelation),meanCorrelation(1,:),steCorrelation(1,:),'lineprops',{'-o','color',myColormap(2,:),'linewidth',2},'patchSaturation',0.2);h.mainLine.DisplayName = ' with the first day';     
            h = shadedErrorBar(1:length(meanCorrelation),meanCorrelation(round(length(meanCorrelation)/2),:),steCorrelation(1,:),'lineprops',{'-o','color',myColormap(3,:),'linewidth',2},'patchSaturation',0.2);h.mainLine.DisplayName = ' with the middle day';   
            h = shadedErrorBar(1:length(meanCorrelation),meanCorrelation(end,:),steCorrelation(end,:),'lineprops',{'-o','color',myColormap(4,:),'linewidth',2},'patchSaturation',0.2);h.mainLine.DisplayName = ' with the last day';     
            legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'));
            grid on;axis square;title('mean PF correlation');
        end
        
        function PF = buildSessionPF(obj,varargin)
            %Builds PF of a session with spikes, traces and big and small time bins
            %OPTIONAL INPUTS: rField, session, binN, minVel, pxPerCm, dtCamera, dT
                        
            %Get parameters
            param = obj.getDefaultParams({'binN','minVel','pxPerCm','dtCamera','dT','traceField','spikeField'},varargin);
            
            %Get position and neuron responses
            [pos,r,nSessions] = obj.getXR(varargin);
            
            %Compute PF
            PF = cell([1,nSessions]);
            for k = 1:nSessions
                [PF{k}.mean,PF{k}.std,PF{k}.ste] = PCAnalysis.buildPF(r{k}, pos{k}, param);    
            end
            if nSessions == 1;PF = PF{1};end
            
        end
                
        function PFinf = computeCellInformationSession(obj,varargin)
            %OPTIONAL INPUTS: session, rField, bigN, smallN, smallShift, bigShift, smallBigDistanceShift, binN, minVel, pxPerCm, dtCamera, dT, 
            
            %Get parameters
            param = obj.getDefaultParams({'binN','minVel','pxPerCm','dtCamera','dT'},varargin);
            
            param.smallN = PCAnalysis.parseInput({varargin,'smallN',10});
            param.bigN = PCAnalysis.parseInput({varargin,'bigN',50});
            param.smallShift = PCAnalysis.parseInput({varargin,'smallShift',param.dtCamera}); %In seconds
            param.bigShift = PCAnalysis.parseInput({varargin,'bigShift',param.dtCamera*10}); %In seconds
            param.smallBigDistanceShift = PCAnalysis.parseInput({varargin,'smallBigDistanceShift',param.dtCamera*100}); %In seconds
            
            %Get responses
            [pos,r,nSessions] = obj.getXR(varargin);
            
            %Compute cells information
            PFinf = cell([1,nSessions]);
            for k = 1:nSessions
                PFinf{k} = PCAnalysis.computeCellinformation(r{k}, pos{k}, param);    
            end
            if nSessions == 1;PFinf = PFinf{1};end
            
        end
        
        function output = loadOutput(obj,sessions)
            if nargin < 2
                sessions = obj.chooseSessions;
            end
            output = cell([1,length(sessions)]);
            for k = 1:length(sessions)
                disp(['...loading ',sessions(k).tracesEventsFileName,' ...'])
                fullPath = fullfile(sessions(k).rootFolder,sessions(k).folderName,sessions(k).tracesEventsFileName);
                output{k} = load(fullPath);dummyField = fields(output{k});output{k} = output{k}.(dummyField{1});%load regardless of struct name 
                %If positions missing, get them from behav
                if ~isfield(output{k},'position')
                    PCAnalysis.appendTrajectoriesToOutput(sessions(k))
                    output{k} = load(fullPath);dummyField = fields(output{k});output{k} = output{k}.(dummyField{1});%load regardless of struct name 
                end
            end
            if length(sessions) == 1;output = output{1};end
        end
        
        function [x,r,nSessions] = getXR(obj,varargin)
            %OPTIONAL INPUTS: session, rField, decideRForEachFile
            
            sessions = obj.chooseSessions(PCAnalysis.parseInput({varargin,'sessions',[]}));
            nSessions = length(sessions);
            
            decideRForEachFile = PCAnalysis.parseInput({varargin,'decideRForEachFile',0});
            
            [x,r] = deal(cell([1,length(sessions)]));
                        
            for k = 1:length(sessions)
                %Get position and responses from output file           
                output = obj.loadOutput(sessions(k));
                
                %Get response field
                if decideRForEachFile || ~exist('rField','var')
                    rField = PCAnalysis.parseInput({varargin,'rField',@() obj.selectOutput(output)});
                    if iscell(rField);rField = rField{k};end
                end

                %Get traces and events and convert spikes into one-hot embedding   
                x{k} = output.position;
                r{k} = PCAnalysis.ind2OneHot(eval(['output.',rField]),length(x));

                %Check that length calcium == length behav
                if length(x{k}) ~= size(r{k},1)
                    warning(['Calcium traces have a length of ',num2str(size(r{k},1)),'. However, behavioral trajectories have a length of ',num2str(length(x{k}))])
                end
            end
        end
                
    end
    
    methods(Static)
        
        function compareSpikeOutputs(output,varargin)
            
            traces = output.(PCAnalysis.parseInput({varargin,'traces','rawProb'}));
            spikes = PCAnalysis.parseInput({varargin,'spikes',{'spikeDeconv','spikeML','tresholdEvents'}});
            
            %Compare spike trains
            figure;
            myColormap = lines;
            for celli = 1:size(traces,2)
               clf;
               subplot(1,2,1);plot(0.05*(1:length(traces)),traces(:,celli));hold on;grid on
               h = plot(nan,nan,'Color',myColormap(1,:));
               subplot(1,2,2);plot(0.05*(1:length(traces)),traces(:,celli));hold on;grid on
               k = 1;
               for spike = spikes
                  myS = output.(spike{1})(:,celli);
                  mySFind = find(myS);mySProb = myS(myS~=0);
                  miniShift = (k-1)*.005; %For visualization purposes
                  subplot(1,2,1)
                  plot(miniShift+0.05*[mySFind,mySFind]',([zeros([1,length(mySProb)])',mySProb])'+repmat([0,0],[length(mySFind),1])','lineWidth',1.7,'Color',myColormap(k+1,:))
                  subplot(1,2,2)
                  plot(0.05*[mySFind,mySFind]',(-[(.1*(k-1))+zeros([1,length(mySProb)])',.1*k+zeros([1,length(mySProb)])'])'+ repmat([0,0],[length(mySFind),1])','lineWidth',1.7,'Color',myColormap(k+1,:))                  
                  h = [h,plot(nan,nan,'Color',myColormap(k+1,:),'lineWidth',1.7)]; %dummy legend
                  k = k + 1;
               end
               subplot(1,2,1);xlabel('Time(s)');
               subplot(1,2,2);xlabel('Time(s)');
               linkaxes;
               legend(h,'Original',spikes{:});
               pause;
            end
        end
        
        function [r,x,decisions] = isPF(r,x,param,varargin)
            %Searchs for significance with entropy of PF against shifts for
            %default shift parameters
            
            smallN = PCAnalysis.parseInput({varargin,'smallN',10});
            bigN = PCAnalysis.parseInput({varargin,'bigN',20});
            smallShift = PCAnalysis.parseInput({varargin,'smallShift',param.dtCamera}); %In seconds
            bigShift = PCAnalysis.parseInput({varargin,'bigShift',param.dtCamera*10}); %In seconds
            smallBigDistanceShift = PCAnalysis.parseInput({varargin,'smallBigDistanceShift',param.dtCamera*100}); %In seconds

            dtCamera = param.dtCamera;
            
            %Build shift vectors
            shiftSlowSec = -smallShift*smallN/2:smallShift:smallShift*smallN/2;
            shiftSlowF = round(shiftSlowSec./dtCamera);
            
            shiftFastSec = smallBigDistanceShift + (0:bigShift:bigN*bigShift);
            shiftFastF = round(shiftFastSec./dtCamera);    
            
            %Compute inf of shifted
            if ~iscell(r)
                r = {r};x = {x};
            end
            decisions = cell([1,length(r)]);
            textprogressbar(['Computing information of sessions... ']);
            for k = 1:length(r)
                pfEntropyCell = zeros([length([shiftSlowF,shiftFastF]),size(r{k},2)]);
                
                textprogressbar(100*k/length(r));
                shiftVect = [shiftSlowF,shiftFastF];
                parfor i = 1:length(shiftVect)                    
                    %Shift responses
                    rShift = circshift(r{k},shiftVect(i),1);

                    %Compute PF
                    PF = PCAnalysis.buildPF(rShift, x{k}, param);
                    normalizedPF = PF./nansum(nansum(PF,2),1);

                    %If any PF is negative, shift everything up
                    [~,~,zi] = ind2sub(size(normalizedPF),find(normalizedPF<0));
                    normalizedPF(:,:,zi) = normalizedPF(:,:,zi) - nanmin(nanmin(normalizedPF(1:2,1,zi),[],2),[],1);

                    %Compute entropy
                    entropy = -squeeze(nansum(nansum(normalizedPF.*log2(normalizedPF),2),1));
                    pfEntropyCell(i,:) = entropy;
                end
                
                %Check for significance with ttest
                decisions{k} = zeros([1,size(r{k},2)]);
                for celli = 1:size(r{k},2)
                    decisions{k}(celli) = ttest2(pfEntropyCell(1:length(shiftSlowSec),celli)',pfEntropyCell(length(shiftSlowSec)+1:end,celli)',0.05,'left');
                end
                r{k} = r{k}(:,logical(decisions{k}));
            end
            textprogressbar(' done');
            disp([num2str(100*sum(cellfun(@sum,decisions))/sum(cellfun(@length,decisions))),'% of cells had significant PF']);
            if length(r) == 1
                r = r{1};x = x{1};decisions = decisions{1};
            end
        end
        
        function info = computeFullCellinformation(r,x,param)         
            %Get parameters
            smallN = param.smallN;
            bigN = param.bigN;
            smallShift = param.smallShift;
            bigShift = param.bigShift;
            smallBigDistanceShift = param.smallBigDistanceShift;
            dtCamera = param.dtCamera;
            
            %Build shift vectors
            shiftSlowSec = -smallShift*smallN/2:smallShift:smallShift*smallN/2;
            shiftSlowF = round(shiftSlowSec./dtCamera);
            
            shiftFastSec = smallBigDistanceShift + (0:bigShift:bigN*bigShift);
            shiftFastF = round(shiftFastSec./dtCamera);
            
            %Compute inf of shifted
            pfEntropyCell = [];
            miCell = [];
            textprogressbar('Computing information... ');
            i = 0;
            for shift = [shiftSlowF,shiftFastF]
                textprogressbar(100*i/length([shiftSlowF,shiftFastF]));i = i + 1;
                %Shift responses
                rShift = circshift(r,shift,1);
                 
                %Compute PF
                PF = PCAnalysis.buildPF(rShift, x, param);
                normalizedPF = PF./nansum(nansum(PF,2),1);
                
                %If any PF is negative, shift everything up
                [~,~,zi] = ind2sub(size(normalizedPF),find(normalizedPF<0));
                normalizedPF(:,:,zi) = normalizedPF(:,:,zi) - nanmin(nanmin(normalizedPF(1:2,1,zi),[],2),[],1);
                
                %Compute entropy
                entropy = -squeeze(nansum(nansum(normalizedPF.*log2(normalizedPF),2),1));
%                 figure;for k = 1:size(normalizedPF,3);bar(normalizedPF(:,1,k));title(-nansum(nansum(normalizedPF(:,1,k).*log2(normalizedPF(:,1,k)),2),1));pause;end
%                 figure;
%                 myShift = [-5,0,100];
%                 for k = 1:(i^2)
%                     subplot(i,i,k);bar(myPF(:,:,k));xlim([1,param.binN(1)])
%                     myLegend = cell([3,1]);
%                     for j = 1:3
%                        myLegend{j} =  ['Shift: ',num2str(myShift(j)),' frames',char(10),'Entropy: ',num2str(-nansum(myPF(:,j,k).*log2(myPF(:,j,k)),1))];
%                     end
%                     legend(myLegend);
%                 end
                pfEntropyCell = [pfEntropyCell,entropy];
                %ComputeMI
                miThisCell = zeros([size(rShift,2),size(x,2)]);
                for celli = 1:size(rShift,2)
                    for dim = 1:size(x,2)
                        miThisCell(celli,dim) = mutualinfo(rShift(:,celli),x(:,dim));
                    end
                end
                miCell = cat(3,miCell,miThisCell);
            end
            textprogressbar(' done');
            info.miCell = miCell;
            info.pfEntropyCell = pfEntropyCell;
            [pFmean,~,pFerror] = PCAnalysis.buildPF(r, x, param);
            info.zeroShiftPF.mean = pFmean;
            info.zeroShiftPF.ste = pFerror;
            info.shiftSlowSec = shiftSlowSec;
            info.shiftFastSec = shiftFastSec;
        end
        
        function viewPFinfo(info)
            
            tempFields = fields(info);
            if ~isstruct(info.(tempFields{1}))
                tempInfo.response = info;
                info = tempInfo;
            end
            myFields = fields(info);
            
            cellN = size(info.(myFields{1}).zeroShiftPF.mean,3); 
            
            fig = figure;
            for celli = 1:cellN
                clf;
                for k = 1:length(myFields)
                    binY = size(info.(myFields{k}).zeroShiftPF.mean,2);
                    pFmean = info.(myFields{k}).zeroShiftPF.mean;
                    pFerror = info.(myFields{k}).zeroShiftPF.ste;
                    pfEntropyCell = info.(myFields{k}).pfEntropyCell;
                    shiftSlowSec = info.(myFields{k}).shiftSlowSec;
                    shiftFastSec = info.(myFields{k}).shiftFastSec;

                    if binY == 1
                        miCellx = squeeze(info.(myFields{k}).miCell(:,1,:));
                        [validInf, pValInf] = ttest2(pfEntropyCell(celli,1:length(shiftSlowSec))',pfEntropyCell(celli,length(shiftSlowSec)+1:end)',0.05,'left');
                        [validMI, pValMI] = ttest2(miCellx(celli,1:length(shiftSlowSec))',miCellx(celli,length(shiftSlowSec)+1:end)',0.05,'right');
                        if validInf; pColorInf = 'green';else;pColorInf = 'red';end
                        if validMI; pColorMI = 'green';else;pColorMI = 'red';end
                        subplot(length(myFields),3,1+(k-1)*3)
                        plot([shiftSlowSec,shiftFastSec],pfEntropyCell(celli,:),'-o','color',pColorInf);xlabel('delay(s)');ylabel('entropy');title(['pvalue ',num2str(pValInf)]);axis square;grid on;ylim([0,inf])
                        subplot(length(myFields),3,2+(k-1)*3)
                        shadedErrorBar(1:size(pFmean,1),pFmean(:,1,celli),pFerror(:,1,celli),'lineprops',{'-o','linewidth',2,'color',pColorInf},'patchSaturation',0.2);title(['cell',num2str(celli)]);axis square
                        subplot(length(myFields),3,3+(k-1)*3)
                        plot([shiftSlowSec,shiftFastSec],miCellx(celli,:),'-o','color',pColorMI);xlabel('delay(s)');ylabel('mutual info');title(['pvalue ',num2str(pValMI)]);axis square;grid on;ylim([0,inf])
                        if length(myFields)>1;annotation('textbox',[0,0,1,(1+1/length(myFields))-k*(1-1/length(myFields))/(length(myFields)-1)],'string',myFields{k},'fontSize',20,'FontWeight','bold');end
                    else
                       error('Not yet implemented')
                    end
                end
                pause;
            end
        end
        
        function [PFmeanRate,PFstdRate,PFsteRate,binCounts] = buildPF(r, x, param)
%             disp('Computing PF...')
%             disp(repmat('.',[1,100]))

            
            %Get parameters
            dT = param.dT; %Temporal bin size (s)
            binN = param.binN; %Number of bins
            dtCamera = param.dtCamera; %Period of camera (s)
            pxPerCm = param.pxPerCm; %Pixels/cm ratio
            minVel = param.minVel; %Minimal speed (cm/s)
            
            %Check for outliers
            outliers = x>(mean(x)+3*std(x));
            if sum(outliers(:,1)) ~= 0 %only for x
                warning('OUTLIERS DETECTED, MANUAL ACTION REQUIERED')
                %pause
                x = x(100:end,:);
                r = r(100:end,:);
            end
            
            %Put pos in cm units:
            x = x./pxPerCm;
            
            %Get absolute velocity
            v = abs(diff(x)/dtCamera);
            
            %Get frames under min velocity
            framesUnderMinVel = unique([find(v(:,1) < minVel(1)); find(v(:,2) < minVel(2))]);
            
            %Remove frames under min velocity
            x(framesUnderMinVel,:) = [];
            r(framesUnderMinVel,:) = [];
            
%             %Get only directional PF (linear track only)
%             framesDirection = find(diff(x(:,1))<0);
%             x = x(framesDirection,:);
%             r = r(framesDirection,:);
            
            %Integrate in time
            [rTB,xTB] = PCAnalysis.integrateInTime(r,x,dT,dtCamera);
%             figure;plot(x,'-*');hold on;plot(linspace(5,length(x)+5,length(xTB)),xTB,'o');title('>4 cm/s');plot(repmat(1:length(x),[20,1])',(ones([20,length(x)]).*linspace(10,110,20)')','k--')
%             figure;plot(linspace(5,length(x)+5,length(xTB)),PCAnalysis.binPosition(xTB,binN),'-o');hold on;plot(PCAnalysis.binPosition(x,binN),'-*');plot(repmat(1:length(x),[20,1])',(ones([20,length(x)]).*linspace(1,20,20)')','k--')
            
            %Bin positions
            [xTB,binCounts] = PCAnalysis.binPosition(xTB,binN);
            
            %Get activity per bin
            PF = cell([max(xTB),size(rTB,2)]);
            for t = 1:length(xTB)
                for celli = 1:size(rTB,2)
                    PF{xTB(t,1),xTB(t,2),celli} = [PF{xTB(t,1),xTB(t,2),celli},rTB(t,celli)];    
                end
            end
            
            %Get mean count and std
            PFstdCount = cellfun(@std,PF);
            PFmeanCount = cellfun(@mean,PF);
            PFsteCount = PFstdCount./sqrt(binCounts);
            
            %Normalize per dT
            PFmeanRate = PFmeanCount./dT;
            PFstdRate = PFstdCount./dT;
            PFsteRate = PFsteCount./dT;
        end
        
        function viewPF(PF, varargin)
            
            normalize = PCAnalysis.parseInput({varargin,'normalize',0});
            
            if isstruct(PF)
                myFields = fields(PF);
                binX = size(PF.(myFields{1}).mean,1);
                binY = size(PF.(myFields{1}).mean,2);
                cellN = size(PF.(myFields{1}).mean,3);
            else
                binX = size(PF,1);
                binY = size(PF,2);
                cellN = size(PF,3);
            end
            
            figure;
            for cell = 1:cellN
                if binY == 1
                    if isstruct(PF)
                        clf;
                        myColormap = lines;
                        i = 1;
                        for myField = myFields'
                            pFmean = PF.(myField{1}).mean(:,1,cell);if normalize;pFmean = pFmean./max(pFmean);end
                            pFerror = PF.(myField{1}).ste(:,1,cell);
                            hold on
                            h = shadedErrorBar(1:length(pFmean),pFmean,pFerror,'lineprops',{'-o','color',myColormap(i,:),'linewidth',2},'patchSaturation',0.2);
                            h.mainLine.DisplayName = myField{1};
                            i = i + 1;
                        end
                        legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'));
                    else
                        plot(PF(:,1,cell),'-o');
                    end
                    title(['Place field of cell ',num2str(cell)])
                    xlim([0,binX+1])
                    pause;
                else
                    error('Not yet implemented')
                end   
            end           
        end
               
        function [rTB,xTB] = integrateInTime(r,x,dT,dtCamera)
            timeBin = dT/dtCamera;
            sizeTB = floor(length(x)/timeBin);
            rTB = zeros([sizeTB,size(r,2)]);
            xTB = zeros([sizeTB,2]);
            for k = 1:sizeTB
                lb = (k-1)*timeBin + 1;
                ub = (k)*timeBin;
                rTB(k,:) = sum(r(lb:ub,:),1);
                xTB(k,:) = mean(x(lb:ub,:),1);
            end            
        end
        
        function [posB,binCounts] = binPosition(position,binN)
            
            %Normalize position
            posN = (position - min(position)) ./ (max(position) - min(position));
            
            %Find bin edges
            [binCounts,binEdgesX,binEdgesY] = histcounts2(posN(:,1),posN(:,2),binN);
            
            %Bin positions
            posB = zeros(size(posN));
            for dim = 1:2
                if dim == 1;dimBinEdges = binEdgesX;else;dimBinEdges = binEdgesY;end
                for k = 1:(length(dimBinEdges)-1)
                    posInBin = find(posN(:,dim)>=dimBinEdges(k) & posN(:,dim)<=dimBinEdges(k+1));
                    posB(posInBin,dim) = k;
                end
            end
        end
        
        function scores = getClassifierConfidence(sessions,PCdecisions)           
           %Run decoder 
           scores = {};
           for k = 1:length(sessions)
                decisions = load(fullfile(sessions(k).rootFolder,sessions(k).folderName,sessions(k).decisionsFileName));
                emFile = load(fullfile(sessions(k).rootFolder,sessions(k).folderName,sessions(k).analysisFileFileName));
                [~,fileName,fileExt] = fileparts(emFile.emAnalysisOutput.movieFilename);
                moviePath = fullfile(sessions(k).rootFolder,sessions(k).folderName,[fileName,fileExt]);
                cInf = cellInfo(moviePath,permute(emFile.emAnalysisOutput.cellImages,[3,1,2]),emFile.emAnalysisOutput.scaledProbability,[]);
                sessionScores = cInf.predictClassifier;
                if nargin < 2 || isempty(PCdecisions)
                    scores{k} = sessionScores(logical(decisions.validCellMax));                    
                else
                    scores{k} = sessionScores(PCdecisions{k}(logical(decisions.validCellMax)));
                end
           end
        end
        
        function globalIDs = alignSessions(sessions, PCdecisions)

            if nargin < 2
                PCdecisions = [];
            end
            
           %Load session cellmaps
           allCellMaps = {};
           for k = 1:length(sessions)
                decisions = load(fullfile(sessions(k).rootFolder,sessions(k).folderName,sessions(k).decisionsFileName));
                emFile = load(fullfile(sessions(k).rootFolder,sessions(k).folderName,sessions(k).analysisFileFileName));
                allCellMaps{k} = emFile.emAnalysisOutput.cellImages(:,:,logical(decisions.validCellMax));
           end
           allCellMaps = cellfun(@(x) thresholdImages(x),allCellMaps,'UniformOutput',0);
           
           %Make sure that all cellmaps are the same size
           maxSize = [max(cellfun(@(x) size(x,1),allCellMaps)),max(cellfun(@(x) size(x,2),allCellMaps))];
           for k = 1:length(allCellMaps)
               maxSizeTemplate = zeros([maxSize,size(allCellMaps{k},3)]);
               maxSizeTemplate(1:size(allCellMaps{k},1),1:size(allCellMaps{k},2),:) = allCellMaps{k};
               allCellMaps{k} = maxSizeTemplate;
           end
           
           %Put all cells in the same image
           allCellMapsMax = cellfun(@(x) max(x,[],3),allCellMaps,'UniformOutput',0);

           %Align cells separately
           templateIdx = round(length(allCellMaps)/2);

           d1 = size(allCellMapsMax{1},1);
           d2 = size(allCellMapsMax{1},2);

           options = NoRMCorreSetParms('d1',d1,'d2',d2,'grid_size',[d1/2,d2/2,1],'bin_width',1,'min_patch_size',[5,5,1],'upd_template',0,'correct_bidir',0,'use_parallel',0,'max_shift',200,'shifts_method','linear');
           template = allCellMapsMax{templateIdx};
           alignedImgNC = cell([1,length(allCellMaps)]);
           for j = 1:length(allCellMapsMax)
               [~,shifts] = normcorre(allCellMapsMax{j},options,template);
               for cellIdx = 1:size(allCellMaps{j},3)
                   alignedImgNC{j} = cat(3,alignedImgNC{j},double(apply_shifts(allCellMaps{j}(:,:,cellIdx),shifts,options)));
               end
           end
           %Find cell centroids
           coords = {};
           for j=1:length(allCellMapsMax)
               [xCoords, yCoords] = findCentroid(alignedImgNC{j});
               coords{j} = [xCoords; yCoords]';
           end
           
           %Run cell correspondence finder (Biafra's)
           globalIDoptions.analysisType = 'pairwise';
           globalIDoptions.trialToAlign = templateIdx;
           globalIDoptions.maxDistance = 5;
           globalIDoptions.nCorrections = 3;
           globalIDoptions.trialIDs = [];
           [OutStruct] = computeGlobalIdsPairwise([],coords,globalIDoptions,length(allCellMapsMax));
           globalIDs = OutStruct.globalIDs;
           
           %Remove unwanted cells
           if ~isempty(PCdecisions)
               for k = 1:length(PCdecisions)
                   removeIdx = setdiff(1:size(globalIDs,1),PCdecisions{k});
                   globalIDs(removeIdx,k) = 0;
%                    if isempty(setdiff(PCdecisions{k},[0,1])) 
%                        PCidx = find(~PCdecisions{k}); %Is a logical
%                    else
%                        PCidx = setdiff(1:max(globalIDs(:,k)),PCdecisions{k}); %Is an idx
%                    end
%                     for PCdecision = PCidx
%                         globalIDs(globalIDs(:,k) == PCdecision,k) = 0;
%                     end
               end
           end
           
        end
        
        function matHot = ind2OneHot(idx,totalSize)
            %Converts cells of indices into one hot encoding matrix
            if iscell(idx)
                matHot = zeros([totalSize,length(idx)]);
                for cell = 1:length(idx)
                    matHot(idx{cell},cell) = 1;
                end
            else
                matHot = idx;
            end
        end
        
        function imagescNanAsWhite(X)
            figure;imagesc(X);colorbar;title('AutoCorr')
            alphaMap = (X == 0);
            white = cat(3, ones(size(X)), ones(size(X)), ones(size(X))); 
            hold on; h = imagesc(white); hold off; set(h, 'AlphaData', alphaMap);
        end
        
    end
    
    methods(Static, Access = private)
               
        function match = checkRegexp(sourceText,regexpList)
            match = zeros(size(sourceText));
            k = 1;
            for text = sourceText
                regexpMatch = 0;
                for regexpSeq = regexpList
                    if regexp(text{1},regexpSeq{1})
                        regexpMatch = regexpMatch + 1;
                    end
                end
                if regexpMatch > 0
                    match(k) = 1;
                end
                k = k + 1;
            end
        end  
        
        function myOutput = selectOutput(output)
            myFields = fields(output);
            s = listdlg('PromptString','Select cell response variable:',...
                'SelectionMode','single',...
                'ListString',myFields);
            myOutput = myFields{s};
        end
        
        function [varargout] = parseInput(varargin)
            %{{inputs},optional(arg1 default1, arg2 default2, ...)}
            tmpVar = varargin{1};
            if iscell(tmpVar{1})
                argumentList = tmpVar(2:2:end);
                defaultVal = tmpVar(3:2:end);
                input = tmpVar{1};
                inputArg = input(1:2:end);
                inputVal = input(2:2:end);
                for k = 1:length(argumentList)
                    idx = strcmp(argumentList{k},inputArg);
                    if sum(idx) == 0
                        varargout{k} = defaultVal{k}();
                    else
                        varargout{k} = inputVal{idx};
                    end
                end                
            else
                for k = 2:2:length(tmpVar)
                    varargout{k} = tmpVar{k};
                end
            end
        end
        
        function appendTrajectoriesToOutput(session)
            %get behavior files
            fullPath = fullfile(session.rootFolder,session.folderName,session.tracesEventsFileName);
            allItems = dir(fullfile(session.rootFolder,session.folderName));
            allItems = allItems(3:end); %remove WINDOWS metadata folders
            allFiles = allItems(~[allItems.isdir]);
            behaviorFiles = cell([1,size(session.date,1)]);
            for k = 1:size(session.date,1)
                thisSess = session.date(k,:);
                day = datestr(thisSess,PCAnalysis.dateParser);
                time = datestr(thisSess,PCAnalysis.timeParser);
                behavFiles = {allFiles(find(session.checkRegexp({allFiles.name},{[day,'.*',time,'.*behavior']}))).name};
                behaviorFiles{k} = fullfile(session.rootFolder,session.folderName,behavFiles{1});
            end
            appendPosToTracesEventsFile(fullPath,behaviorFiles)
        end
        
    end
    
    methods(Access = private)
        
        function param = getDefaultParams(obj,myFields,input)
            for field = myFields
                param.(field{1}) = PCAnalysis.parseInput({input,field{1},obj.(field{1})});
            end
        end
        
        function getAnalyzableSessions(obj)
            allItems = dir(obj.rootFolder);
            allItems = allItems(3:end); %remove WINDOWS metadata folders
            allFolders = allItems([allItems.isdir]);
            i = 1;
            for dayFolder = {allFolders.name}
                regexpOut = regexp(dayFolder{1}, '(\d+)','tokens');
                day = regexpOut{end}{1};%Assumess date is the last sequence of numbers in the folder name
                fullPath = fullfile(obj.rootFolder, dayFolder{1});
                obj.dayStruct(i).path = fullPath;
                obj.dayStruct(i).folderName = dayFolder{1};
                obj.dayStruct(i).date = datevec(day,obj.dateParser);
                obj.dayStruct(i).mouse = obj.mouse;
                
                allInsideitems = dir(fullPath);
                allInsideitems = allInsideitems(3:end); %remove WINDOWS metadata folders
                allFiles = allInsideitems(~[allInsideitems.isdir]);
                tracesEventsIdx = find(obj.checkRegexp({allFiles.name},obj.tracesEventsRegexp));
                if isempty(tracesEventsIdx)
                    warning(['cannot find traces/events files for ' dayFolder{1}])
                else
                    tracesEventsFiles = {allFiles(tracesEventsIdx).name};
                    k = 1;
                    for tracesEventsFile = tracesEventsFiles
                        noExtensionName = split(tracesEventsFile{1},'.');
                        regexpOut = regexp(noExtensionName{1}, '(\d+&\d+)','tokens');
                        if isempty(regexpOut)
                            regexpOut = regexp(noExtensionName{1}, '(\d+)','tokens');     
                        end
                        session = regexpOut{end}{1};
                        sessionTimes = regexp(session, '(\d+)','tokens');
                        sessionTimes = cat(2,sessionTimes{:});
                        formatTimes = cellfun(@(x) datevec([day,x],[obj.dateParser,obj.timeParser]),sessionTimes,'UniformOutput',0);
                        obj.dayStruct(i).sessions(k).date = cat(1,formatTimes{:});
                        obj.dayStruct(i).sessions(k).folderName = dayFolder{1};
                        obj.dayStruct(i).sessions(k).rootFolder = obj.rootFolder;
                        obj.dayStruct(i).sessions(k).tracesEventsFileName = tracesEventsFile{1};
                        
                        %Find also cellmap and decision files
                        analysisFileRegexp = cellfun(@(x) [session,'.*',x],obj.analysisFileRegexp,'UniformOutput',0);
                        analysisFileIdx = find(obj.checkRegexp({allFiles.name},analysisFileRegexp));
                        obj.dayStruct(i).sessions(k).analysisFileFileName = allFiles(analysisFileIdx).name;
                        
                        decisionsFileRegexp = cellfun(@(x) [session,'.*',x],obj.decisionsRegexp,'UniformOutput',0);
                        decisionsFileIdx = find(obj.checkRegexp({allFiles.name},decisionsFileRegexp));
                        obj.dayStruct(i).sessions(k).decisionsFileName = allFiles(decisionsFileIdx).name;                        
                        %
                        
                        k = k + 1;
                    end
                    i = i+1; %If no concat files were found, ignore the folder
                end
            end 
        end
        
    end
    
end