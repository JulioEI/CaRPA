classdef PFAnalysis < singleAnimalAnalysis

    properties
    end
    
    properties (Constant)
    end
    
    methods
        
        function obj = PFAnalysis(folder)
            if nargin == 0
                folder = [];
            end
            obj@singleAnimalAnalysis(folder)
        end
                
        function turnoverAlignmentMat = quantifyTurnoverAlignment(obj,varargin)
            %OPTIONAL INPUTS: filterPF, sessions, rField, binN, minVel, pxPerCm, dtCamera, dT, 
            
            %Get parameters
            param = obj.getDefaultParams({'binN','minVel','pxPerCm','dtCamera','dT','filterPF'},varargin);
                        
            %Get the sessions
            sessions = obj.chooseSessions(PFAnalysis.parseInput({varargin,'sessions',[]}));
            
            %Remove cells with non significant place fields
            if length(param.filterPF) > 1
                PCdecisions = param.filterPF;
            elseif param.filterPF
                [pos,r] = obj.getXR(varargin{:});
                [~,~,PCdecisions] = PFAnalysis.isPF(r,pos,param,varargin{:});                
            else
                PCdecisions = [];
            end
            
            %Get alignment for sessions
            globalIDs = PFAnalysis.alignSessions(sessions,PCdecisions);
            
            
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
            sessions = obj.chooseSessions(PFAnalysis.parseInput({varargin,'sessions',[]}));
            
            %Get traces and positions
            [pos,r] = obj.getXR(varargin{:});
            
            %Remove cells with non significant place fields
            if length(param.filterPF) > 1
                PCdecisions = param.filterPF;
            elseif param.filterPF
                [r,pos,PCdecisions] = PFAnalysis.isPF(r,pos,param,varargin);               
            else
                PCdecisions = cellfun(@(x) ones([1,size(x,2)]),r,'UniformOutput',false);
            end
            
            %Get alignment for sessions
            globalIDs = PFAnalysis.alignSessions(sessions,PCdecisions);
            
            for k = 1:length(sessions)
                
                halfT = round(length(r{k})/2);
                %Compute PF for first half of the session
                PF.firstHalf = PFAnalysis.buildPF(r{k}(1:halfT,:), pos{k}(1:halfT,:), param);
                %Compute PF for second half of the session
                PF.secondHalf = PFAnalysis.buildPF(r{k}(halfT+1:end,:), pos{k}(halfT+1:end,:), param); 
                       
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
            sessions = obj.chooseSessions(PFAnalysis.parseInput({varargin,'sessions',[]}));
            
            %Get traces and positions
            [pos,r] = obj.getXR(varargin{:});
            
            %Remove cells with non significant place fields
            if param.filterPF
                [r,pos,PCdecisions] = PFAnalysis.isPF(r,pos,param,varargin);
            else
                PCdecisions = cellfun(@(x) ones([1,size(x,2)]),r,'UniformOutput',false);
            end
            
            %Get alignment for sessions
            globalIDs = PFAnalysis.alignSessions(sessions,PCdecisions);

            %Build all PF
            allPF = arrayfun(@(i) PFAnalysis.buildPF(r{i}, pos{i}, param),1:length(sessions),'UniformOutput',0);
            
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
            [pos,r,nSessions,sessions] = obj.getXR(varargin);
            [~,v] = obj.getXR('sessions',sessions,'rField','velocity');     
            %Compute PF
            PF = cell([1,nSessions]);
            for k = 1:nSessions
                param.v = v{k};
                [PF{k}.mean,PF{k}.std,PF{k}.ste] = PFAnalysis.buildPF(r{k}, pos{k}, param);    
            end
            if nSessions == 1;PF = PF{1};end
            
        end
                
        function PFvar = computeVarianceExplainedUponShiftSession(obj,varargin)
            
        %OPTIONAL INPUTS: session, rField, bigN, smallN, smallShift, bigShift, smallBigDistanceShift, binN, minVel, pxPerCm, dtCamera, dT, 
            
            %Get parameters
            param = obj.getDefaultParams({'binN','minVel','pxPerCm','dtCamera','dT'},varargin);
            
            param.smallN = PFAnalysis.parseInput({varargin,'smallN',10});
            param.bigN = PFAnalysis.parseInput({varargin,'bigN',50});
            param.smallShift = PFAnalysis.parseInput({varargin,'smallShift',param.dtCamera}); %In seconds
            param.bigShift = PFAnalysis.parseInput({varargin,'bigShift',param.dtCamera*10}); %In seconds
            param.smallBigDistanceShift = PFAnalysis.parseInput({varargin,'smallBigDistanceShift',param.dtCamera*100}); %In seconds
            
            %Get responses
            [pos,r,nSessions] = obj.getXR(varargin);
            
            %Compute cells information
            PFvar = cell([1,nSessions]);
            for k = 1:nSessions
                PFvar{k} = PFAnalysis.computeVarianceExplainedUponShift(r{k}, pos{k}, param);    
            end
            if nSessions == 1;PFvar = PFvar{1};end
        end
        
        function PFinf = computeCellInformationSession(obj,varargin)
            %OPTIONAL INPUTS: session, rField, bigN, smallN, smallShift, bigShift, smallBigDistanceShift, binN, minVel, pxPerCm, dtCamera, dT, 
            
            %Get parameters
            param = obj.getDefaultParams({'binN','minVel','pxPerCm','dtCamera','dT'},varargin);
            
            param.smallN = PFAnalysis.parseInput({varargin,'smallN',10});
            param.bigN = PFAnalysis.parseInput({varargin,'bigN',50});
            param.smallShift = PFAnalysis.parseInput({varargin,'smallShift',param.dtCamera}); %In seconds
            param.bigShift = PFAnalysis.parseInput({varargin,'bigShift',param.dtCamera*10}); %In seconds
            param.smallBigDistanceShift = PFAnalysis.parseInput({varargin,'smallBigDistanceShift',param.dtCamera*100}); %In seconds
            
            %Get responses
            [pos,r,nSessions] = obj.getXR(varargin);
            
            %Compute cells information
            PFinf = cell([1,nSessions]);
            for k = 1:nSessions
                PFinf{k} = PFAnalysis.computeFullCellinformation(r{k}, pos{k}, param);    
            end
            if nSessions == 1;PFinf = PFinf{1};end
            
        end
                        
    end
    
    methods(Static)
        
        function compareSpikeOutputs(output,varargin)
            
            traces = output.(PFAnalysis.parseInput({varargin,'traces','rawProb'}));
            spikes = PFAnalysis.parseInput({varargin,'spikes',{'spikeDeconv','spikeML','tresholdEvents'}});
            
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
            
            smallN = PFAnalysis.parseInput({varargin,'smallN',10});
            bigN = PFAnalysis.parseInput({varargin,'bigN',20});
            smallShift = PFAnalysis.parseInput({varargin,'smallShift',param.dtCamera}); %In seconds
            bigShift = PFAnalysis.parseInput({varargin,'bigShift',param.dtCamera*10}); %In seconds
            smallBigDistanceShift = PFAnalysis.parseInput({varargin,'smallBigDistanceShift',param.dtCamera*100}); %In seconds

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
                    PF = PFAnalysis.buildPF(rShift, x{k}, param);
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
        
        function varOut = computeVarianceExplainedUponShift(r,x,param)
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
            varExplained = [];
            textprogressbar('Computing var explained... ');
            i = 0;
            PFLog = [];
            for shift = [shiftSlowF,shiftFastF]
                textprogressbar(100*i/length([shiftSlowF,shiftFastF]));i = i + 1;
                %Shift responses
                rShift = circshift(r,shift,1);
                 
                %Compute PF
                [PF,PFstd] = PFAnalysis.buildPF(rShift, x, param);
                
                %Save pf 
                PFLog = cat(4,PFLog,PF);
                
                %Compute variance explained
                varOfTheMean = nanstd(reshape(PF,[size(PF,1)*size(PF,2),size(PF,3)])).^2;
                meanOfTheVars = squeeze(nanmean(nanmean(PFstd.^2)))';                
                varExp = 1 - meanOfTheVars./varOfTheMean; 
                
                varExplained = [varExplained;varExp];
            end
            textprogressbar(' done');
            varOut.varExplained = varExplained;
            varOut.PFLog = PFLog;
            [pFmean,~,pFerror] = PFAnalysis.buildPF(r, x, param);
            varOut.zeroShiftPF.mean = pFmean;
            varOut.zeroShiftPF.ste = pFerror;
            varOut.shiftSlowSec = shiftSlowSec;
            varOut.shiftFastSec = shiftFastSec;
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
                PF = PFAnalysis.buildPF(rShift, x, param);
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
            [pFmean,~,pFerror] = PFAnalysis.buildPF(r, x, param);
            info.zeroShiftPF.mean = pFmean;
            info.zeroShiftPF.ste = pFerror;
            info.shiftSlowSec = shiftSlowSec;
            info.shiftFastSec = shiftFastSec;
        end
        
        function viewPFvar(PFvar)
            tempFields = fields(PFvar);
            if ~isstruct(PFvar.(tempFields{1}))
                tempInfo.response = PFvar;
                PFvar = tempInfo;
            end
            myFields = fields(PFvar);

            cellN = size(PFvar.(myFields{1}).zeroShiftPF.mean,3); 
            fig = figure;
            fig2 = figure;
            for celli = 1:cellN
                clf;
                for k = 1:length(myFields)
                    figure(fig)
                    pFmean = PFvar.(myFields{k}).zeroShiftPF.mean;
                    pfVarCell = PFvar.(myFields{k}).varExplained';
                    shiftSlowSec = PFvar.(myFields{k}).shiftSlowSec;
                    shiftFastSec = PFvar.(myFields{k}).shiftFastSec;
                    pfLog = PFvar.(myFields{k}).PFLog;
                    [validVar, pValVar] = ttest2(pfVarCell(celli,1:length(shiftSlowSec))',pfVarCell(celli,length(shiftSlowSec)+1:end)',0.05);
                    if validVar; pColorVar = 'green';else;pColorVar = 'red';end
                    subplot(length(myFields),2,1+(k-1)*2)
                    plot([shiftSlowSec,shiftFastSec],pfVarCell(celli,:),'-o','color',pColorVar);xlabel('delay(s)');ylabel('Var Explained');title(['pvalue ',num2str(pValVar)]);axis square;grid on;ylim([-inf,0])
                    subplot(length(myFields),2,2+(k-1)*2)
                    %shadedErrorBar(1:size(pFmean,1),pFmean(:,1,celli),pFerror(:,1,celli),'lineprops',{'-o','linewidth',2,'color',pColorInf},'patchSaturation',0.2);title(['cell',num2str(celli)]);axis square
                    %imagesc(imgaussfilt(pFmean(:,:,celli),0.5));title(['cell',num2str(celli)]);axis square;colorbar;colormap('jet')
                    imagesc(pFmean(:,:,celli));title(['cell',num2str(celli)]);axis square;colorbar;%colormap('jet')
                    if length(myFields)>1;annotation('textbox',[0,0,1,(1+1/length(myFields))-k*(1-1/length(myFields))/(length(myFields)-1)],'string',myFields{k},'fontSize',20,'FontWeight','bold');end
                    figure(fig2)
                    
                    for shft = 1:size(pfLog,4)
                        if shft < length(shiftSlowSec)
                            subplot(ceil(sqrt(size(pfLog,4))),ceil(sqrt(size(pfLog,4))),shft)
                        else
                            subplot(ceil(sqrt(size(pfLog,4))),ceil(sqrt(size(pfLog,4))),shft+1)
                        end
                        imagesc(pfLog(:,:,celli,shft));axis off;
                    end
                end
                pause;
            end
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
                    
%                     if binY == 1
                        miCellx = squeeze(info.(myFields{k}).miCell(:,1,:));
                        [validInf, pValInf] = ttest2(pfEntropyCell(celli,1:length(shiftSlowSec))',pfEntropyCell(celli,length(shiftSlowSec)+1:end)',0.05,'left');
                        [validMI, pValMI] = ttest2(miCellx(celli,1:length(shiftSlowSec))',miCellx(celli,length(shiftSlowSec)+1:end)',0.05,'right');
                        if validInf; pColorInf = 'green';else;pColorInf = 'red';end
                        if validMI; pColorMI = 'green';else;pColorMI = 'red';end
                        subplot(length(myFields),3,1+(k-1)*3)
                        plot([shiftSlowSec,shiftFastSec],pfEntropyCell(celli,:),'-o','color',pColorInf);xlabel('delay(s)');ylabel('entropy');title(['pvalue ',num2str(pValInf)]);axis square;grid on;ylim([0,inf])
                        subplot(length(myFields),3,2+(k-1)*3)
                        %shadedErrorBar(1:size(pFmean,1),pFmean(:,1,celli),pFerror(:,1,celli),'lineprops',{'-o','linewidth',2,'color',pColorInf},'patchSaturation',0.2);title(['cell',num2str(celli)]);axis square
                        imagesc(imgaussfilt(pFmean(:,:,celli),1));title(['cell',num2str(celli)]);axis square;colorbar;colormap('jet')
                        subplot(length(myFields),3,3+(k-1)*3)
                        plot([shiftSlowSec,shiftFastSec],miCellx(celli,:),'-o','color',pColorMI);xlabel('delay(s)');ylabel('mutual info');title(['pvalue ',num2str(pValMI)]);axis square;grid on;ylim([0,inf])
                        if length(myFields)>1;annotation('textbox',[0,0,1,(1+1/length(myFields))-k*(1-1/length(myFields))/(length(myFields)-1)],'string',myFields{k},'fontSize',20,'FontWeight','bold');end
%                     else
%                        error('Not yet implemented')
%                     end
                end
                pause;
            end
        end
        
        function [PFmeanRate,PFstdRate,PFsteRate,binCounts] = buildPF(r, x, param)
%             disp('Computing PF...')
%             disp(repmat('.',[1,100]))
            [xTB,rTB,binCounts] = dataAnalysis.curateXandR(x,r,param);
            
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
            PFmeanRate = PFmeanCount;%./dT;
            PFstdRate = PFstdCount;%./dT;
            PFsteRate = PFsteCount;%./dT;
            
           %Constraint PF inside circle if nescessary
            if isfield(param,'boxRadi')
                %%Replace nan with 0
                PFmeanRate(isnan(PFmeanRate)) = 0;
                PFstdRate(isnan(PFstdRate)) = 0;
                PFsteRate(isnan(PFsteRate)) = 0;
                nanMask = PFAnalysis.maskPF(param.boxRadi,param.boxDims,x,xTB,param);
                PFmeanRate = PFmeanRate + nanMask;
                PFstdRate = PFstdRate + nanMask;
                PFsteRate = PFsteRate + nanMask;
            end

        end
        
        function nanMask = maskPF(boxRadi,boxDims,x,xTB,param,doPlot)
            if nargin < 6
                doPlot = 0;
            end
            xCenter = boxDims(1);
            yCenter = boxDims(2);
            theta = 0 : 0.01 : 2*pi;
            radius = boxRadi;
            xP = radius * cos(theta) + xCenter;
            yP = radius * sin(theta) + yCenter;
            param.minVel = 0;
            xPTB = dataAnalysis.curateXandR([xP',yP'],[],param);
            nanMask = nan(param.binN);
            for i = 1:size(nanMask,1)
                for j = 1:size(nanMask,2)
                    if inpolygon(i,j,xPTB(:,1),xPTB(:,2))
                        nanMask(i,j) = 0;
                    end
                end
            end

            if doPlot
                figure;
                subplot(1,3,1)
                plot(x(:,1),x(:,2));
                hold on;
                plot(xP,yP,'LineWidth',3,'Color','r')
                %rectangle('Position',[boxDims(1)-boxRadi boxDims(2)-boxRadi 2*boxRadi 2*boxRadi],'Curvature',[1 1],'LineWidth',3,'EdgeColor','r')
                axis equal;  
                subplot(1,3,2)
                plot(xTB(:,1),xTB(:,2));
                hold on
                plot(xPTB(:,1),xPTB(:,2),'LineWidth',3,'Color','r');
                axis equal
                subplot(1,3,3)
                imagesc(nanMask);
                axis equal
            end
%             
%             for celli = 1:size(PF,3)
%                 for i = 1:size(PF,1)
%                     for j = 1:size(PF,2)
%                         if isnan(nanMask(i,j))
%                             PF{i,j,celli} = nan;
%                         end
%                     end
%                 end
%             end
            
        end
        
        function viewPF(PF, varargin)
            
            normalize = PFAnalysis.parseInput({varargin,'normalize',0});
            
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
        
        function sessions = sortSessionsByTime(sessions)
           dates_mat = cell2mat(cellfun(@(x) x(1,:)', {sessions.date},'UniformOutput',0))';
           dates = datetime(dates_mat,'InputFormat','yyyy MM dd HH mm');
           [~,idx] = sort(dates);
           sessions = sessions(idx);
        end
        
        function viewGlobalIDs(sessions,globalIDs)
            sessions = PFAnalysis.sortSessionsByTime(sessions);
            labels = {sessions.folderName};
            %labels = cellfun(@(x) ['\begin{tabular}{c}', strrep(x,'-',' \\ '), '\end{tabular}'], labels, 'UniformOutput',0);
            figure;
            img = imagesc(globalIDs);
            colormap(jet)
            set(img,'AlphaData',globalIDs~=0)
            xticks(1:length(sessions))
            set(gca, 'XTickLabel', labels,'XTickLabelRotation',90);%, 'TickLabelInterpreter', 'latex'
            colorbar;
            ylabel('Global ID')
        end
        
        function globalIDs = alignSessions(sessions, PCdecisions)

            if nargin < 2
                PCdecisions = [];
            end
            
           %Sort sessions by time
           sessions = PFAnalysis.sortSessionsByTime(sessions);
           
           %Load session cellmaps
           allCellMaps = {};
           for k = 1:length(sessions)
                decisions = load(fullfile(sessions(k).rootFolder,sessions(k).folderName,sessions(k).decisionsFileName));
                emFile = load(fullfile(sessions(k).rootFolder,sessions(k).folderName,sessions(k).analysisFileFileName));
                allCellMaps{k} = emFile.emAnalysisOutput.cellImages(:,:,logical(decisions.validCellMax));%Take good cells
%                 allCellMaps{k} = emFile.emAnalysisOutput.cellImages; %Take all cells
           end
           disp(' O ');
           allCellMaps = cellfun(@(x) thresholdImages(x),allCellMaps,'UniformOutput',0);%Treshold cells
           
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
%                disp([' %%%%%%%%%%%%% Error is here ? ==> j=',num2str(j)])
                %figure;subplot(1,2,1);imagesc(allCellMapsMax{j});subplot(1,2,2);imagesc(template)%See
                %the alighments
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
           
           %Look at results
           %viewAlignmentResults(allCellMaps,alignedImgNC,globalIDs)
           
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
                
        function imagescNanAsWhite(X)
            figure;imagesc(X);colorbar;title('AutoCorr')
            alphaMap = (X == 0);
            white = cat(3, ones(size(X)), ones(size(X)), ones(size(X))); 
            hold on; h = imagesc(white); hold off; set(h, 'AlphaData', alphaMap);
        end
        
    end
    
    
end