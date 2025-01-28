classdef cellInfo < handle 
    properties 
        %Global
        timeAverageEvents = true;
        timeWindowSize = 10;
        timeWeightDecay = 2;
        smoothImages = true;
        cropSizeLength = 10;
        filterSpaceRatio = 0.1;
        filterHeight = 5;
        timeSeqGlobalSNR = -3:3;
        spikeROI = -20:20;
        slopeFrameWindow = 10;
        maxShowEvents = 50;
        calHeight = 3;
        eventStdDetection = 2;
        validCells
    end
    
    properties(SetAccess = private)
        %Global
        inputMovie
        inputMovieStr
        inputImages
        inputSignals
        signalPeakIdx
        cellN
        minValMovie
        maxValMovie 
        inputImagesT
        %Per cell
        croppedPeakImages
        inputImage
        inputSignal
        signalPeakId
        cellImg
        PSF
        cellImgF
        cellMask
        cellLine
        i
        %Per cell (events)
        eventMask
        cellEvents
        eventPeaks
        cellEventsF
        cellEventsTimeAvg
        meanValVect
        tresholdedImg
        tresholdedOutside
    end
    
    properties (Access = private)
         optimization;
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%CONSTRUCTOR%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        function obj = cellInfo(inputMovie,inputImages,inputSignals,signalPeakIdx,varargin)
            if nargin < 1
               [movieName,pathMovie] = uigetfile('*.h5','Enter downsampled calcium movie');
               inputMovie = [pathMovie,filesep,movieName];
               [cellsFile,pathCellsFile] = uigetfile([pathMovie,filesep,'*.mat'],'Enter the found cells file');
               load([pathCellsFile,filesep,cellsFile]);
               try 
                   inputImages = permute(emAnalysisOutput.cellImages,[3,1,2]);
                   inputSignals = emAnalysisOutput.scaledProbability;
                   signalPeakIdx = emAnalysisOutput.eventTimes;
               catch
                   inputImages = permute(pcaicaAnalysisOutput.IcaFilters,[3,1,2]);
                   inputSignals = pcaicaAnalysisOutput.IcaTraces;
                   signalPeakIdx = [];
               end
               answer = questdlg('Enter decisions?', ...
                   'Do you want to enter a decisions file?', ...
                   'Yes','No','No');
               if strcmp(answer,'Yes')
                   [inputDecisions,pathDecisions] = uigetfile([pathCellsFile,filesep,'*.mat'],'Enter the decisions cells file');
                   decisions = load([pathDecisions,filesep,inputDecisions]);dummyField = fields(decisions);decisions = decisions.(dummyField{1});
                   obj.validCells = decisions;
               end
               if ~isempty(signalPeakIdx)
               answer = questdlg('Recalculate peaks?', ...
                   'Do you want to recompute the signal peaks', ...
                   'Yes','No','Yes');
               if strcmp(answer,'Yes');signalPeakIdx = [];end    
               end        
            end
            %Get arguments
            if isstr(inputMovie)
                obj.inputMovieStr = inputMovie;
                disp('Loading movie...')
                hinf = hdf5info(inputMovie);
                inputMovie = hdf5read(hinf.GroupHierarchy.Datasets);
            end
            if size(inputImages,1) == size(inputMovie,1)%Images need to be in nImageXwidthXheight
                inputImages = permute(inputImages,[3,1,2]);
            end
            obj.inputMovie = inputMovie;
            if size(inputImages,2) ~= size(inputImages,3) && length(unique(size(inputImages)))~=3
                warning('First dimension of inputImages must be number of cells')
            end
            obj.inputImages = inputImages;
            obj.inputSignals = double(inputSignals);
            if isempty(signalPeakIdx)
                robustInputSignals = obj.inputSignals;
                robustInputSignals(isnan(robustInputSignals)) = 0;
                [~, signalPeakIdx] = computeSignalPeaks(robustInputSignals,'numStdsForThresh',obj.eventStdDetection);
                
                %Force event finding in no event cells, take max (temporary
                %fix)
                noEventCells = find(cellfun(@isempty,signalPeakIdx));
                for noEventCell = noEventCells
                    [~,peakIdx] = max(robustInputSignals(noEventCell,:));
                    signalPeakIdx{noEventCell} = peakIdx;
                end          
            end
            obj.signalPeakIdx = signalPeakIdx;
            
            if size(inputSignals,2) ~= size(inputMovie,3)
                error('Movie and signals have different lengths')
            end
            if size(inputSignals,1) ~= size(inputImages,1)
                error('Images and signals have different cell count')
            end
            %Set defaults
            obj.cellN = size(inputImages,1);
            obj.minValMovie = min(inputMovie(:));
            obj.maxValMovie = max(inputMovie(:));
            if isempty(obj.validCells);obj.validCells = 3*ones([1,obj.cellN]);end

%             randVect = randperm(obj.cellN);
%             obj.validCells(randVect(1:floor(obj.cellN/3)))=0;
%             obj.validCells(randVect(ceil(obj.cellN/3):2*floor(obj.cellN/3)))=1;
            
            %Initialize optimitztion obj
            obj.optimization = obj.createOptimizationObj;
            
            %Parse varargin arguments, as ('name',value)
            for k = 1:2:length(varargin)
                obj.(varargin{k}) = varargin{k+1};
            end  
            obj.inputImagesT = permute(thresholdImages(permute(inputImages,[2,3,1])),[3,1,2]);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%PUBLIC METHODS%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%  
        function gui(obj)
            stay = 1;
            fn = {'Manual sort','Run decoder','Save decisions','Set remaining cells to bad','Reset all decisions','Quit menu'};
            while stay
            [indx,res] = listdlg('PromptString','Select what to do:',...
                       'SelectionMode','single',...
                       'ListString',fn); 
            if ~res 
                indx = length(fn);
            end
               switch indx
                   case 1
                       obj.manualSort('tresholds',{'getOverlap','>=0.5','getCellShapePropieties.Eccentricity','<=0.9','getCellShapePropieties.Solidity','>0.8'})%,{'getOverlap','>=.4','getglobalSNR','>1'});
                   case 2
                       obj.predictClassifier;
                   case 3
                       answer = questdlg('Enter method', ...
                       'Please enter the method used to extract the cells:', ...
                       'em','ica','ica');
                       obj.saveDecisions(answer)
                   case 4
                       answer = questdlg('This will set all ambiguous cells to bad and cannot be undone, proceed?', ...
                       'Confirmation?', ...
                       'yes','no','no');
                        if strcmp(answer,'yes')
                            obj.validCells(obj.validCells == 3) = 0;
                        end  
                   case 5
                       answer = questdlg('This will forever erase all the current decisions, proceed?', ...
                       'Confirmation?', ...
                       'yes','no','no');  
                        if strcmp(answer,'yes')
                            obj.validCells = 3*ones(size(obj.validCells));
                        end  
                   case length(fn)
                       stay = 0;
               end
            end
        end
        
        function saveDecisions(obj,type,varargin)
            avaliableType = {'em','ica'};
            if sum(strcmp(type,avaliableType)) == 0
                disp(cell2table(avaliableType))
                error(['Type "',type,'" not reconized.'])
            end
            
            switch type
                case 'ica'
                    valid = obj.validCells;
                    variableName = 'valid';
                    extensionName = 'ICdecisions.mat';
                case 'em'
                    validCellMax = obj.validCells;
                    variableName = 'validCellMax';
                    extensionName = 'emAnalysisSorted.mat';
                otherwise
                    error(['Sorry, type: ',type,' not yet implemented'])
            end
            
            decisionsName = [];
            if ~isempty(obj.inputMovieStr)%if movie name exists, parse input from it.
                splitArray = split(obj.inputMovieStr,'\');
                if length(splitArray{end}) > 1
                    fileName = splitArray{end};
                    pathNameCell = splitArray(1:end-1);
                else
                    fileName = splitArray{end-1}; %if the archive ends in '\'
                    pathNameCell = splitArray(1:end-2);
                end
                splitArray = split(fileName,{'turboreg','crop','dfof','downsasmple'});
                if length(splitArray) > 1
                    decisionsName = splitArray{1};
                end
                pathName = '';
                for folderName = pathNameCell'
                    pathName = [pathName,folderName{1},'\'];
                end
            end
            if ~isempty(decisionsName)
                if ~obj.parseInput({varargin,'skipConfirmation',0})
                    choice = questdlg(['The file: ',newline,newline,char(strcat(pathName,decisionsName,extensionName)),newline,newline,' is going to be created. Is this correct?'], ...
                    'Attention', ...
                    'Yes','Let me change the folder','Let me change the folder and name','Yes');
                % Handle response
                    switch choice
                        case 'Yes'
                        case 'Let me change the folder'
                            pathName = [uigetdir,'\'];
                        case 'Let me change the folder and name'
                            decisionsName = [];
                        otherwise
                            disp('Aborting save')
                            return;
                    end
                end
            end
            if isempty(decisionsName) %manually assemble the name otherwise.
                manualDate = inputdlg('Enter date of file (yyyy_mm_dd)');
                manualProtocol = inputdlg('Enter protocol of file (leave empty for default)');
                manualMouse = inputdlg('Enter animal number of file');
                manualSession = inputdlg('Enter session of file (leave empty for default)');
                if isempty(manualDate{1})
                    manualDate{1} = '0000_00_00';
                end
                if isempty(manualProtocol{1})
                    manualProtocol{1} = 'p000';
                end
                if isempty(manualMouse{1})
                    manualMouse{1} = '0000';
                end
                if isempty(manualSession{1})
                    manualSession{1} = 'NULL000';
                end
                decisionsName = [manualDate{1},'_',manualProtocol{1},'_mouse',manualMouse{1},'_',manualSession{1},'_'];
                pathName = [uigetdir,'\'];
            end
            save(strcat(pathName,decisionsName,extensionName),variableName)
        end
        
        function [sScore] = getsScore(obj,i,varargin)    
            createCellVariables(obj,i);
            
%             K = 0.125*ones(3);    
            K = 1*ones(1);
            sIMG = conv2(obj.cellImg,K,'same');
            
            sScore = sum(sum(imregionalmax(sIMG)))/numel(sIMG);
%             sScore = sum(sum(abs(del2(obj.cellImg))))/7.36; %7.36 is mean noise value
            if obj.parseInput({varargin,'doPlot',0})
%                 if obj.parseInput({varargin,'openFig',1});figure;end
%                 subplot(1,2,1)
%                 meshc(obj.cellImg)
%                 title('Cell shape')
                subplot(1,2,1)
                meshc(sIMG')
                xlabel(sScore)
                [tempX,tempY] = find(imregionalmax(sIMG));
                hold on;
                for i = 1:length(tempX)
                    plot3(tempX(i),tempY(i),sIMG(tempX(i),tempY(i)),'.','markerSize',50)
                end
                hold off;
                subplot(1,2,2)
                imagesc(obj.cellImg);
%                 subplot(1,2,2)
%                 meshc(abs(del2(obj.cellImg)))
%                 title('Absolute value of Laplacian')
            end
        end
        
        function [centroidDist] = getcentroidDist(obj,i,varargin)
            createCellVariables(obj,i);
            cellMap = squeeze(obj.inputImage);
            [x,y] = find(cellMap==max(cellMap(:)));
            if length(x)>1 || length(y)>1; x=x(1); y=y(1); disp('@@@@ ATTENTION This is a manual correction (Pasha)');end
            cellMapCenter = size(cellMap)/2;
            centroidDist = sqrt((x-cellMapCenter(1))^2 + (y-cellMapCenter(2))^2)/sqrt((size(cellMap,1)-cellMapCenter(1))^2 + (size(cellMap,2)-cellMapCenter(2))^2);
            if obj.parseInput({varargin,'doPlot',0})
                if obj.parseInput({varargin,'openFig',1});figure;end
                imagesc(cellMap')
                hold on;
                line([x,cellMapCenter(1)],[y,cellMapCenter(2)],'Color',[.8 .8 .8])
                plot(x,y,'rx')
                plot(cellMapCenter(1),cellMapCenter(2),'kx')
            end
        end
        
        function [globalSNR] = getglobalSNR(obj,i,varargin)
            createCellVariables(obj,i);
            createEventVariables(obj,varargin);
            peakIdx = obj.computePeakIdx;
            noiseSignal = obj.inputSignal(setdiff(1:length(obj.inputSignal),peakIdx));
            noiseSignal(noiseSignal > median(obj.inputSignal)+2*std(obj.inputSignal)) = [];
%             globalSNR = nanmean(obj.inputSignal(obj.eventPeaks)/nanstd(noiseSignal));
%             globalSNR = nanmean(obj.inputSignal(obj.eventPeaks))/nanmean(noiseSignal);
            meanPeakVal = nanmean(abs(obj.inputSignal(obj.eventPeaks)));
            percentilVal = quantile(noiseSignal,0.95);
            globalSNR = (meanPeakVal-percentilVal)/percentilVal;
%             globalSNR = nanmean(abs(obj.inputSignal(obj.eventPeaks)))/abs(nanmean(obj.inputSignal(:)));
            if obj.parseInput({varargin,'doPlot',0})
                if obj.parseInput({varargin,'openFig',1});figure;end
                stem(obj.inputSignal,'k','marker','none')
                hold on
                stem(peakIdx,abs(obj.inputSignal(peakIdx)),'.','marker','none')
                plot(obj.eventPeaks,obj.inputSignal(obj.eventPeaks),'.','markersize',10)
                xlabel('idx')
                ylabel('Signal')
                legend('Noise','WindowedPeaks','Peaks')
            end
        end
        
        function [RMS] = getRMS(obj,i,varargin)
            createCellVariables(obj,i);
            createEventVariables(obj,varargin);
            peakIdx = obj.computePeakIdx;
            RMS = mean(abs(obj.inputSignal(peakIdx)));
            if obj.parseInput({varargin,'doPlot',0})
                if obj.parseInput({varargin,'openFig',1});figure;end
                stem(obj.inputSignal,'k','marker','none')
                hold on
                stem(peakIdx,abs(obj.inputSignal(peakIdx)),'.','marker','none')
                xlabel('idx')
                ylabel('Signal')
                plot(obj.eventPeaks,obj.inputSignal(obj.eventPeaks),'.','markersize',10)
                legend('Noise','WindowedPeaks','Peaks')
            end
        end
        
        function [peakSNR] = getpeakSNR(obj,i,varargin)
            createCellVariables(obj,i);
            createEventVariables(obj,varargin);
            spikeCenterTrace = obj.computeSpikeCenterTrace;
            peakSNR = std(mean(spikeCenterTrace,1))/mean(std(spikeCenterTrace));
            if obj.parseInput({varargin,'doPlot',0})
                if obj.parseInput({varargin,'openFig',1});figure;end
                plot(spikeCenterTrace','color',[.5,.5,.5]);
                hold on
                shadedErrorBar(1:length(mean(spikeCenterTrace,1)),mean(spikeCenterTrace,1),std(spikeCenterTrace))
                plot(mean(spikeCenterTrace,1),'r--','linewidth',3)
                title(['peakSNR ',num2str(peakSNR)])
            end
        end
        
        function expRatio = getexpRatio(obj,i,varargin)
            createCellVariables(obj,i);
            createEventVariables(obj,varargin);
            spikeCenterTrace = obj.computeSpikeCenterTrace;
            prePeakIdx = find(obj.spikeROI==-(obj.slopeFrameWindow)):find(obj.spikeROI==0);
            postPeaklIdx = find(obj.spikeROI==0):find(obj.spikeROI==obj.slopeFrameWindow);
            
%             expRatio = zeros([1,length(obj.eventPeaks)]);
            if obj.parseInput({varargin,'doPlot',0})
                if obj.parseInput({varargin,'openFig',1});figure;end
            end
%             for k = 1:length(expRatio)
                
%                 meanSignalPre =  spikeCenterTrace(k,prePeakIdx);
%                 meanSignalPost =  spikeCenterTrace(k,postPeaklIdx);
                meanSignalPre =  mean(spikeCenterTrace(:,prePeakIdx),1);
                meanSignalPost =  mean(spikeCenterTrace(:,postPeaklIdx),1);  
                
                f = {NaN,NaN};
                ind = 0;
                for meanSignal = [flip(meanSignalPre)',meanSignalPost']
                    ind = ind + 1;
%                     [~,locs] = findpeaks(diff(meanSignal));
%                     if isempty(locs) %meanSignal has no minima
%                         [~,locs] = max(diff(meanSignal));
%                         locs = flip(locs);
%                     end
                    trimMeanSignal = meanSignal;%(1:locs(1));
                    x = 1:size(trimMeanSignal);
                    
                    w = ((x+length(x)/2).^2)/max((x+length(x)/2).^2);
                    coeffInit = polyfit(x',log(trimMeanSignal),1);
                    if sum(isinf(coeffInit)) || sum(~isreal(coeffInit))
                        coeffInit = [-1,log(1)];
                    end
                    warning('off')
                    try
                        coeff = fitnlm(x',trimMeanSignal,@(coeff,x) coeff(1).*exp(coeff(2)*x),[exp(coeffInit(2)),coeffInit(1)],'Weight',w);
                        f{ind} = struct('a' ,coeff.Coefficients.Estimate(1), 'b',coeff.Coefficients.Estimate(2));
                    catch
                        warning('Could not fit trace')
                        f{ind} = struct('a' , NaN, 'b', NaN);
                    end
                    warning('on')
%                     try
%                         warning('off')
%                         f{ind} = fit(x',trimMeanSignal,'exp1');
%                         warning('on')
%                     catch
%                         try
%                             coeff = polyfit(x',log(trimMeanSignal),1);
%                             f{ind} = struct('a' , exp(coeff(2)), 'b', coeff(1));
%                         catch
%                             warning('Could not fit the trace')
%                             f{ind} = struct('a' , NaN, 'b', NaN);
%                         end
%                     end
                end
                expRatio = f{1}.b/f{2}.b;
                
                if obj.parseInput({varargin,'doPlot',0})
%                     subplot(ceil(sqrt(length(expRatio))),ceil(sqrt(length(expRatio))),k)            
                    fittedPre = f{1}.a*exp(f{1}.b*(1:0.01:size(meanSignalPre,2)));
                    fittedPost = f{2}.a*exp(f{2}.b*(1:0.01:size(meanSignalPost,2)));            
                    myStd = [std(spikeCenterTrace(:,prePeakIdx),1),std(spikeCenterTrace(:,postPeaklIdx),1)];
                    shadedErrorBar([prePeakIdx,postPeaklIdx],[meanSignalPre,meanSignalPost],myStd,{'b','linewidth',2})
                    hold on;
                    plot(linspace(prePeakIdx(1),prePeakIdx(end),size(fittedPre,2)),flip(fittedPre),'r','linewidth',2);
                    plot(linspace(postPeaklIdx(1),postPeaklIdx(end),size(fittedPost,2)),fittedPost,'r','linewidth',2);
%                     title(['expRatio: ',num2str(expRatio(k))])
                    title(['expRatio: ',num2str(expRatio)])
                    plot([prePeakIdx(end),prePeakIdx(end)],[0,max(meanSignalPre)],'b--','linewidth',2)
                    hold off;
                    axis square
                    grid on
                end
%             end
        end
        
        function [overlap] = getOverlap(obj,i,varargin)
            createCellVariables(obj,i);
            createEventVariables(obj,varargin);
            maxOverlapEvents = obj.maxShowEvents;
            eventsNum = size(obj.tresholdedImg,3);if eventsNum > maxOverlapEvents; eventsNum = maxOverlapEvents;end
            overlap = zeros([1,eventsNum]);
            if obj.parseInput({varargin,'doPlot',0})
                if obj.parseInput({varargin,'openFig',1});figure;end
                subplot(ceil(sqrt(1+size(obj.tresholdedImg,3))),ceil(sqrt(1+size(obj.tresholdedImg,3))),1)
                imagesc((obj.cellImgF>nanmean(obj.cellImgF(:))+1*nanstd(obj.cellImgF(:)))')
%                 imagesc(obj.cellMask)
            end
            outWeight = .2; %penalization for having the event outbound the mask [0,1]. 
            
            for k = 1:eventsNum
                missed = sum(sum(0~=(obj.cellImgF>nanmean(obj.cellImgF(:))+1*nanstd(obj.cellImgF(:))))) - sum(sum(0~=obj.tresholdedImg(:,:,k)));
                overlap(k) = sum(sum(0~=obj.tresholdedImg(:,:,k)))/(outWeight*missed+sum(sum(0~=obj.tresholdedOutside(:,:,k)))); 
                if obj.parseInput({varargin,'doPlot',0})
                    subplot(ceil(sqrt(1+eventsNum)),ceil(sqrt(1+eventsNum)),k+1)
                    imagesc(( obj.tresholdedOutside(:,:,k)*-min(obj.meanValVect(1:eventsNum)) + obj.tresholdedImg(:,:,k))')
                    caxis([-abs(min(obj.meanValVect(1:eventsNum))),max(obj.meanValVect(1:eventsNum))])
                    title(['Overlap ',num2str(overlap(k))])
                end
            end
            if obj.parseInput({varargin,'weighted',0})
                weights = (obj.meanValVect(1:eventsNum)-min(obj.meanValVect(1:eventsNum)))/sum(obj.meanValVect(1:eventsNum)-min(obj.meanValVect(1:eventsNum)));
                if obj.parseInput({varargin,'doPlot',0})
                    subplot(ceil(sqrt(1+size(obj.tresholdedImg,3))),ceil(sqrt(1+size(obj.tresholdedImg,3))),1)
                    title(['W. mean overlap: ',num2str(weights*shapingFun(overlap)')]);
                    figure;bar([overlap;shapingFun(overlap)]');legend('normal','sigmoid')
                end
                overlap = weights*shapingFun(overlap)';
            end
            function x = shapingFun(x)
                x = x+2*x.^4;
%                 x(x>.5) = 2./(1+exp(-10*(x(x>.5))+5));
            end
        end
        
        function regProps = getCellShapePropieties(obj,i,varargin)
            createCellVariables(obj,i);
            regProps = regionprops('table',obj.cellMask,'Area','Eccentricity','Solidity');
            [~,index] = max([regProps.Area]);
            regProps = regProps(index,:);
            
            if obj.parseInput({varargin,'doPlot',0})
                if obj.parseInput({varargin,'openFig',1});figure;end
                subplot(1,2,1)
                imagesc(obj.cellImg)
                title('Image')
                subplot(1,2,2)
                imagesc(obj.cellMask)
                title('eventMask')
            end
        end
        
        function regPropsAll = getCellEventPropierties(obj,i,varargin)
            createCellVariables(obj,i);
            createEventVariables(obj,varargin);
            if obj.parseInput({varargin,'doPlot',0});figure;end
            regPropsAll = cell([1,size(obj.tresholdedImg,3)]);
            for k = 1:size(obj.tresholdedImg,3)
                try
                    regProps = regionprops('table',logical(obj.tresholdedImg(:,:,k)),'Area','Eccentricity','Solidity'); 
                catch
                    [Area,Eccentricity,Solidity] = deal(NaN);
                    regProps = table(Area,Eccentricity,Solidity);
                end
               [~,index] = max([regProps.Area]);
               regProps = regProps(index,:);
               if obj.parseInput({varargin,'doPlot',0})
                   subplot(ceil(sqrt(size(obj.tresholdedImg,3))),ceil(sqrt(size(obj.tresholdedImg,3))),k)
                   imagesc(obj.tresholdedImg(:,:,k))
               end
               regPropsAll{k} = regProps;
            end
        end
        
        function interCorr = getinterCorr(obj,i,varargin)
            createCellVariables(obj,i);
            createEventVariables(obj,varargin);
            patchTemp = reshape(obj.cellEvents,[size(obj.cellEvents,1)*size(obj.cellEvents,2),size(obj.cellEvents,3)]);
            patchTemp(isnan(patchTemp)) = 0;
            corrMat = corrcoef(patchTemp);
            corrMatUpper = triu(corrMat,1);
            interCorr = mean(corrMatUpper(corrMatUpper~=0));
            if obj.parseInput({varargin,'doPlot',0})
                figure;
                imagesc(corrMat)
                title('Correlation matrix')
                figure;
                for k = 1:size(obj.cellEvents,3)
                    subplot(ceil(sqrt(size(obj.cellEvents,3))),ceil(sqrt(size(obj.cellEvents,3))),k)
                    imagesc(obj.cellEvents(:,:,k))
                    caxis([0,max(obj.cellEvents(:))])
                end
            end
        end
        
        function imMovCorr = getimMovCorr(obj,i,varargin)
            createCellVariables(obj,i);
            createEventVariables(obj,varargin);
            imMovCorr = ones([1,size(obj.cellEvents,3)]);
            cellImgNormalized = -1 + 2.*(obj.cellImg - min(min(obj.cellImg)))./(max(max((obj.cellImg ))) - min(min(obj.cellImg)));
            cellImgNormalized(~obj.cellMask) = NaN;
            if obj.parseInput({varargin,'doPlot',0})
                if obj.parseInput({varargin,'openFig',1});figure;end
                subplot(ceil(sqrt(size(obj.cellEvents,3)+1)),ceil(sqrt(size(obj.cellEvents,3)+1)),1)
                imagesc(cellImgNormalized)
            end
            for k = 1:size(obj.cellEvents,3)
                tmp = obj.cellEvents(:,:,k);
                tmp2 = tmp(:);
                tmp2(~obj.cellMask) = NaN;
                tmp2 = -1 + 2.*(tmp2 - min(tmp2))./(max(tmp2) - min(tmp2)); %Normalize betweeen 0 and 1.
                normalizedImg = reshape(tmp2,size(tmp));
                imMovCorr(k) = obj.myCorr2(cellImgNormalized,normalizedImg);
                
                if obj.parseInput({varargin,'doPlot',0})
                    subplot(ceil(sqrt(size(obj.cellEvents,3)+1)),ceil(sqrt(size(obj.cellEvents,3)+1)),k+1)
                    imagesc(normalizedImg)
                    caxis([0,1])
                    title(['Movie corr: ', num2str(imMovCorr(k))])
                end
            end
        end
        
        function cellStruct = runAllMethods(obj,idx,i,cellStruct,varargin)
            cellStruct.sScore(idx) = obj.getsScore(i,varargin);
            cellStruct.centroidDist(idx) = obj.getcentroidDist(i,varargin);
            cellStruct.globalSNR(idx) = obj.getglobalSNR(i,varargin);
            cellStruct.RMS(idx) = obj.getRMS(i,varargin);
            cellStruct.interCorr(idx) = obj.getinterCorr(i,varargin);
            cellStruct.peakSNR(idx) = obj.getpeakSNR(i,varargin);
            cellStruct.imMovCorr(idx) = {obj.getimMovCorr(i,varargin)};
            CellShapePropieties = obj.getCellShapePropieties(i,varargin);
            cellStruct.CellArea(idx) = CellShapePropieties.Area;
            cellStruct.CellEccentricity(idx) = CellShapePropieties.Eccentricity;
            cellStruct.CellSolidity(idx) = CellShapePropieties.Solidity;
            CellEventPropierties = obj.getCellEventPropierties(i,varargin);
            for propierty = {'Area','Eccentricity','Solidity'}
                valTemp = cellfun(@(x) x.(propierty{1}),CellEventPropierties,'UniformOutput',false);
                emptyVal = cellfun(@isempty ,valTemp);
                val = nan([1,length(valTemp)]);
                val(~emptyVal) = cell2mat(valTemp);
                cellStruct.(['Event',propierty{1}])(idx) = {val};
            end
        end
        
        function [prediction,goodEvents,goodEventsPercent] = predict(obj,i,varargin)
            eventPercent = obj.parseInput({varargin,'eventPercent',0});
            eventMin = obj.parseInput({varargin,'eventMin',0}); 
            sourceMask = obj.parseInput({varargin,'sourceMask',1});
            createCellVariables(obj,i);
            createEventVariables(obj,varargin);
            if eventPercent > 1
                error('neededEvents is a percentage and spans [0,1]')
            end
            goodEvents = size(obj.cellEvents,3);
            goodEventsPercent = goodEvents/length(obj.signalPeakId);
            if goodEventsPercent >= eventPercent && sourceMask && goodEvents >= eventMin
                prediction = 1;
            else
                prediction = 0;
            end
            if obj.parseInput({varargin,'doPlot',0})
                plotEvents(obj,i,[varargin,'prediction',prediction])
            end
        end
        
        function cellStruct = runAllCells(obj,varargin)
            cellList = obj.parseInput({varargin,'cellList',1:obj.cellN});
            cellsWithEvents = find(~cellfun(@isempty,obj.signalPeakIdx));
            cellList = intersect(cellList,cellsWithEvents);
            if strcmp(obj.parseInput({varargin,'methods','all'}),'all')
                cellStruct.sScore = nan([length(cellList),1]);
                cellStruct.centroidDist = nan([length(cellList),1]);
                cellStruct.globalSNR = nan([length(cellList),1]);
                cellStruct.RMS = nan([length(cellList),1]);
                cellStruct.interCorr = nan([length(cellList),1]);
                cellStruct.peakSNR = nan([length(cellList),1]);
                cellStruct.imMovCorr = cell([length(cellList),1]);
                cellStruct.CellArea = nan([length(cellList),1]);
                cellStruct.CellEccentricity = nan([length(cellList),1]);
                cellStruct.CellSolidity = nan([length(cellList),1]);
                cellStruct.EventArea = cell([length(cellList),1]);
                cellStruct.EventEccentricity = cell([length(cellList),1]);
                cellStruct.EventSolidity = cell([length(cellList),1]);
                
                if obj.parseInput({varargin,'showProgress',0});disp('Getting features...');dispTextPrev = '';end
                for i = 1:length(cellList) %#ok<PROPLC>
                    if obj.parseInput({varargin,'showProgress',0})
                          dispText = [num2str(i),' of ',num2str(length(cellList)),' (',num2str(100*i/length(cellList)),'%)'];
                          fprintf(repmat('%c',[1,length(dispTextPrev)]),repmat(8,[1,length(dispTextPrev)]));
                          fprintf(repmat('%c',[1,length(dispText)]),dispText);
                          dispTextPrev = dispText;
                    end
                    if isempty(varargin)
                        varargin = {{'doPlot',0}};
                    end
                    cellStruct = runAllMethods(obj,i,cellList(i),cellStruct,varargin{1}); %#ok<PROPLC>
                end
                if obj.parseInput({varargin,'showProgress',0})
                    fprintf(repmat('%c',[1,length(dispText)]),repmat(8,[1,length(dispText)]));
                end
            else
                error('Method selection for all cell processing not yet implemented')
            end
        end
        
        function predictClassifier(obj,varargin)
            
            obj.setTresholds('tresholds',{'getOverlap','>=0.5','getCellShapePropieties.Eccentricity','<=0.9','getCellShapePropieties.Solidity','>0.8'},'eventPercent',.3);
            valid_strict = obj.validCells;
            obj.setTresholds('tresholds',{'getOverlap','>=0.5','getCellShapePropieties.Eccentricity','<=0.9','getCellShapePropieties.Solidity','>0.8'},'eventPercent',.1);
            valid_gentle = obj.validCells;
            obj.validCells(valid_gentle == 1) = 3;
            obj.validCells(valid_strict == 1) = 1;
            return 
            
            cellList = obj.parseInput({varargin,'cellList',1:length(obj.signalPeakIdx)});
            probTresh = obj.parseInput({varargin,'probTresh',.5});%.9});
            
            cellsWithNoEvents = find(cellfun(@isempty,obj.signalPeakIdx));
            cellList = setdiff(cellList,cellsWithNoEvents);
            obj.validCells(cellsWithNoEvents) = 0;
            
            try load svmModel; catch; error('Trained SVM model not found in path'); end
            
            if obj.parseInput({varargin,'predictNew',0})
                cellList = intersect(find(obj.validCells==3),cellList);    
                if isempty(cellList)
                    disp('None of the specified cells are undecided')
                    return;
                end
            end      
            
            features = obj.runAllCells('cellList',cellList,'showProgress',1);
            
            for field = fields(features)'
                if iscell(features.(field{1}))
                    features.(field{1}) = cellfun(@nanmean,features.(field{1}));
                end
            end
            
            predBasic = rulePrd(features);
            for field = fields(features)'
                features.(field{1}) = features.(field{1})(predBasic);
            end

            featFields = fields(features);
            compactData = zeros([size(features.(featFields{1}),1),length(featFields)]);
            for k = 1:length(featFields)
                compactData(:,k) = features.(featFields{k});
            end
            
            disp('Runing prediction...')
            [~,predProb] = ScoreCVSVMModel.predict(compactData);
            predProb = predProb(:,2);
            
            basicGoodIdx = find(predBasic);            
            goodTreshIdx = basicGoodIdx(predProb >= probTresh)';
            badTreshIdx = basicGoodIdx(predProb < (1-probTresh))';
            ambiguousTreshIdx = basicGoodIdx(predProb < probTresh & predProb >= (1-probTresh))';

            fullDecisions = zeros([1,length(predBasic)]);
            fullDecisions(goodTreshIdx) = 1;
            fullDecisions(ambiguousTreshIdx) = 3;
            obj.validCells(cellList) = fullDecisions;
            disp(['Done: ',num2str(length(goodTreshIdx)),' good cells ('...
                ,num2str(100*length(goodTreshIdx)/length(cellList)),'%), '...
                num2str(length(badTreshIdx) + length(cellList) - length(basicGoodIdx)), ...
                ' bad cells (' num2str(100*(length(badTreshIdx) + length(cellList) - length(basicGoodIdx))/length(cellList))...
                '%) and ', num2str(length(ambiguousTreshIdx)),...
                ' ambiguous (',num2str(100*length(ambiguousTreshIdx)/length(cellList))...
                '%) of a total of ',num2str(length(cellList)),' cells'])
        end
        
        function predProb = predictClassifier_Pasha(obj,varargin)

            cellList = obj.parseInput({varargin,'cellList',1:length(obj.signalPeakIdx)});
            probTresh = obj.parseInput({varargin,'probTresh',.5});%.9});

            cellsWithNoEvents = find(cellfun(@isempty,obj.signalPeakIdx));
            cellList = setdiff(cellList,cellsWithNoEvents);
            obj.validCells(cellsWithNoEvents) = 0;
            
            try load svmModel; catch; error('Trained SVM model not found in path'); end
            
            if obj.parseInput({varargin,'predictNew',0})
                cellList = intersect(find(obj.validCells==3),cellList);    
                if isempty(cellList)
                    disp('None of the specified cells are undecided')
                    return;
                end
            end      
            
            features = obj.runAllCells('cellList',cellList,'showProgress',1);
            
            for field = fields(features)'
                if iscell(features.(field{1}))
                    features.(field{1}) = cellfun(@nanmean,features.(field{1}));
                end
            end
            
%             predBasic = rulePrd(features);
%             predBasic = ones(size(features.sScore));
%             for field = fields(features)'
%                 features.(field{1}) = features.(field{1})(predBasic);
%             end

            featFields = fields(features);
            compactData = zeros([size(features.(featFields{1}),1),length(featFields)]);
            for k = 1:length(featFields)
                compactData(:,k) = features.(featFields{k});
            end
            
            disp('Runing prediction...')

            [~,predProb] = ScoreCVSVMModel.predict(compactData);
            predProb = predProb(:,2);
            
        end
        
        function setTresholds(obj,varargin)
            %Accepts input like obj.predictAll('tresholds',{'getOverlap','>=0.6','getCellShapePropieties.Eccentricity','<=0.9','getCellShapePropieties.Solidity','>0.8'},'eventPercent',.2)
            tresholds = obj.parseInput({varargin,'tresholds',[]});
            if obj.parseInput({varargin,'showProgress',0});dispTextPrev = '';end
            if obj.parseInput({varargin,'predictNew',0})
                c2Pred = find(obj.validCells==3);
            else
                c2Pred = 1:obj.cellN;
            end
            for k = 1:length(c2Pred)
                i = c2Pred(k);
                if ~isempty(obj.signalPeakIdx{i})
                    if obj.parseInput({varargin,'showProgress',0})
                      dispText = [num2str(i),' of ',num2str(length(c2Pred)),' (',num2str(100*i/length(c2Pred)),'%)'];
                      fprintf(repmat('%c',[1,length(dispTextPrev)]),repmat(8,[1,length(dispTextPrev)]));
                      fprintf(repmat('%c',[1,length(dispText)]),dispText);
                      dispTextPrev = dispText;
                    end
                    createCellVariables(obj,i);
                    createEventVariables(obj);
                    if ~isempty(tresholds)
                        eventMask = 1;
                        sourceMask = 1;
                        for k = 1:2:length(tresholds)
                            fun = tresholds{k};
                            field = '';
                            tresh = tresholds{k+1};
                            cutPt = strfind(fun,'.');
                            if cutPt
                                field = fun(cutPt:end);
                                fun = fun(1:cutPt-1);
                            end
                            command = ['obj.',fun,'(',num2str(i),')',field,tresh];
                            mask = eval(command);
                            if length(mask)>1
                                eventMask = eventMask & mask;
                            else
                                sourceMask = sourceMask & mask;
                            end
                        end
                        obj.validCells(i) = obj.predict(i,'eventMask',eventMask,'sourceMask',sourceMask,varargin{:},'doPlot',0);
                    else
                        error('No tresholds specified')
                    end
                else
                    obj.validCells(i) = 0;
                end
            end
            %disp(['Prediction finished: ',num2str(sum(1==obj.validCells)),' good, ',num2str(sum(0==obj.validCells)),' bad, for a good/all ratio of ',num2str(100*sum(1==obj.validCells)/length(obj.validCells)),'%'])
        end
        
        function manualSort(obj,varargin)
            %Accepts input like obj.manualSort('tresholds',{'getOverlap','>=0.5','getCellShapePropieties.Eccentricity','<=0.9','getCellShapePropieties.Solidity','>0.8'},'eventPercent',.2)

            if obj.parseInput({varargin,'openFig',1});fig = figure;else;fig = gcf;end;
            set(fig,'KeyPressFcn',@myCallback);
            i = 1;
            myKey = '';
            stay = 1;
            tresholds = obj.parseInput({varargin,'tresholds',[]});
            cellList = obj.parseInput({varargin,'cellList',1:obj.cellN});
            prediction = [];
            
            videoName = 'filtering';
            frames = 100;
            frameRate = 1;
            v = VideoWriter(videoName);
            v.FrameRate = frameRate;
            v.Quality = 100;
            open(v)
            firstTemp = 1;
            ambiguous = 0;
            good = 0;
            instructions = ({'q: Quit',...
                'g: go to cell',...
                'a: Select ambiguous',...
                't: Select good',...
                'p: Select cell with the cursor',...
                'up: specify cell as good',...
                'down: specify cell as bad',...
                'left: previous cell',...
                'right: next cell'});
            
            
            while stay

                i = mod(i,length(cellList));if i == 0;i=length(cellList);end
                
                if ambiguous
                    cellListAmbiguous = intersect(cellList,find(obj.validCells==3));%1)); %% TO MOVE BETWEEN GOODCELS ONLY
                    if isempty(cellListAmbiguous)
                        warning('No ambiguous cells left')
                        ambiguous = 0;
                    else
                        if ambiguous == 3 || i > length(cellListAmbiguous)
                           ambiguous = 1;
                           minIdx = find(abs(cellListAmbiguous-cellList(i)) == min(abs(cellListAmbiguous-cellList(i))),1,'first');
                           i = minIdx;
                        end
                        j = cellListAmbiguous(i);
                    end
                elseif good
                    cellListGood = intersect(cellList,find(obj.validCells==1));%1)); %% TO MOVE BETWEEN GOODCELS ONLY
                    if isempty(cellListGood)
                        warning('No good cells left')
                        good = 0;
                    else
                        if good == 3 || i > length(cellListGood)
                           good = 1;
                           minIdx = find(abs(cellListGood-cellList(i)) == min(abs(cellListGood-cellList(i))),1,'first');
                           i = minIdx;
                        end
                        j = cellListGood(i);
                    end
                else
                    j = cellList(i);
                end
                createCellVariables(obj,j);
                createEventVariables(obj);
                if ~isempty(tresholds)
                    eventMask = 1;
                    sourceMask = 1;
                    for k = 1:2:length(tresholds)
                        fun = tresholds{k};
                        field = '';
                        tresh = tresholds{k+1};
                        cutPt = strfind(fun,'.');
                        if cutPt
                            field = fun(cutPt:end);
                            fun = fun(1:cutPt-1);
                        end
                        command = ['obj.',fun,'(',num2str(j),')',field,tresh];
                        try
                            mask = eval(command);
                        catch
                            error([fun,' ',tresh,' is not a valid treshold'])
                        end
                        if length(mask)>1
                            eventMask = eventMask & mask;
                        else
                            sourceMask = sourceMask & mask;
                        end
                    end
                    if obj.parseInput({varargin,'predict',0})
                        prediction = obj.predict(j,'eventMask',eventMask,'sourceMask',sourceMask,varargin{:},'doPlot',0);
                    end
                else
                   eventMask = [];
                   sourceMask = [];
                end
            
%                 set(0,'CurrentFigure',fig)
                obj.plotEvents(j,'openFig',0,'eventMask',eventMask,'sourceMask',sourceMask,'prediction',prediction);drawnow
                if ambiguous
                    annotation('textbox',[0,0,1,1],'String','Only ambiguous','FitBoxToText','on','Color','red');
                end
                if good
                    annotation('textbox',[0,0,1,1],'String','Only good','FitBoxToText','on','Color','red');
                end
                annotation('textbox',[0,0,.6,.6],'String',instructions,'FitBoxToText','on');
                if ~firstTemp
                    try
                        set(gcf,'units','normalized','outerposition',[0 0 1 1])                                      
                        v.writeVideo(getframe(gcf));
                    catch
                        warning('Resizing window...')
                    end
                end
                firstTemp = 0;

                while true
                    switch myKey
                        case ''
                            pause(0.001)
                        case 'rightarrow'
                            i = i + 1;
                            break
                        case 'uparrow'
                            obj.validCells(j) = 1;
                            if ~ambiguous && ~good
                                i = i + 1;
                            end
                            break
                        case 'downarrow'
                            obj.validCells(j) = 0;
                            if ~ambiguous && ~good
                                i = i + 1;
                            end
                            break
                        case 'leftarrow'
                            i = i - 1;
                            break
                        case 'p'
                            %figure;imagesc(squeeze(max(obj.inputImages)));
                            [x,y] = ginput(1);
                            window_size = 10;
                            window_x = (floor(x) - window_size):(floor(x) + window_size);
                            window_y = (floor(y) - window_size):(floor(y) + window_size);
                            bad_idx = window_x<1 | window_y<1 | window_x>size(obj.inputImages,3) | window_y>size(obj.inputImages,2);
                            window_x(bad_idx) = [];
                            window_y(bad_idx) = [];
                            [X,Y] = meshgrid(window_x,window_y);
                            try
                                [~,i] = max(squeeze(mean(mean(obj.inputImagesT(:,Y(:),X(:)),2),3)));
                            catch
                                warning('Issue setting the cell')
                            end
                            %imagesc(squeeze(obj.inputImages(i,:,:)));hold on;plot(x,y,'r.','markersize',20)
                            %title([num2str(val),' ',num2str(i),' ',num2str(round(x)),' ' ,num2str(round(y))])
                            break
                        case 'q'
                            stay = 0;
                            break
                        case 't'
                            if ambiguous
                                ambiguous = 0;
                            end
                            if good
                                good = 0;
                                i = find(cellList==cellListGood(i));
                            else
                                good = 3;
                            end
                            break
                        case 'a'
                            if good
                                good = 0;
                            end
                            if ambiguous
                                ambiguous = 0;
                                i = find(cellList==cellListAmbiguous(i));
                            else
                                ambiguous = 3;
                            end
                            break
                        case 'g'
                            options.Interpreter = 'tex';
                            tempI = str2double(inputdlg('Enter cell number','Celln',[1 40],{''},options));
                            if ~isempty(tempI)
                                i = tempI;
                            end
                            if ambiguous
                                if ~sum(cellListAmbiguous == i)
                                    ambiguous = 0;
                                else
                                    ambiguous = 3;
                                end
                            end
                            if good
                                if ~sum(cellListGood == i)
                                    good = 0;
                                else
                                    good = 3;
                                end
                            end
                            break
                        otherwise
                            pause(0.001)
                    end
                end
                myKey = '';
                cla('reset')
            end
            close(v)
            close(gcf)
            function myCallback(~,chinf)
                myKey = chinf.Key;
            end
        end
        
        function plotEvents(obj,i,varargin)
            createCellVariables(obj,i);
            createEventVariables(obj);
            if ~isempty(varargin) && iscell(varargin{1})
                varargin = varargin{1};
            end
            eventMask = obj.parseInput({varargin,'eventMask',[]});
            sourceMask = obj.parseInput({varargin,'sourceMask',[]});
            prediction = obj.parseInput({varargin,'prediction',obj.validCells(i)});
            
            %Create the image montage
            [croppedPeakImages2,infoGrid2] = obj.getMontage(cat(3,obj.cellImg,obj.cellEvents),obj.maxValMovie,obj.maxShowEvents,obj.cellLine,eventMask,sourceMask);
            infoGridR = infoGrid2;
            infoGridR(infoGrid2==1) = 173/255;
            infoGridR(infoGrid2==2) = 174/255;
            infoGridR(infoGrid2==3) = 0;
            infoGridG = infoGrid2;
            infoGridG(infoGrid2==1) = 1;
            infoGridG(infoGrid2==2) = 0;
            infoGridG(infoGrid2==3) = 191/255;
            infoGridB = infoGrid2;
            infoGridB(infoGrid2==1) = 47/255;
            infoGridB(infoGrid2==2) = 45/255;
            infoGridB(infoGrid2==3) = 1;
            infoGrid = cat(3,infoGridR,infoGridG,infoGridB);
            infoGrid(isnan(infoGrid(:))) = 0;
            infoGridMask = ones(size(infoGrid2));
            infoGridMask(isnan(infoGrid2)) = 0;
            
            if obj.parseInput({varargin,'openFig',1});h = figure('units','normalized','outerposition',[0 0 1 1]);
                ;end
            clf
            
            s1 = subplot(3,4,[1,2,5,6]);
            imagesc(croppedPeakImages2);colormap(s1,'default');hold on;
            h = imshow(infoGrid);
            set(h, 'AlphaData', infoGridMask*.7);
            hold off;
            if obj.parseInput({varargin,'showInfo',1})
                cellProp = table2cell(obj.getCellShapePropieties(i));
                title(['Area: ',num2str(cellProp{1}),'    Eccentricity: ',num2str(cellProp{2}),'    Solidity: ',num2str(cellProp{3}),'    Event correlation: ',num2str(obj.getinterCorr(i))])
            end
            
            s2 = subplot(3,4,[3,4,7,8]);

            redC = obj.validCells==0;redC(i) = 0;
            greenC = obj.validCells==1;greenC(i) = 0;
            grayC = obj.validCells==3;grayC(i) = 0;
            blueC = i;

            redImg =(squeeze(max(obj.inputImagesT(redC,:,:),[],1)));
            greenImg = (squeeze(max(obj.inputImagesT(greenC,:,:),[],1)));
            grayImg = (squeeze(max(obj.inputImagesT(grayC,:,:),[],1)));
            blueImg = (squeeze(obj.inputImagesT(blueC,:,:)));
            if isempty(redImg);redImg = zeros([size(obj.inputImagesT,2),size(obj.inputImagesT,3)]);end
            if isempty(greenImg);greenImg = zeros([size(obj.inputImagesT,2),size(obj.inputImagesT,3)]);end
            if isempty(grayImg);grayImg = zeros([size(obj.inputImagesT,2),size(obj.inputImagesT,3)]);end
            
            redCellMap = 173/255*greenImg + 174/255*redImg + 0.7*grayImg;
            greenCellMap = greenImg + 0.7*grayImg + 191/255*blueImg;
            blueCellMap = 47/255*greenImg + 45/255*redImg + blueImg + 0.7*grayImg;
            
            colorCellMap = cat(3,redCellMap,greenCellMap,blueCellMap);
            colorCellMap = colorCellMap/max(colorCellMap(:));
%             whiteFrame = cat(3,ones(size(redCellMap)),ones(size(redCellMap)),ones(size(redCellMap)));
%             whiteMask = sum(colorCellMap,3);
%             whiteMask(whiteMask~=0) = 1;
            imshow(colorCellMap);hold on;
            halfSide = obj.cropSizeLength/2;
%             [y,x] = find(squeeze(obj.inputImage)==max(max(squeeze(obj.inputImage))));
            [x y] = findCentroid(squeeze(obj.inputImage),'waitbarOn',0);
            rectangle('Position',[x-halfSide,y-halfSide,obj.cropSizeLength,obj.cropSizeLength],'linewidth',1,'EdgeColor',[0,191/255,1]);       

            
%             h = imshow(whiteFrame);
%             set(h,'AlphaData',1-whiteMask)
            
            if obj.parseInput({varargin,'showInfo',1})
                title(['Cell ',num2str(i),' of ',num2str(obj.cellN),' (',num2str(sum(obj.validCells==1)),' good, ',num2str(sum(obj.validCells==0)),' bad, ',num2str(sum(obj.validCells==3)),' unknown)'])
%                 title(['Center-dist: ', num2str(obj.getcentroidDist(i))])
            end
            
            s3 = subplot(3,4,[9,10]);
            plot(obj.inputSignal);
            yMax = max(obj.inputSignal(:));
%             yMax = obj.calHeight;
%             if yMax < max(obj.inputSignal(:)); yMax = max(obj.inputSignal(:));end
            ylim([0,yMax])
            hold on
            if ~isempty(eventMask)
                plot(obj.eventPeaks(eventMask),obj.inputSignal(obj.eventPeaks(eventMask)),'.','color',[173/255 1 47/255],'markersize',15)
                plot(obj.eventPeaks(~eventMask),obj.inputSignal(obj.eventPeaks(~eventMask)),'.','color',[174/255,0,45/255],'markersize',15)
            else
                plot(obj.eventPeaks,obj.inputSignal(obj.eventPeaks),'c.','markersize',15)
            end
            if obj.parseInput({varargin,'showInfo',1})
            title(['Power of the events: ', num2str(obj.getRMS(i)),'    Event SNR: ',num2str(obj.getglobalSNR(i))])
            end
            hold off
            
            s4 = subplot(3,4,[11,12]);
            spikeCenterTrace = obj.computeSpikeCenterTrace;
            midleIdx = size(spikeCenterTrace,2)/2;
            lLim = floor(midleIdx-obj.slopeFrameWindow);
            rLim = ceil(midleIdx+obj.slopeFrameWindow);
            if lLim < 1; lLim = 1;end
            if rLim > size(spikeCenterTrace,2);rLim = size(spikeCenterTrace,2);end
            spikeWindow = spikeCenterTrace(:,lLim:rLim);
            hold on
            
            for k = size(spikeWindow,1):-1:1
                try
                    plot(1:size(spikeWindow,2),spikeWindow(k,:),'color',(log(k)./log(size(spikeWindow,1))).*[1,1,1]);
                catch
                    warning('Too few events detected')
                end
            end
            plot(1:size(spikeWindow,2),mean(spikeWindow,1)+std(spikeWindow),'r--','linewidth',1);
            plot(1:size(spikeWindow,2),mean(spikeWindow,1)-std(spikeWindow),'r--','linewidth',1);
            plot(1:size(spikeWindow,2),mean(spikeWindow,1),'r','linewidth',2);
            yMax = max(obj.inputSignal(:));
%             yMax = obj.calHeight;
%             if yMax < max(spikeWindow(:)); yMax = max(spikeWindow(:));end
            line([1+size(spikeWindow,2)/2 1+size(spikeWindow,2)/2],[0,yMax],'linewidth',2)
            xlim([1,size(spikeWindow,2)])
            ylim([0,yMax])
            hold off
            if obj.parseInput({varargin,'showInfo',1})
%             title(['Up/downrising ratio: ', num2str(obj.getexpRatio(i)),'    Waveform SNR: ',num2str(obj.getpeakSNR(i))])
            title(['Waveform SNR: ',num2str(obj.getpeakSNR(i))])
            end         
%             warning ('off')
%             if obj.parseInput({varargin,'openFig',1});figure;end
%             ax1 = subplot(2,2,1);
%             imagesc(croppedPeakImages2);colormap(ax1,[colormap;[0,1,0];[1,0,0]]);
%             subplot(ceil(sqrt(size(obj.cellEvents,3)+1)),ceil(sqrt(size(obj.cellEvents,3)+1)),1)
%             imagesc(obj.cellImg);hold on
%             if ~isempty(sourceMask)
%                 if sourceMask
%                     rectangle('Position',[1,1,size(obj.cellImg,1)-1,size(obj.cellImg,2)-1],'linewidth',5,'EdgeColor','g');       
%                 else
%                     rectangle('Position',[1,1,size(obj.cellImg,1)-1,size(obj.cellImg,2)-1],'linewidth',5,'EdgeColor','r'); 
%                 end
%             end
%             for k = 1:size(obj.cellEvents,3)
%                 subplot(ceil(sqrt(size(obj.cellEvents,3)+1)),ceil(sqrt(size(obj.cellEvents,3)+1)),k+1)
%                 imagesc(obj.cellEvents(:,:,k)')
%                 caxis([0,max(obj.cellEvents(:))])
%                 hold on;
%                 if ~isempty(eventMask)
%                     if eventMask(k)
%                         rectangle('Position',[1,1,size(obj.cellEvents,1)-1,size(obj.cellEvents,2)-1],'linewidth',5,'EdgeColor','g');
%                     else
%                         rectangle('Position',[1,1,size(obj.cellEvents,1)-1,size(obj.cellEvents,2)-1],'linewidth',5,'EdgeColor','r'); 
%                     end   
%                 end
%             end
            fig = gcf;
            if obj.validCells(i)~=3
                prediction = obj.validCells(i);
            else
                fig.Color = [1,1,1];
            end
            
            if ~isempty(prediction)
                if prediction
                    fig.Color = [154/255,205/255,50/255];
                else
                    fig.Color = [240/255,128/255,128/255];
                end
            end
            
            warning ('on')
        end
                              
    end  
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%PRIVATE STUFF%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%    
        
    methods(Access = private)
        
        function obj = createCellVariables(obj,i)
            if isempty(obj.i) || i ~= obj.i  %Dont redo this if the same (one) cell is asked
                obj.i = i;
                obj.resetEvents;
                obj.createOptimizationObj;
                %Set selected input/objects
                obj.inputImage = obj.inputImages(i,:,:);
                obj.inputSignal = obj.inputSignals(i,:);
                obj.signalPeakId = obj.signalPeakIdx{i};

                %Get cropped images at peak value
                obj.croppedPeakImages = obj.viewMontageCustom(obj.inputMovie,obj.inputImage,obj.inputSignal,obj.signalPeakId,obj.minValMovie,obj.maxValMovie,obj.cropSizeLength,true);

                %Get smoothing filter
                obj.PSF = fspecial('gaussian',ceil(size(obj.croppedPeakImages,1)*obj.filterSpaceRatio),obj.filterHeight);

                %Get filter and treshold cell
                obj.cellImg = obj.croppedPeakImages(:,:,1);
                if obj.smoothImages
                    obj.cellImgF = filter2(obj.PSF,obj.cellImg);
                else
                    obj.cellImgF = obj.cellImg;
                end
                obj.cellMask = obj.cellImgF>nanmedian(obj.cellImgF(:))+0*nanstd(obj.cellImgF(:));%obj.cellImgF~=0;
                obj.cellLine = bwperim(obj.cellMask);
                obj.cellLine(obj.cellImgF==max(obj.cellImgF(:))) = 1;
            end
        end  
        
        function obj = resetEvents(obj)
            [obj.cellEvents,obj.eventPeaks,obj.cellEventsF,obj.cellEventsTimeAvg,obj.meanValVect,obj.tresholdedImg,obj.tresholdedOutside] = deal([]);
        end
        
        function obj = createEventVariables(obj,varargin)
            
            if ~isempty(varargin)
                newMask = obj.parseInput({varargin{1},'eventMask',[]});
            else
                newMask = [];
            end
            
            if isempty(obj.cellEvents) || ~isequal(obj.eventMask,newMask) %Events are always cleared when cell cahnges
                obj.eventMask = newMask;
                allEvents = obj.croppedPeakImages(:,:,2:end);
                obj.createOptimizationObj;
                if obj.timeAverageEvents
                    %weightVectR = exp((1:obj.windowSize)-1/(obj.windowSize)-1)/(1.5*max(exp((1:obj.windowSize)-1/(obj.windowSize)-1)));
                    weightVectR = ((1:obj.timeWindowSize).^obj.timeWeightDecay)/max(1.2*((1:obj.timeWindowSize).^obj.timeWeightDecay));
                    weightVect = [weightVectR,1,flip(weightVectR)];

                    windowFrames = bsxfun(@plus,obj.signalPeakId',-obj.timeWindowSize:obj.timeWindowSize);
                    windowFrames(windowFrames<=0)=1;
                    windowFrames(windowFrames>=size(obj.inputMovie,3))=size(obj.inputMovie,3);

                    allFramesLine = obj.viewMontageCustom(obj.inputMovie,obj.inputImage,obj.inputSignal,reshape(windowFrames,[1,numel(windowFrames)]),obj.minValMovie,obj.maxValMovie,obj.cropSizeLength,false);
                    allFrames = reshape(allFramesLine(:,:,2:end),[size(allFramesLine,1),size(allFramesLine,2),length(obj.signalPeakId),length(weightVect)]);

                    eventsTimeAvgMethod = ones(size(obj.croppedPeakImages(:,:,2:end)));
                    %figure;
                     for k = 1:size(allFrames,3)
                        %subplot(ceil(sqrt(size(allFrames,3))),ceil(sqrt(size(allFrames,3))),k)
                        averageMatTemp = repmat(reshape(weightVect,[1,1,numel(weightVect)]),[size(allFrames,1),size(allFrames,2),1]).*squeeze(allFrames(:,:,k,:));
                        averageMat = max(averageMatTemp(:))*(sum(averageMatTemp,3)/max(max(sum(averageMatTemp,3))));
                        eventsTimeAvgMethod(:,:,k) = averageMat;
                        %imagesc(averageMat)
                        %caxis([0,max(croppedPeakImages(:))])
                     end

                     allEvents = eventsTimeAvgMethod;
                end

                if obj.smoothImages
                    allEventsFMethod = zeros(size(allEvents));
                    for j = 1:size(allEvents,3)
                        allEventsFMethod(:,:,j) = filter2(obj.PSF,allEvents(:,:,j));
                    end
                    allEvents = allEventsFMethod;
                end
                
                %Sort and filter events
                if isempty(obj.eventMask)
                    mask = 1:size(allEvents,3);
                else
                    mask = obj.eventMask;
                end
                [allEvents,meanValVectMethod] = obj.sortEventsByIntensity(allEvents,obj.cellMask);
                [tresholdedImgMethod,tresholdedOutsideMethod] = obj.getEventTreshold(obj.cellImgF,allEvents,meanValVectMethod);
                obj.meanValVect = meanValVectMethod(mask);
                obj.cellEvents = allEvents(:,:,mask);
                obj.eventPeaks = obj.signalPeakId(mask);
                if exist('allEventsFMethod','var')
                    obj.cellEventsF = allEventsFMethod(:,:,mask);
                end
                if exist('eventsTimeAvgMethod','var')
                    obj.cellEventsTimeAvg = eventsTimeAvgMethod(:,:,mask);
                end
                obj.tresholdedImg = tresholdedImgMethod(:,:,mask);
                obj.tresholdedOutside = tresholdedOutsideMethod(:,:,mask);
            end
        end
        
        function obj = createOptimizationObj(obj)
            obj.optimization = struct('peakIdx', [],'spikeCenterTrace',[]);
        end
        
        function peakIdx = computePeakIdx(obj)
            if ~isempty(obj.optimization.peakIdx)
                peakIdx = obj.optimization.peakIdx;
            else
                peakIdx = bsxfun(@plus,obj.timeSeqGlobalSNR',obj.eventPeaks);
                peakIdx = unique(peakIdx(:));
                peakIdx = peakIdx(peakIdx<=length(obj.inputSignal) & peakIdx>0);
                obj.optimization.peakIdx = peakIdx;
            end
        end  
        
        function spikeCenterTrace = computeSpikeCenterTrace(obj)
            if ~isempty(obj.optimization.spikeCenterTrace)
                spikeCenterTrace = obj.optimization.spikeCenterTrace;
            else
                extractMatrix = bsxfun(@plus,obj.eventPeaks',obj.spikeROI);
                extractMatrix(extractMatrix<=0)=1;
                extractMatrix(extractMatrix>=size(obj.inputSignal,2))=size(obj.inputSignal,2);
                spikeCenterTrace = reshape(obj.inputSignal(extractMatrix),size(extractMatrix)); 
                obj.optimization.spikeCenterTrace = spikeCenterTrace;
            end
        end
        
    end
    
    methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%SET METHODS %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
        function obj = set.timeAverageEvents(obj,value)
            obj.timeAverageEvents = value;
            obj.variableChanged;
        end
        function obj = set.timeWindowSize(obj,value)
            obj.timeWindowSize = value;
            obj.variableChanged;
        end
        function obj = set.timeWeightDecay(obj,value)
            obj.timeWeightDecay = value;
            obj.variableChanged;
        end       
        function obj = set.smoothImages(obj,value)
            obj.smoothImages = value;
            obj.variableChanged;
        end
        function obj = set.cropSizeLength(obj,value)
            obj.cropSizeLength = value;
            obj.variableChanged;
        end
        function obj = set.filterSpaceRatio(obj,value)
            obj.filterSpaceRatio = value;
            obj.variableChanged;
        end        
        function obj = set.filterHeight(obj,value)
            obj.filterHeight = value;
            obj.variableChanged;
        end           
        function obj = set.timeSeqGlobalSNR(obj,value)
            obj.timeSeqGlobalSNR = value;
            obj.variableChanged;
        end          
        function obj = set.spikeROI(obj,value)
            obj.spikeROI = value;
            obj.variableChanged;
        end   
        function obj = set.slopeFrameWindow(obj,value)
            obj.slopeFrameWindow = value;
            obj.variableChanged;
        end                    
        function obj = variableChanged(obj)
            if ~isempty(obj.i)
                obj.resetEvents;
                createCellVariables(obj,obj.i);
                createEventVariables(obj);
            end
        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%STATIC METHODS %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods(Static, Access = private)
        
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
                        varargout{k} = defaultVal{k};
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
        
        function [allPeakImages,meanValVect] = sortEventsByIntensity(allPeakImages,cellMask)
            meanValVect = zeros([1,size(allPeakImages,3)]);
            meanOutVect = zeros([1,size(allPeakImages,3)]); 
            for j = 1:size(allPeakImages,3)
                tempEvent = allPeakImages(:,:,j);
                meanValVect(j) = nanmean(tempEvent(cellMask));
                meanOutVect(j) = nanmean(tempEvent(~cellMask));
            end
            [meanValVect,sortingIdx] = sort(meanValVect,'descend');
            allPeakImages = allPeakImages(:,:,sortingIdx);           
        end
        
        function [tresholdedImg,tresholdedOutside] = getEventTreshold(cellImgF,allPeakImages,meanValVect)
            tresholdedImg = zeros(size(allPeakImages));
            tresholdedOutside = zeros(size(allPeakImages));
            for k = 1:size(allPeakImages,3)
                tempImg = allPeakImages(:,:,k);
                tresholdedImg(:,:,k) = meanValVect(k)*((tempImg>nanmean(tempImg(:))+1*nanstd(tempImg(:)))).*(cellImgF>nanmean(cellImgF(:))+1*nanstd(cellImgF(:)));
                tresholdedOutside(:,:,k) = (tempImg>nanmean(tempImg(:))+1*nanstd(tempImg(:))) .* (1-tresholdedImg(:,:,k));
            end
        end
            
        function c2 = myCorr2(A,B)
            AMean = nanmean(A(:));
            BMean = nanmean(B(:));
            c2 = nansum((A(:)-AMean).*(B(:)-BMean))./sqrt(nansum((A(:)-AMean).^2).*nansum((B(:)-BMean).^2)); 
        end

        function [croppedPeakImages] = viewMontageCustom(inputMovie,inputImage,thisTrace,signalPeakArray,minValMovie,maxValMovie,cropSizeLength,doSort)
            if isempty(signalPeakArray)
                warning('No events detected')
%                 inputImageT = inputImage(inputImage>mean(inputImage)+std(inputImage));
                croppedPeakImages = inputImage;
                return
            end
            % signalPeakArray
            maxSignalsToShow = 1e16;
            if doSort
                peakSignalAmplitude = thisTrace(signalPeakArray(:));
            % peakSignalAmplitude
                [peakSignalAmplitude, peakIdx] = sort(peakSignalAmplitude,'descend');
                % peakSignalAmplitude
                signalPeakArray = signalPeakArray(peakIdx);
            end
            if length(signalPeakArray)>maxSignalsToShow
                % choose a random subset
                signalPeakArray = signalPeakArray(1:maxSignalsToShow);
            end

        %     signalPeakArray((signalPeakArray-31)<1) = [];
        %     signalPeakArray((signalPeakArray+31)>length(thisTrace)) = [];

            signalPeakArray = {signalPeakArray};
            % signalPeakArray
            croppedPeakImages = compareSignalToMovieCustom(inputMovie, inputImage, thisTrace,'getOnlyPeakImages',1,'waitbarOn',0,'extendedCrosshairs',0,'signalPeakArray',signalPeakArray,'cropSize',cropSizeLength,'crosshairs',0);
            
            function [croppedPeakImages] = compareSignalToMovieCustom(inputMovie, inputFilters, inputSignal, varargin)
                % shows a cropped version of inputMovie for each inputFilters and aligns it to inputSignal peaks to make sure detection is working
                % biafra ahanonu
                % started: 2013.11.04 [18:40:45]
                % inputs
                    % inputMovie - matrix dims are [X Y t] - where t = number of time points
                    % inputFilters - matrix dims are [n X Y] - where n = number of filters, NOTE THE DIFFERENCE
                    % inputSignal - matrix dims are [n t] - where n = number of signals, t = number of time points
                % outputs
                    % none, this is a display function
                % changelog
                    % 2014.01.18 [12:24:29] fully implemented, cut out from controllerAnalysis, need to improve handling at beginning of movie, but that's a playMovie function issue
                % TODO
                    %

                %========================
                % old way of saving, only temporary until full switch
                options.oldSave = 0;
                % size in pixels to show signal image
                options.cropSize = 20;
                % frames before/after to show
                options.timeSeq = -10:10;
                % waitbar
                options.waitbarOn = 1;
                % whether to just get the peak images and ignore showing the movie
                options.getOnlyPeakImages = 0;
                % 1 = plus shaped crosshairs, 0 = dot
                options.extendedCrosshairs = 1;
                %
                options.crosshairs = 1;
                %
                options.signalPeakArray = [];
                % set to 1 if input images should be normalized
                options.normalizeMovieImages = 0;
                % get options
                options = getOptions(options,varargin);
                % unpack options into current workspace
                % fn=fieldnames(options);
                % for i=1:length(fn)
                %     eval([fn{i} '=options.' fn{i} ';']);
                % end
                %========================

                if isempty(options.signalPeakArray)
                    [signalPeaks, signalPeakArray] = computeSignalPeaks(inputSignal, 'makePlots', 0,'makeSummaryPlots',0,'waitbarOn',options.waitbarOn);
                else
                    signalPeakArray = options.signalPeakArray;
                end

                % get the centroids and other info for movie
                [xCoords yCoords] = findCentroid(squeeze(inputFilters),'waitbarOn',options.waitbarOn);
                cropSize = options.cropSize;
                nSignals = size(inputFilters,1);
                nPoints = size(inputMovie,3);
                movieDims = size(inputMovie);
                timeSeq = options.timeSeq;

                % inputMovie(inputMovie>1.3) = NaN;
                % inputMovie(inputMovie<0.8) = NaN;

                % loop over all signals and visualize their peaks side-by-side with movie
                exitSignal = 0;
                for signalNo=1:nSignals
                    peakLocations = signalPeakArray{signalNo};

                    peakIdxs = bsxfun(@plus,timeSeq',peakLocations);
                    peakIdxs(find(peakIdxs<1)) = 1;
                    peakIdxs(find(peakIdxs>nPoints)) = 1;
                    % get region to crop
                    warning off;
                    xLow = xCoords(signalNo) - cropSize;
                    xHigh = xCoords(signalNo) + cropSize;
                    yLow = yCoords(signalNo) - cropSize;
                    yHigh = yCoords(signalNo) + cropSize;
                    % check that not outside movie dimensions
                    xMin = 0;
                    xMax = movieDims(2);
                    yMin = 0;
                    yMax = movieDims(1);

                    % adjust for the difference in centroid location if movie is cropped
                    xDiff = 0;
                    yDiff = 0;
                    if xLow<=xMin; xDiff = xLow-xMin; xLow = xMin+1; end
                    if xHigh>=xMax; xDiff = xHigh-xMax; xHigh = xMax-1; end
                    if yLow<=yMin; yDiff = yLow-yMin; yLow = yMin+1; end
                    if yHigh>=yMax; yDiff = yHigh-yMax; yHigh = yMax-1; end

                    % need to add a way to adjust the cropped movie target point if near the boundary

                    % get the cropped movie at peaks
                    % yLow
                    % yHigh
                    % xLow
                    % xHigh
                    % peakLocations
                    croppedPeakImages = inputMovie(yLow:yHigh,xLow:xHigh,peakLocations);

                    globalCoord = [xCoords(signalNo),yCoords(signalNo)];
                    %croppedPeakImages = shiftEvents(croppedPeakImages,inputMovie,peakLocations,cropSize,movieDims,globalCoord);

                    firstImg = squeeze(inputFilters(signalNo,yLow:yHigh,xLow:xHigh));
                    firstImg = normalizeVector(firstImg,'normRange','zeroToOne');
                    maxVec = nanmax(croppedPeakImages(:));
                    minVec = nanmin(croppedPeakImages(:));
                    firstImg = firstImg*maxVec;
                    % firstImg = (firstImg-minVec)./(maxVec-minVec);
                    % if options.normalizeMovieImages==1
                    % end
                    croppedPeakImages(:,:,end+1) = firstImg;
                    % move inputImage to the front
                    croppedPeakImages = circshift(croppedPeakImages,[0 0 1]);
                    % croppedPeakImagesTmp = croppedPeakImages(:,:,end);
                    % croppedPeakImagesTmp(:,:,end+1:end+(length(croppedPeakImages)-1)) = croppedPeakImages(:,:,1:(end-1));
                    % croppedPeakImages = croppedPeakImagesTmp;
                    for frameNo=1:size(croppedPeakImages,3)
                        cropImg = squeeze(croppedPeakImages(:,:,frameNo));
                        if options.normalizeMovieImages==1
                            cropImg = normalizeVector(cropImg,'normRange','zeroToOne');
                        end
                        croppedPeakImages(:,:,frameNo) = cropImg;
                    end
                    % croppedPeakImages = normalizeMovie(croppedPeakImages,'normalizationType','meanDivision');
                    cDims = size(croppedPeakImages);
                    if options.getOnlyPeakImages==0
                        % get cropped version of the movie
                        croppedMovie = inputMovie(yLow:yHigh,xLow:xHigh,peakIdxs);
                        cDims = size(croppedMovie);
                        exitSignal = playMovie(inputMovie(:,:,peakIdxs),'extraMovie',croppedMovie,...
                            'extraLinePlot',inputSignal(signalNo,peakIdxs),...
                            'windowLength',30,...
                            'colormapColor','jet',...
                            'extraTitleText',['signal #' num2str(signalNo) '/' num2str(nSignals) '    peaks: ' num2str(length(peakLocations))],...
                            'primaryPoint',[xCoords(signalNo) yCoords(signalNo)],...
                            'secondaryPoint',crossHairLocation);
                            % 'recordMovie','test.avi',...
                    end
                    warning on;
                    if exitSignal==1
                        break;
                    end
                end
            end  
        end
        
        function [croppedPeakImagesM,infoGridM] = getMontage(croppedPeakImages,maxValMovie,maxSignalsToShow,shapeLine,eventMask,sourceMask)
            if size(croppedPeakImages,3)-1>maxSignalsToShow
                warning(['More than ',num2str(maxSignalsToShow),' events, ignoring weaker'])
                croppedPeakImages = croppedPeakImages(:,:,1:maxSignalsToShow+1);
            end
            meanTransientImage = nanmean(croppedPeakImages,3);
            croppedPeakImages = cat(3,croppedPeakImages(:,:,1),meanTransientImage,croppedPeakImages(:,:,2:end));
            infoGrid = croppedPeakImages;
            if ~isempty(sourceMask)
                if sourceMask
                    infoGrid(:,:,1) = padarray(croppedPeakImages(2:end-1,2:end-1,1),[1 1],maxValMovie);
                    infoGrid(:,:,1) = infoGrid(:,:,1).*(1-shapeLine) + shapeLine*maxValMovie;
                else
                    infoGrid(:,:,1) = padarray(croppedPeakImages(2:end-1,2:end-1,1),[1 1],maxValMovie*2);
                    infoGrid(:,:,1) = infoGrid(:,:,1).*(1-shapeLine) + shapeLine*maxValMovie*2;
                end
            else
                infoGrid(:,:,1) = padarray(croppedPeakImages(2:end-1,2:end-1,1),[1 1],maxValMovie*3);
                infoGrid(:,:,1) = infoGrid(:,:,1).*(1-shapeLine) + shapeLine*maxValMovie*3;                
            end
            
            infoGrid(:,:,2) = padarray(croppedPeakImages(2:end-1,2:end-1,2),[1 1],maxValMovie*3);
            infoGrid(:,:,2) = infoGrid(:,:,2).*(1-shapeLine) + shapeLine*maxValMovie*3;
            if ~isempty(eventMask)
                for k = 3:size(infoGrid,3)
                    if eventMask(k-2)
                        infoGrid(:,:,k) = padarray(croppedPeakImages(2:end-1,2:end-1,k),[1 1],maxValMovie);
                        infoGrid(:,:,k) = infoGrid(:,:,k).*(1-shapeLine) + shapeLine*maxValMovie;
                    else
                        infoGrid(:,:,k) = padarray(croppedPeakImages(2:end-1,2:end-1,k),[1 1],maxValMovie*2); 
                        infoGrid(:,:,k) = infoGrid(:,:,k).*(1-shapeLine) + shapeLine*maxValMovie*2;
                    end
                end
            else
                for k = 3:size(infoGrid,3)
                    infoGrid(:,:,k) = padarray(croppedPeakImages(2:end-1,2:end-1,k),[1 1],maxValMovie*3);
                    infoGrid(:,:,k) = infoGrid(:,:,k).*(1-shapeLine) + shapeLine*maxValMovie*3;  
                end
            end 
            infoGrid(mod(infoGrid(:),maxValMovie)~=0) = nan;
            infoGrid = infoGrid./maxValMovie;
            
%             if size(croppedPeakImages,3)>maxSignalsToShow
%                 dimDiff = maxSignalsToShow-size(croppedPeakImages,3);
%                 croppedPeakImagesTmp = NaN([size(croppedPeakImages,1) size(croppedPeakImages,2) dimDiff]);
%                 croppedPeakImages = cat(3,croppedPeakImages,croppedPeakImagesTmp);
%             end
            %Separating the mean and source form events
            blackFrames = ceil(sqrt(size(croppedPeakImages,3)))-2;
            croppedPeakImages = cat(3,croppedPeakImages(:,:,1:2),zeros([size(croppedPeakImages,1),size(croppedPeakImages,2),blackFrames]),croppedPeakImages(:,:,3:end));
            infoGrid = cat(3,infoGrid(:,:,1:2),zeros([size(infoGrid,1),size(infoGrid,2),blackFrames]),infoGrid(:,:,3:end));
            warning off
            montage(permute(croppedPeakImages(:,:,:,1),[1 2 4 3]))
            croppedPeakImagesM = getimage;
            montage(permute(infoGrid(:,:,:,1),[1 2 4 3]))
            infoGridM = getimage;
            % change zeros to ones, fixes range of image display
%             croppedPeakImages2(croppedPeakImages2==0)=NaN;
            % croppedPeakImages2(1,1) = minValMovie;
%             croppedPeakImages2(1,1) = 0;
%             croppedPeakImages2(1,2) = maxValMovie*0.4;
%             imagesc(croppedPeakImages2);
%             customColors = customColormap([]);
%             colormap(customColors);
%             axis off;
            % title('frames at signal peaks, press any key to exit');
            % ginput(1);
            % close(2);figure(mainFig);
            % clear croppedPeakImages2
%             warning on
        end
    end
end





