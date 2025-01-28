movieFolders = {'C:\Users\csc\Desktop\caImagg\Behavior\Linear-track\Mouse-2028\Behavioral-videos'};%;{'C:\Users\csc\Desktop\caImagg\Behavior\Linear-track\Mouse-2028\Behavioral-videos'};%
folders = {'E:\Processing\Mouse2028 - 200itEM'};%{'E:\Processing\newTurboreg\Mouse2028 - reg - linear'};
names = {'mouse2028'}; %{'mouse2028'};%
type = {'em'};%{'em'};
movieRegexp ={'downsampleTime'};% {'downsampleTime'};%{
daysExcluded = {'201502281-icx','201502282-icx','201503011-icx','201503012-icx','20150321-icx'};%{'20150228-icx','20150301-icx'};%{'201502281-icx','201502282-icx','201503011-icx','201503012-icx'};%
savePredictions = false;
errorPerAnimalsCell = {};
errorPerAnimalsFromFirstCell = {};
for k = 1:length(folders)
    fprintf('\n')
    fprintf('\n')
    disp(['Processing ',names{k}])
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    allFiles = dir(folders{k});
    
    load([folders{k},'\globalIDs.mat'])
            
    errorPerDays = [];
    errorPerDaysFromFirst = [];
    dayN = 1;
    for dayFolder = {allFiles.name}
        if ~isempty(strfind(dayFolder{1},'-icx')) && 1~=sum(strcmp(dayFolder{1},daysExcluded))%Avoids only specified days
            fprintf('\n')
            disp(['Day ',dayFolder{1}])  
            insideFiles = dir([folders{k},'\',dayFolder{1}]);
            doPrediction = 1;
            for file = {insideFiles.name}
                if strcmp(type{k},'em')
                    if regexp(file{1},'emAnalysis\>')
                        fileName = [folders{k},'\',dayFolder{1},'\',file{1}];
                        load(fileName)
                        inputImages = emAnalysisOutput.cellImages;
                        inputSignals = double(emAnalysisOutput.scaledProbability);
                    end
                    if regexp(file{1},'emAnalysisSorted\>')
                        doPrediction = 0;
                        sortedName = [folders{k},'\',dayFolder{1},'\',file{1}];
                        load(sortedName)
                    end
                elseif strcmp(type{k},'ica')
                    if regexp(file{1},'pcaicaAnalysis\>')
                        error('not implemented')
                    end                
                end
                
                if regexp(file{1},movieRegexp{k})
                    hinfo = hdf5info([folders{k},'\',dayFolder{1},'\',file{1}]);
                    inputMovie5hz = hdf5read(hinfo.GroupHierarchy.Datasets);
                end
                if regexp(file{1},'dfof_1.h5')
                    hinfo = hdf5info([folders{k},'\',dayFolder{1},'\',file{1}]);
                    inputMovie20hz = hdf5read(hinfo.GroupHierarchy.Datasets);
                end
            end

            %Get the behavioral movie traces
            disp('Computing behavior traces...')
            movieDir = dir(movieFolders{k});
            movieSearchName = dayFolder{1}(1:end-4);
            if length(movieSearchName) == 9
                movieSearchName = dayFolder{1}(1:end-5);
                movieSession = str2double(dayFolder{1}(end-4));
            elseif length(movieSearchName) == 10
                movieSearchName = dayFolder{1}(1:end-6);
                movieSession = str2double(dayFolder{1}(end-5:end-4));               
            end
            movieFileIdx = find(cellfun(@(x) ~isempty(x),strfind({movieDir.name},movieSearchName)));
            
            if length(movieFileIdx)>1
                movieFileIdx = movieFileIdx(movieSession); %This supposes time ordering!
            end
            moviePath = [movieFolders{k},'\',movieDir(movieFileIdx).name];
            [position,~] = getMouseTrajectory(moviePath);
            
            %Get the predictions
            if doPrediction
            warning('Predictions not found');
            cInf = cellInfo(inputMovie5hz,permute(inputImages,[3,1,2]),inputSignals,[]);
            disp('Computing predictions...')
            cInf.predictAll('showProgress',1,'tresholds',{'getOverlap','>=0.6','getCellShapePropieties.Eccentricity','<=0.9','getCellShapePropieties.Solidity','>0.8'},'eventPercent',.2);            
            validCellMax = cInf.validCells;
                if savePredictions
                    save([fileName(1:end-4),'Sorted.mat'],'validCellMax')
                end
            end
            
            %Compute upscaled traces
            try
                load([folders{k},'\',dayFolder{1},'\upScaledPhi'])
                if length(validCellMax) == size(upScaledPhi,1)
                    upScaledPhi = upScaledPhi(logical(validCellMax),:);
                end
            catch
                disp('Computing upscaled traces...')
                robustMovie = inputMovie20hz;
                robustMovie(isnan(inputMovie20hz))=min(inputMovie20hz(:));
                [upScaledPhi, upFilteredTraces] = recalcPhiAndDetectEvents(robustMovie,inputImages(:,:,logical(validCellMax)),emAnalysisOutput.dsCellTraces(logical(validCellMax),:),emAnalysisOutput.CELLMaxoptions,'runEventDetection',0);
                upScaledPhi = upScaledPhi*10^5;
                save([folders{k},'\',dayFolder{1},'\upFilteredTraces'],'upFilteredTraces')
                save([folders{k},'\',dayFolder{1},'\upScaledPhi'],'upScaledPhi')
            end
            
            normalizedSP = upScaledPhi./repmat(max(upScaledPhi')',[1,size(upScaledPhi,2)]);
            x = normalizedSP';
            
            %Decode position
            normalizedPosX = (position(:,1)-min(position(:,1)))/(max(position(:,1))-min(position(:,1)));
            binN = 20;
            yTemp = zeros([1,length(normalizedPosX)]);
            for j = 1:binN
                if j ~= binN
                    idx = find(normalizedPosX>=(j-1)/binN&normalizedPosX<(j/binN));
                else
                    idx = find(normalizedPosX>=(j-1)/binN&normalizedPosX<=(j/binN));
                end
                yTemp(idx) = j;
            end
            noPopulationBins = find(hist(yTemp,linspace(1,20,20)) < sum(hist(yTemp,linspace(1,20,20))/20/5)); %Find bins where population is 1/5 from uniform and recompute them
            for j = noPopulationBins
                if j ~= binN
                    idx = find(normalizedPosX>=(j-1)/binN&normalizedPosX<(j/binN));
                else
                    idx = find(normalizedPosX>=(j-1)/binN&normalizedPosX<=(j/binN));
                end
                normalizedPosX(idx) = [];
                x(idx,:) = [];
            end
            
            normalizedPosX = (normalizedPosX-min(normalizedPosX))/(max(normalizedPosX)-min(normalizedPosX));
            y = zeros([1,length(normalizedPosX)]);
            for j = 1:binN
                if j ~= binN
                    idx = find(normalizedPosX>=(j-1)/binN&normalizedPosX<(j/binN));
                else
                    idx = find(normalizedPosX>=(j-1)/binN&normalizedPosX<=(j/binN));
                end
                y(idx) = j;
            end
                       
            if size(x,1) ~= length(y)
                disp('Behavior and calcium have different lengths. Manual action requiered')
                pause()
            end
            %Bootstrap with length(y)/n data points, in uniform chunks
            bootstrapError = spatialBootstrap(x,y,10);
            errorPerDays = [errorPerDays,mean(bootstrapError)];
            
            if dayN == 1
                firstDayMdl = fitcnb(x,y);
                firstDayCellN = size(x,2);
                idxFirstDay = find(globalIDs(:,1));
                in1CellWas = globalIDs(idxFirstDay,1); %CHANGE WHEN GLOBALID IS FIXED
                relativeError = 0; %CHANGE WHEN GLOBALID IS FIXED
            else
                               
            sameCellsAsFirst = globalIDs(idxFirstDay,dayN);
            xPad = zeros([length(y),firstDayCellN]);
            for j = 1:length(in1CellWas)
                if sameCellsAsFirst(j)~=0
                    xPad(:,in1CellWas(j)) = x(:,sameCellsAsFirst(j));
                end
            end
            
            relativeError = mean(sqrt((firstDayMdl.predict(xPad)-y').^2));
            end
            
            errorPerDaysFromFirst = [errorPerDaysFromFirst,relativeError];

%             dayName = dayFolder{1};            
            clear inputImages inputSignals inputSignals signalPeakIdx hinfo inputMovie cInfo
            dayN = dayN + 1;
        end
    end
    errorPerAnimalsCell = [errorPerAnimalsCell;errorPerDays];
    errorPerAnimalsFromFirstCell = [errorPerAnimalsFromFirstCell;errorPerDaysFromFirst];
%     dataPerAnimal.(names{k}) = dataPerDay;
%     clear dataPerDay
end
%%
daysInAnimal = cellfun(@length,errorPerAnimalsCell);
timePerAnimal = nan([length(errorPerAnimalsCell),max(daysInAnimal)]);
for k = 1:length(errorPerAnimalsCell)
    timePerAnimal(k,1:daysInAnimal(k)) = errorPerAnimalsCell{k};
end
figure;
ax=subplot(1,2,1);bar(timePerAnimal);xticklabels(names);ax.XTickLabelRotation=45;title('total time per day');grid on
ax=subplot(1,2,2);bar(nanmean(timePerAnimal,2));xticklabels(names);ax.XTickLabelRotation=45;title('average time per day');grid on