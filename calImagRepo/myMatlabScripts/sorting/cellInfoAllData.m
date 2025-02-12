%CellInfoAllData
%clear all

folders = {'E:\Processing\Mouse2028'}%,'E:\Processing\Mouse2028 - 125itEM','E:\Processing\Mouse2028 - 200itEM','E:\Processing\Mouse2028 - 500itEM','E:\Processing\Mouse2028'};
names = {'m28EM'}%,'m28EMit125','m28EMit200','m28EMit500','m28ICA500'};
type = {'em'}%,'em','em','em','ica'};
movieRegexp = {'downsampled'}%,'downsampleTime','downsampleTime','downsampleTime','downsampleTime'};
daysExcluded = {'Mouse-2028-20150228-linear-track','Mouse-2028-20150301-linear-track'}%{'201502281-icx','201502282-icx','201503011-icx','201503012-icx','20150303-icx','20150305-icx','20150307-icx','20150309-icx'};
timePerAnimalsCell = {};
predPerAnimalsCell = {};
for k = 1:length(folders)
    fprintf('\n')
    fprintf('\n')
    disp(['Processing ',names{k}])
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    allFiles = dir(folders{k});
    
    timePerDays = [];
    predNPerDays = [];
    for dayFolder = {allFiles.name}
        if ~isempty(strfind(dayFolder{1},'-linear-track')) && 1~=sum(strcmp(dayFolder{1},daysExcluded))%Avoids only specified days
            fprintf('\n')
            disp(['Day ',dayFolder{1}])
            insideFiles = dir([folders{k},'\',dayFolder{1}]);
            for file = {insideFiles.name}
                
                if strcmp(type{k},'em')
                    if regexp(file{1},'emAnalysis\>')
                        fileName = [folders{k},'\',dayFolder{1},'\',file{1}];
                        load(fileName)
                        timePerDays = [timePerDays,emOptions.time.endTime];
                        inputImages = permute(emAnalysisOutput.cellImages,[3,1,2]);
                        inputSignals = double(emAnalysisOutput.scaledProbability);
%                         signalPeakIdx = emAnalysisOutput.eventTimes;
                    end
                elseif strcmp(type{k},'ica')
                    if regexp(file{1},'pcaicaAnalysis\>')
                        fileName = [folders{k},'\',dayFolder{1},'\',file{1}];
                        load(fileName)
                        timePerDays = [timePerDays,pcaicaAnalysisOutput.time.endTime];
                        inputImages = permute(pcaicaAnalysisOutput.IcaFilters,[3,1,2]);
                        inputSignals = double(pcaicaAnalysisOutput.IcaTraces);
%                         [~, signalPeakIdx] = computeSignalPeaks(inputSignals,'makePlots', 0,'makeSummaryPlots',0,'waitbarOn',0);
                    end                
                end
                
                if regexp(file{1},movieRegexp{k})
                    hinfo = hdf5info([folders{k},'\',dayFolder{1},'\',file{1}]);
                    inputMovie = hdf5read(hinfo.GroupHierarchy.Datasets);
                end
            end
            
            disp('Computing data...')
            cInf = cellInfo(inputMovie,inputImages,inputSignals,[]);
            dayName = dayFolder{1};
%             dataPerDay.(['d',dayName(1:end-4)]) = cInf.runAllCells('showProgress',1);
            dayNum = regexp(dayName,'(201\d\d*)','tokens');
            dataPerDay.(['d',dayNum{1}{1}]) = cInf.runAllCells('showProgress',1);
            
            disp('Computing predictions...')
%             cInf.predictAll('showProgress',1,'tresholds',{'getOverlap','>=0.6','getCellShapePropieties.Eccentricity','<=0.9','getCellShapePropieties.Solidity','>0.8'},'eventPercent',.2);            
            cInf.setTresholds('showProgress',0,'tresholds',{'getOverlap','>0.5','getglobalSNR','>1','getCellShapePropieties.Eccentricity','<=0.9','getCellShapePropieties.Area','>100','getsScore','<0.02'});
            validCells = cInf.validCells;
            predNPerDays = [predNPerDays,100*sum(validCells)/length(validCells)];
%             predPerDay.(['d',dayName(1:end-4)]) = find(validCells);
            predPerDay.(['d',dayNum{1}{1}]) = validCells;
%             switch type{k}
%                 case 'em'
%                     validCellMax = logical(validCells);
%                     save([fileName(1:end-4),'Sorted.mat'],'validCellMax')                   
%                 case 'ica'
%                     valid = logical(validCells);
%                     save([fileName(1:end-18),'ICdecisions'],'valid') 
%                 otherwise
%                     error('type not understood')
%             end
            
            clear inputImages inputSignals inputSignals signalPeakIdx hinfo inputMovie cInf validCells
        end
    end
    timePerAnimalsCell = [timePerAnimalsCell;timePerDays];
    predPerAnimalsCell = [predPerAnimalsCell;predNPerDays];
    dataPerAnimal.(names{k}) = dataPerDay;
    predPerAnimalAll.(names{k}) = predPerDay;
    clear dataPerDay predPerDay
end
%%
daysInAnimal = cellfun(@length,timePerAnimalsCell);
timePerAnimal = nan([length(timePerAnimalsCell),max(daysInAnimal)]);
for k = 1:length(timePerAnimalsCell)
    timePerAnimal(k,1:daysInAnimal(k)) = timePerAnimalsCell{k};
end
figure;
ax=subplot(1,2,1);bar(timePerAnimal);xticklabels(names);ax.XTickLabelRotation=45;title('total time per day');grid on
ax=subplot(1,2,2);bar(nanmean(timePerAnimal,2));xticklabels(names);ax.XTickLabelRotation=45;title('average time per day');grid on

predPerAnimal = nan([length(predPerAnimalsCell),max(daysInAnimal)]);
for k = 1:length(predPerAnimalsCell)
    predPerAnimal(k,1:daysInAnimal(k)) = predPerAnimalsCell{k};
end
figure;
ax=subplot(1,2,1);bar(predPerAnimal);xticklabels(names);ax.XTickLabelRotation=45;title('good cells per day');grid on
ax=subplot(1,2,2);bar(nanmean(predPerAnimal,2));xticklabels(names);ax.XTickLabelRotation=45;title('average good cells');grid on
