%folders = {'E:\Processing\Mouse2028 - defaultEM','E:\Processing\Mouse2028 - 125itEM','E:\Processing\Mouse2028 - 200itEM','E:\Processing\Mouse2028 - 500itEM','E:\Processing\Mouse2028'};
% names = {'m28EMit50','m28EMit125','m28EMit200','m28EMit500','m28ICA500'};
% type = {'em','em','em','em','ica'};
% movieRegexp = {'downsampleTime','downsampleTime','downsampleTime','downsampleTime','downsampleTime'};
% daysExcluded = {};%{'20150303-icx'}%{'201502281-icx','201502282-icx','201503011-icx','201503012-icx','20150303-icx','20150305-icx','20150307-icx','20150309-icx'};
folders = {'E:\Processing\newTurboreg\Mouse2028 - reg - linear','E:\Processing\Mouse2028 - 200itEM'};%
names = {'mouse2028','mouse2028'};%
type = {'em','em'};
movieRegexp ={'downsampleTime','downsampleTime'};%{
daysExcluded = {'20150228-icx','20150301-icx','201502281-icx','201502282-icx','201503011-icx','201503012-icx','20150321-icx'};%{'201502281-icx','201502282-icx','201503011-icx','201503012-icx'};%

predPerAnimalsCell = {};
for k = 1:length(folders)
    fprintf('\n')
    fprintf('\n')
    disp(['Predicting ',folders{k}])
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    allFiles = dir(folders{k});
    predPerDays = [];
    for dayFolder = {allFiles.name}
%         try
        if ~isempty(strfind(dayFolder{1},'-icx')) && 1~=sum(strcmp(dayFolder{1},daysExcluded))%Avoids only specified days
            fprintf('\n')
            disp(['Day ',dayFolder{1}])
            insideFiles = dir([folders{k},'\',dayFolder{1}]);
            for file = {insideFiles.name}
                
                if strcmp(type{k},'em')
                    if regexp(file{1},'emAnalysis\>')
                        fileName = [folders{k},'\',dayFolder{1},'\',file{1}];
                        load(fileName)
                        inputImages = permute(emAnalysisOutput.cellImages,[3,1,2]);
                        inputSignals = double(emAnalysisOutput.scaledProbability);
%                         signalPeakIdx = emAnalysisOutput.eventTimes;
                    end
                elseif strcmp(type{k},'ica')
                    if regexp(file{1},'pcaicaAnalysis\>')
                        fileName = [folders{k},'\',dayFolder{1},'\',file{1}];
                        load(fileName)                      
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

            cInf = cellInfo(inputMovie,inputImages,inputSignals,[]);
            disp('Computing predictions...')
            cInf.predictAll('showProgress',1,'tresholds',{'getOverlap','>=0.5','getCellShapePropieties.Eccentricity','<=0.9','getCellShapePropieties.Solidity','>0.8'},'eventPercent',.2);            
            validCells = cInf.validCells;
            predPerDays = [predPerDays,100*sum(validCells)/length(validCells)];
            validCellMax = validCells;
            dayName = dayFolder{1};
            allPredPerDay.(['d',dayName(1:end-4)]) = validCellMax;
%             save([fileName(1:end-4),'Sorted.mat'],'validCellMax')
%             save([folders{k},'\',dayFolder{1},'\validCells.mat'],'validCells');
            clear inputImages inputSignals inputSignals hinfo inputMovie cInf validCells
        end
%         catch
%             warning(['Error in processing day ',num2str(dayFolder{1})])
%         end
    end
    predPerAnimalsCell = [predPerAnimalsCell;predPerDays];
    allPredPerAnimal.(names{k}) = allPredPerDay;
    clear allPredPerDay
end

%%
daysInAnimal = cellfun(@length,predPerAnimalsCell);
predPerAnimal = nan([length(predPerAnimalsCell),max(daysInAnimal)]);
for k = 1:length(predPerAnimalsCell)
    predPerAnimal(k,1:daysInAnimal(k)) = predPerAnimalsCell{k};
end
figure;
ax=subplot(1,2,1);bar(predPerAnimal);xticklabels(names);ax.XTickLabelRotation=45;title('good cells per day');grid on
ax=subplot(1,2,2);bar(nanmean(predPerAnimal,2));xticklabels(names);ax.XTickLabelRotation=45;title('average good cells');grid on
