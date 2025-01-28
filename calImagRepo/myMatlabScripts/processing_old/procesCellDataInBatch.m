%Process all cell data

days = {'282','281','31','29','27','25'};%{'06','11','16','21','26','31','281','282','011','012'};
myDir = 'C:\Users\David\Desktop\temp_CellINfo\';

allFiles = dir(myDir);
cellsPerBatch = [];
allValidCellMax = [];
myData = [];
for day = days
    disp(['Loading day ',day{1}])
    regList = regexp({allFiles(:).name},['^201\d\d\d',day{1},'-icx$']);

    
    indexAll = find(~cellfun(@isempty,regList));
    
    if length(indexAll) == 1
        loadFrom = [myDir,allFiles(indexAll).name,'\cellInfo'];
        variableNames = {'inputSignals','testpeaksArray','cellImages','validCellMax','inputMovie'};
        for variable = variableNames
            load([loadFrom,'\',variable{1},'.mat'])
        end
        cellsPerBatch = [cellsPerBatch,size(inputSignals,1)];
        allValidCellMax = [allValidCellMax,validCellMax];
       % myTempData = getCellInfo(inputSignals,testpeaksArray,permute(cellImages,[2,3,1]),inputMovie);
        disp([num2str(sum(cellfun(@isempty,testpeaksArray))/length(testpeaksArray)*100),'% cells with no events'])
        myTempData = simpleGetCellInfo(inputSignals,testpeaksArray,permute(cellImages,[2,3,1]),inputMovie);
        if isempty(myData)
            myData = myTempData;
        else
            dataFields = fieldnames(myData);
            for i = 1:length(dataFields)
                myData.(dataFields{i}) = [myData.(dataFields{i});myTempData.(dataFields{i})];
            end
        end
        clear(variableNames{:});
        clear myTempData
        
    elseif length(indexAll) > 1
        error('Regexp returned two files')
    end
        
end
%%
% myFeatures = {'expRatio','expPeak','globalSNR','peakSNR','RMS'};
% myNames = {'Ratio of exponentials','Event peak','SNR global','SNR peak variability','RMS'};
% feature2DScatter(myFeatures,myData,allValidCellMax,myNames)
% myFeatures2 = {'cellEcc','cellSol','pScore','imMovCorr'};
% myNames2 = {'Cell ecentricity','Cell solidity','Ratio between perimeter and perimeter of ellipsis','Event correlation'};
% feature2DScatter(myFeatures2,myData,allValidCellMax,myNames2)