function [ myWholeData,allValidCellMaxWhole,cellsPerBatchWhole] = getData(dataDir,animalNames)
for animal = 1:length(animalNames)
    animalName = animalNames{animal};
    load([dataDir,animalName,'\myData.mat'])
    load([dataDir,animalName,'\allValidCellMax.mat'])  
    load([dataDir,animalName,'\cellsPerBatch.mat']) 
    if animal == 1
        dataFields = fields(myData);
        myWholeData = myData;
        allValidCellMaxWhole = allValidCellMax;
        cellsPerBatchWhole = cellsPerBatch;
    else
        for k = 1:length(dataFields)
            myWholeData.(dataFields{k}) = [myWholeData.(dataFields{k});myData.(dataFields{k})];
        end
        allValidCellMaxWhole = [allValidCellMaxWhole,allValidCellMax];
        cellsPerBatchWhole = [cellsPerBatchWhole,cellsPerBatch];
    end
    clear myData allValidCellMax cellsPerBatch
end
end

