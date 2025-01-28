function simDataStore = genUncorrDataset(nCellsVec,eventRateVec,nTrials,imgSize,minDist,totTime,imgsStorePath)


% init result storage
simDataStore = initSimDataStore(length(nCellsVec),length(eventRateVec),nTrials,imgsStorePath);

% set simulation options
options.extraBorder=6;

% start simulation loops
tic
for nCellsInd=1:length(nCellsVec)
    nCellsDesired=nCellsVec(nCellsInd);
    
    for evRateInd=1:length(eventRateVec)
        
        eventRate=eventRateVec(evRateInd);
        disp(['nCells progress ' num2str(nCellsInd/length(nCellsVec))])
        disp(['spike rate progress ' num2str(evRateInd/length(eventRateVec))])
        toc
    
        for trInd=1:nTrials

            %%%%%%%%%%%%%%%%%%%% simulate data and store simulation results
            [imgs,realTraces,realParams,realEventTimes]=simulateData_makeTestData_v2(nCellsDesired,...
                imgSize,minDist,eventRate,totTime,options);
            simDataStore = storeSimData(simDataStore,nCellsInd,evRateInd,trInd,...
                realTraces,realParams,imgs,realEventTimes,eventRate,[],0,imgsStorePath);
        end
    end
    
    
end
