function simDataStore = genJitterDataset(nCellsVec,jitterVec,eventRate,nTrials,imgSize,minDist,totTime,imgsStorePath)

% init result storage
simDataStore = initSimDataStore(length(nCellsVec),length(jitterVec),nTrials,imgsStorePath);


% start simulation loops
tic
for nCellsInd=1:length(nCellsVec)
    nCellsDesired=nCellsVec(nCellsInd);
    
    for jitterInd=1:length(jitterVec)
        
        maxJitter=jitterVec(jitterInd);
        disp(['Generating simulated data.... | nCells progress ' num2str(nCellsInd/length(nCellsVec))...
            ' | jitter progress ' num2str(jitterInd/length(jitterVec))])
        toc
    
        for trInd=1:nTrials

            disp(['Trial num ' num2str(trInd)])
            
            %%%%%%%%%%%%%%%%%%%% simulate data and store simulation results
            simOptions.maxJitter=maxJitter;
            [imgs,realTraces,realParams,realEventTimes]=simulateData_makeTestData_v2(nCellsDesired,...
                imgSize,minDist,eventRate,totTime,simOptions);
            simDataStore = storeSimData(simDataStore,nCellsInd,jitterInd,trInd,...
                realTraces,realParams,imgs,realEventTimes,eventRate,[],0,imgsStorePath);
            
        end
    end
end