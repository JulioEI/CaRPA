function simDataStore = genJitterDatasetFromJitterVector(nCellsVec,jitterVec,eventRate,nTrials,imgSize,minDist,totTime,imgsStorePath)

% init result storage
simDataStore = initSimDataStore(length(nCellsVec),1,nTrials,imgsStorePath);


% start simulation loops
tic
for nCellsInd=1:length(nCellsVec)
    nCellsDesired=nCellsVec(nCellsInd);
   
    disp(['Generating simulated data.... | nCells progress ' num2str(nCellsInd/length(nCellsVec))])
    toc

    for trInd=1:nTrials

        disp(['Trial num ' num2str(trInd)])

        %%%%%%%%%%%%%%%%%%%% simulate data and store simulation results
        simOptions.jitterVec=jitterVec;
        simOptions.maxJitter=1;
        [imgs,realTraces,realParams,realEventTimes]=simulateData_makeTestData_v2(nCellsDesired,...
            imgSize,minDist,eventRate,totTime,simOptions);
        simDataStore = storeSimData(simDataStore,nCellsInd,1,trInd,...
            realTraces,realParams,imgs,realEventTimes,eventRate,[],0,imgsStorePath);

    end
end