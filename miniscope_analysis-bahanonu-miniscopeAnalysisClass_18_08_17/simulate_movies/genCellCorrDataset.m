function simDataStore = genCellCorrDataset(meanCorrVec,eventRate,nCellsDesired,nTrials,imgSize,minDist,totTime,imgsStorePath)

% init result storage
simDataStore = initSimDataStore(length(meanCorrVec),nTrials,1,imgsStorePath);

% set simulation options
options.extraBorder=6;

% start simulation loops
tic
for meanCorrInd=1:length(meanCorrVec)

    % set mean cell correlation and display progress
    meanCorr=meanCorrVec(meanCorrInd);
    disp(['mean corr. progress ' num2str(meanCorrInd/length(meanCorrVec))])
    toc

    for trInd=1:nTrials

        % display progress
        disp(['Trial num ' num2str(trInd)])

        %%%%%%%%%%%%%%%%%%%% simulate data and store simulation results
        [imgs,realTraces,realParams, realEventTimes, C]=simulateData_cellCellCorr(nCellsDesired,...
            imgSize,minDist,eventRate,totTime,meanCorr);
        simDataStore = storeSimData(simDataStore,meanCorrInd,trInd,1,...
            realTraces,realParams,imgs,realEventTimes,eventRate,[],meanCorr,imgsStorePath);
    end
end