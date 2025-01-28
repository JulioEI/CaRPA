function simDataStore = storeSimData(simDataStore,ind1,ind2,ind3,realTraces,realParams,imgs,eventTimes,eventRate,C,meanCellCorr,imgsStorePath) %#ok<INUSL>

simDataStore.nCells(ind1,ind2,ind3)=size(realParams,1);
simDataStore.eventRate(ind1,ind2,ind3)=eventRate;
simDataStore.traces{ind1,ind2,ind3}=realTraces;
simDataStore.params{ind1,ind2,ind3}=realParams;
simDataStore.eventTimes{ind1,ind2,ind3}=eventTimes;
simDataStore.covMat{ind1,ind2,ind3}=C;
simDataStore.meanCellCorr(ind1,ind2,ind3)=meanCellCorr;
[simDataStore.meanCellDists(ind1,ind2,ind3),...
    simDataStore.meanNearestDists(ind1,ind2,ind3),~]=getCellDistances(realParams);

save([imgsStorePath filesep 'imgs_' num2str(ind1) '_' num2str(ind2) '_' num2str(ind3)],'imgs')