function simDataStore = initSimDataStore(size1, size2, size3, imgsStorePath)

simDataStore.traces=cell(size1, size2, size3);
simDataStore.params=cell(size1, size2, size3);
simDataStore.covMat=cell(size1, size2, size3);
simDataStore.eventTimes=cell(size1, size2, size3);
simDataStore.nCells=zeros(size1, size2, size3);
simDataStore.eventRate=zeros(size1, size2, size3);
simDataStore.meanCellDists=zeros(size1, size2, size3);
simDataStore.meanNearestDists=zeros(size1, size2, size3);
simDataStore.meanCellCorr=zeros(size1, size2, size3);

mkdir(imgsStorePath)