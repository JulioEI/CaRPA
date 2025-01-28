myPath = 'E:\Processing\Mouse2028 - 125itEM\201502282-icx';
load([myPath,'\','2015_02_28_p000_mouse2028_NULL000_emAnalysis'])
hinf = hdf5info([myPath,'\','2015_02_28_p000_mouse2028_NULL000_turboreg_crop_dfof_downsampleTime_1.h5']);
inputMovie = hdf5read(hinf.GroupHierarchy.Datasets);
[signalPeaks, signalPeaksArray, signalSigmas] = computeSignalPeaks(double(emAnalysisOutput.scaledProbability));
%%

cInf = cellInfo(inputMovie,permute(emAnalysisOutput.cellImages,[3,1,2]),emAnalysisOutput.scaledProbability,signalPeaksArray);
nCells = size(emAnalysisOutput.scaledProbability,1);

[predictedLog,goodLog,goodPercentLog] = deal(nan([1,nCells]));
predictedK = [];
for k = 1:nCells
[predicted,good,goodPercent] = cInf.predict(k,'eventMask',cInf.getOverlap(k)>=0.6,'sourceMask',cInf.getCellShapePropieties(k).Eccentricity<0.9 & cInf.getCellShapePropieties(k).Solidity>0.8,'eventPercent',.2,'doPlot',0);
% text(0,0,{['Ecc: ',num2str(cInf.getCellShapePropieties(k).Eccentricity)],['Sol: ',num2str(cInf.getCellShapePropieties(k).Solidity)],['gEv: ',num2str(b),' %: ', num2str(c),' accepted: ',num2str(a)]},'fontsize',14)
% set(gcf, 'Position', get(0,'Screensize'));
% pause()
% close all
predictedLog(k) = predicted;
if predicted
    predictedK = [predictedK,k];
end
goodLog(k) = good;
goodPercentLog(k) = goodPercent;
end

figure;bar([sum(predictedLog==1),sum(predictedLog==0)])
figure;histogram(goodLog)
figure;histogram(goodPercentLog)
%%
myIdx = 1:nCells;%predictedK;%setdiff(1:nCells,predictedK);
myIdx = myIdx(randperm(length(myIdx)));
%%
close all
figure;
for i = 1:length(myIdx)
    k = myIdx(i);
    cInf.plotEvents(k,'eventMask',cInf.getOverlap(k)>=0.6,'sourceMask',cInf.getCellShapePropieties(k).Eccentricity<0.9 & cInf.getCellShapePropieties(k).Solidity>0.8,'prediction',predictedLog(k),'openFig',0);
%     text(0,0,{['Ecc: ',num2str(cInf.getCellShapePropieties(k).Eccentricity)],['Sol: ',num2str(cInf.getCellShapePropieties(k).Solidity)],['gEv: ',num2str(goodLog(k)),' %: ', num2str(goodPercentLog(k)),' accepted: ',num2str(predictedLog(k))]},'fontsize',14)
% 	set(gcf, 'Position', get(0,'Screensize'));
	pause()
%     clf
% 	close all    
end