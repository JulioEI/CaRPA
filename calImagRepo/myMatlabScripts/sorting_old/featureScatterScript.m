% dataDir = 'C:\Users\David\Desktop\idibaps\caImagging\Data\';
% animalNames = {'m2029','m2028'};
% 
% [myData,allValidCellMax] = getData(dataDir,animalNames);

% predConcat = [];dayFields = fields(predPerAnimalAll.m28EM);for k = 1:length(dayFields);predConcat = [predConcat,predPerAnimalAll.m28EM.(dayFields{k})];end;
% dataPerAnimalConcat = concatDataPerAnimal(dataPerAnimal,'doMean',1);
% myData = dataPerAnimalConcat.m28EM;
% allValidCellMax = predConcat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Keep data inside rule boundaries
% [goodDataLgc,myData] = rulePrd(myData);
% allValidCellMax = allValidCellMax(goodDataLgc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%For histograms


figure;
dataFields = fields(myData);
dataRange = [0,2;0,4;0,20;0,4;0,1;0.5,1;0,1.3;0,5;0,0.5;0,150;0,2;-0.1,1;0,20;-0.1,1;0,1;0,1;0.5,1;0,1.2;0,5;-0.1,0.6;0,150;0,0.5];
binN = 20;
for k = 1:length(dataFields)
    subplot(ceil(sqrt(length(dataFields))),ceil(sqrt(length(dataFields))),k)
    histData = myData.(dataFields{k});
    histDataCrop = histData(~isnan(histData));
    validCrop = allValidCellMax(~isnan(histData)); 
    %minAxis = min(histDataCrop(abs(histDataCrop) < mean(histDataCrop)+3*std(histDataCrop)));
    %maxAxis = max(histDataCrop(abs(histDataCrop) < mean(histDataCrop)+3*std(histDataCrop)));
%     histogram(histDataCrop(~logical(validCrop)),[dataRange(k,1):dataRange(k,2)/binN:dataRange(k,2)],'FaceColor',[174/255,0,45/255]);
    histogram(histDataCrop(~logical(validCrop)),'FaceColor',[174/255,0,45/255]);
    hold on
%     histogram(histDataCrop(logical(validCrop)),[dataRange(k,1):dataRange(k,2)/binN:dataRange(k,2)],'FaceColor',[121/255,180/255,0])
     histogram(histDataCrop(logical(validCrop)),'FaceColor',[121/255,180/255,0])
    hold off
%     xlim(dataRange(k,:))%([minAxis,maxAxis])
    title(dataFields{k})
end
% text(.75,1.25,animalNames,'FontSize',14)

%%
%For scatter plots
showData = true;
showBoundary = false;
showEllipsis = true;

myFeatures = {'cellSol','sScore','imMovCorr','interCorr','eventSol','eventC','eventPer'};
feature2DScatter(myFeatures,myData,allValidCellMax,showData,showBoundary,showEllipsis)

% myFeatures = {'expRatio','expPeak','globalSNR','RMS'};%,'peakSNR','RMS'};
% myNames = {'Ratio of exponentials','Event peak','SNR global','RMS'};%,'SNR peak variability','RMS'};
% feature2DScatter(myFeatures,myData,allValidCellMax,showData,showBoundary,showEllipsis,myNames)
% myFeatures2 = {'cellEcc','cellSol','pScore','imMovCorr'};
% myNames2 = {'Cell ecentricity','Cell solidity','Ratio between perimeter and perimeter of ellipsis','Event correlation'};
% feature2DScatter(myFeatures2,myData,allValidCellMax,showData,showBoundary,showEllipsis,myNames2)
