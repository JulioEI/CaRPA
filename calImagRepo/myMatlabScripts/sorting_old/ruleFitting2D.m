dataDir = 'C:\Users\David\Desktop\idibaps\caImagging\Data\';
animalNames = {'m2029','m2028'};

[myData,allValidCellMax] = getData(dataDir,animalNames);

[~,~,wheresNan] = removeNansAnd3FromData(myData,allValidCellMax);
maskVect = ones(size(allValidCellMax));
maskVect(wheresNan) = 0;

[logFP,logFN, logTP, logTN] = deal(nan([100,100]));

%expRatioRange = linspace(min(myData.expRatio),max(myData.expRatio),100);
firstRule = 'RMS';
firstRuleRange = linspace(min(myData.(firstRule)),max(myData.(firstRule)),100);
secondRule = 'imMovCorr';
secondRuleRange = linspace(min(myData.(secondRule)),max(myData.(secondRule)),100);

for i = 1:length(firstRuleRange)
    for j = 1:length(secondRuleRange)
        pred = maskVect' & myData.expRatio < 1 & myData.(firstRule) < firstRuleRange(i) & myData.(secondRule) > secondRuleRange(j);

        logFP(i,j) = sum(pred'==1 & allValidCellMax==0)/sum(allValidCellMax==0)*100;
        logFN(i,j) = sum(pred'==0 & allValidCellMax==1)/sum(allValidCellMax==1)*100;
        logTP(i,j) = sum(pred'==1 & allValidCellMax==1)/sum(allValidCellMax==1)*100;
        logTN(i,j) = sum(pred'==0 & allValidCellMax==0)/sum(allValidCellMax==0)*100;
    end
end

logs = cat(3,logFP,logFN,logTP,logTN);
names = {'False positive','False negative','True positive','True negative'};

figure;
for k = 1:4
    subplot(2,2,k)
    imagesc(firstRuleRange,secondRuleRange,logs(:,:,k))
    caxis([0 100])
    colorbar
    title(names{k})
    xlabel(firstRule)
    ylabel(secondRule)
end
text(-4,1.25,[animalNames,' with expRatio < 1'],'FontSize',14)



