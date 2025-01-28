% dataDir = 'C:\Users\David\Desktop\idibaps\caImagging\Data\';
% animalNames = {'m2029','m2028'};
% 
% [myData,allValidCellMax] = getData(dataDir,animalNames);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Keep data inside rule boundaries
%[goodDataLgc,myData] = rulePrd(myData);
%allValidCellMax = allValidCellMax(goodDataLgc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[logFP,logFN, logTP, logTN] = deal(nan([100,1]));

figure;

firstRule = 'imMovCorr';
firstRuleRange = linspace(0,1,100);

for k = 1:length(firstRuleRange)
    pred = myData.(firstRule) > firstRuleRange(k) & rulePrd(myData);
    logFP(k) = sum(pred'==1 & allValidCellMax==0)/sum(allValidCellMax==0)*100;
    logFN(k) = sum(pred'==0 & allValidCellMax==1)/sum(allValidCellMax==1)*100;
    logTP(k) = sum(pred'==1 & allValidCellMax==1)/sum(allValidCellMax==1)*100;
    logTN(k) = sum(pred'==0 & allValidCellMax==0)/sum(allValidCellMax==0)*100;
end

plot(firstRuleRange,logFP,'linewidth',2)
hold on
plot(firstRuleRange,logFN,'linewidth',2)
plot(firstRuleRange,logTP,'linewidth',2)
plot(firstRuleRange,logTN,'linewidth',2)
hold off
xlabel([firstRule])
ylabel('%')
% title([animalNames])%,' with RMS > 0.5 and expRatio < 1'])
axis square
grid on
legend('False positive','False negative','True positive','True negative')

%%
%Filter mydata
myData = myDataSave;
allValidCellMax = allValidCellMaxSave;
predBasic = rulePrd(myData);
for field = fields(myData)'
    myData.(field{1}) = myData.(field{1})(predBasic);
end
allValidCellMax = allValidCellMax(predBasic);


