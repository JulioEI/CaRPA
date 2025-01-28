dataDir = 'C:\Users\David\Desktop\idibaps\caImagging\Data\';
animalNames = {'m2029','m2028'};

[myData,allValidCellMax] = getData(dataDir,animalNames);
fieldNames = fields(myData);

[logFP,logFN, logTP, logTN] = deal(nan([100,1]));

maskableFields = {'expRatio','expPeak','cellSol','pScore','cScore','sScore','cellPer','RMS','interCorr','centroidDist','eventSol','eventP','eventC','eventS','eventPer','imMovCorr'};

firstRule = 'strength';
firstRuleRange = linspace(0,2,100);

tpSum = nan([1,length(maskableFields)]);
[~,wheres1] = min(abs(firstRuleRange-1));

figure;
for j = 1:length(maskableFields)

maskedFields = fieldNames(~strcmp(fieldNames,maskableFields{j}));

for k = 1:length(fieldNames)
   if sum(strcmp(fieldNames{k},maskedFields))
        strengthMask.(fieldNames{k}) = 0;
   else
       strengthMask.(fieldNames{k}) = 1;
   end
end

for k = 1:length(firstRuleRange)
    strength = firstRuleRange(k);
    pred = rulePrd(myData,strengthMask,strength);
   
    logFP(k) = sum(pred'==1 & allValidCellMax==0)/sum(allValidCellMax==0)*100;
    logFN(k) = sum(pred'==0 & allValidCellMax==1)/sum(allValidCellMax==1)*100;
    logTP(k) = sum(pred'==1 & allValidCellMax==1)/sum(allValidCellMax==1)*100;
    logTN(k) = sum(pred'==0 & allValidCellMax==0)/sum(allValidCellMax==0)*100;
    if k == wheres1;
        tpSum(j) = 100-logTP(k);
    end
end

subplot(ceil(sqrt(length(maskableFields))),ceil(sqrt(length(maskableFields))),j)

plot(firstRuleRange,logFP,'linewidth',2)
hold on
plot(firstRuleRange,logFN,'linewidth',2)
plot(firstRuleRange,logTP,'linewidth',2)
plot(firstRuleRange,logTN,'linewidth',2)
hold off
xlabel([firstRule])
ylabel('%')
title([maskableFields{j}])
%title([animalNames])%,' with RMS > 0.5 and expRatio < 1'])
axis square
grid on
end
%%
legend('False positive','False negative','True positive','True negative')
text(-12.75,1.25,animalNames,'FontSize',14)
figure;stem(tpSum)
title({'% of true positives discarded at strenght == 1.',['Worst case scenario: ',num2str(sum(tpSum)),'% TP loss']})
xlim([0,length(maskableFields)+1])
ax = gca;
ax.XTick = 1:length(maskableFields);
ax.XTickLabel = (maskableFields);
ylabel('% error')


