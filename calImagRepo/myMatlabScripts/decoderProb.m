
dataFields = fields(myData);
compactData = zeros([size(myData.(dataFields{1}),1),length(dataFields)]);
for k = 1:length(dataFields)
    compactData(:,k) = myData.(dataFields{k});
end
%%
holdoutLength = floor(.3 *length(allValidCellMax));
randVect = randperm(length(allValidCellMax));
testIdx = randVect(1:holdoutLength);
trainIdx = randVect((holdoutLength+1):end);

trainX = compactData(trainIdx,:);
trainY = allValidCellMax(trainIdx);
testX = compactData(testIdx,:);
testY = allValidCellMax(testIdx);

CVSVMModel = fitcsvm(trainX,trainY,'Standardize',true);
[ScoreCVSVMModel,transform] = fitSVMPosterior(CVSVMModel);
disp('----------Test set----------')
getConfusion(testY,CVSVMModel.predict(testX)');
disp('----------Train set----------')
getConfusion(trainY,CVSVMModel.predict(trainX)');
%%
[pred,predProb] = ScoreCVSVMModel.predict(testX);
predProb = predProb(:,2);
figure;histogram(predProb,20,'Normalization','probability');title('Distribution of prediction probabilities');axis('square')
%%
tresholdVec = .5:.001:1;
[incorrectR,falseNRatio,falsePRatio,percentCells] = deal(zeros([1,length(tresholdVec)]));
for k = 1:length(tresholdVec)
    subTreshIdx = find((predProb >= tresholdVec(k)) | (predProb <= 1-tresholdVec(k)));
    incorrectR(k) = sum(testY(subTreshIdx)' ~= pred(subTreshIdx))/length(subTreshIdx);
    falseNRatio(k) = sum(testY(subTreshIdx)' == 1 & pred(subTreshIdx) == 0)/length(subTreshIdx);
    falsePRatio(k) = sum(testY(subTreshIdx)' == 0 & pred(subTreshIdx) == 1)/length(subTreshIdx);
    percentCells(k) = length(subTreshIdx)/length(predProb);
end

figure;plot(tresholdVec,incorrectR,'lineWidth',2);hold on;
plot(tresholdVec,falseNRatio,'lineWidth',2)
plot(tresholdVec,falsePRatio,'lineWidth',2)
plot(tresholdVec,percentCells,'lineWidth',2);
xlabel('treshold');legend('error','falseN','falseP','predictedCells');ylabel('%');grid on;
title('Ratios when predicting all cells above the specified treshold probability');axis('square')
%%
cInf = cellInfo('E:\Processing\Mouse2028\Mouse-2028-20150303-linear-track\Mouse-2028-20150303_103622-linear-track-dfof-downsampled.h5',permute(emAnalysisOutput.cellImages,[3,1,2]),emAnalysisOutput.scaledProbability,[],'validCells',validCellMax);
%%
basicGoodIdx = find(predBasic);
treshold = .9;
goodTreshIdx = basicGoodIdx(testIdx(find(predProb >= treshold)))';
badTreshIdx = basicGoodIdx(testIdx(find(predProb < (1-treshold))))';
ambiguousTreshIdx = basicGoodIdx(testIdx(find(predProb < treshold & predProb >= (1-treshold))))';

fullDecisions = zeros([1,length(predBasic)]);
fullDecisions(goodTreshIdx) = 1;
fullDecisions(ambiguousTreshIdx) = 3;

% fullDecisions = zeros([1,length(predBasic)]);
% fullDecisions(basicGoodIdx) = 3;

thisDayDecisions = fullDecisions(1:length(cInf.inputImages));
cInf.validCells = thisDayDecisions;
cInf.manualSort('tresholds',{'getOverlap','>=.4','getglobalSNR','>1'})%basicGoodIdx(nonNanIdx))

