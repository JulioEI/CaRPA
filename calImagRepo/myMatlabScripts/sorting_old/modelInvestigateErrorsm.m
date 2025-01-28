%Train:
clear all
load('C:\Users\David\Desktop\myData.mat')
load('C:\Users\David\Desktop\allValidCellMax.mat')

%Compact the data into a observationsXfeatures matrix
dataFields = fieldnames(myData);
dataFields = dataFields(1:end);

cellN = size(myData.(dataFields{1}),1);

compactData = zeros([size(myData.(dataFields{1}),1),length(dataFields)]);

for k = 1:length(dataFields)
    
    tempData = myData.(dataFields{k});
    tempData(isnan(tempData)) = max(tempData);
    myData.(dataFields{k}) = tempData;
    
    compactData(:,k) = myData.(dataFields{k}); %No normalitzation
    %compactData(:,k) = (myData.(dataFields{k}) - min(myData.(dataFields{k})))./(max(myData.(dataFields{k})) - min(myData.(dataFields{k})));
end
%

ratioTrainTest = 0.7;
randIndx = randperm(size(compactData,1));

elemTrain = ceil(ratioTrainTest*size(randIndx,2));
trainX = compactData(randIndx(1:elemTrain),:);
trainY = allValidCellMax(randIndx(1:elemTrain));
testX = compactData(randIndx(elemTrain+1:end),:); 
testY = allValidCellMax(randIndx(elemTrain+1:end));

%mdl = TreeBagger(60,trainX,trainY','Method','classification');
mdl = fitcsvm(trainX,trainY','Standardize',true,'KernelFunction','linear');
[predYTrain,predYTest] = getPredict(mdl,trainX,testX);

falsePositives = find(predYTest == 1 & testY == 0);
falseNegatives = find(predYTest == 0 & testY == 1);

idx = [1,3,12];%featureIdxSortbyP(1:3);

figure;
plot3(compactData(falsePositives,idx(1)),compactData(falsePositives,idx(2)),compactData(falsePositives,idx(3)),'.','color','k')
hold on
plot3(compactData(falseNegatives,idx(1)),compactData(falseNegatives,idx(2)),compactData(falseNegatives,idx(3)),'.','color','c')

treeCorrectGood = predYTest == testY & testY == 1;
treeCorrectBad = predYTest == testY & testY == 0;
plot3(compactData(treeCorrectGood,idx(1)),compactData(treeCorrectGood,idx(2)),compactData(treeCorrectGood,idx(3)),'.','color','g')
plot3(compactData(treeCorrectBad,idx(1)),compactData(treeCorrectBad,idx(2)),compactData(treeCorrectBad,idx(3)),'.','color','r')


xlabel(dataFields{idx(1)})
ylabel(dataFields{idx(2)})
zlabel(dataFields{idx(3)})
legend('FalsePositives','FalseNegatives','CorrectGood','CorrectBad')


