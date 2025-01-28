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

numMethods = 5;
disp('----------------TreeBagger----------------')
mdlT = TreeBagger(60,trainX,trainY','Method','classification');
%mdlT = fitensemble(trainX,trainY','LSBoost',100,'Tree');
[treePredYTrain,treePredYTest] = getPredict(mdlT,trainX,testX);
[treeTrain, treeTest] = getConfMat(trainY,treePredYTrain,testY,treePredYTest);

% plotconfusion(testY((b(:,1)>conf | b(:,2)>conf)),a((b(:,1)>conf | b(:,2)>conf))',[num2str(conf*100),'% confidence'])
disp('----------------SVM----------------')
mdlS = fitcsvm(trainX,trainY','Standardize',true,'KernelFunction','linear');
[svmPredYTrain,svmPredYTest] = getPredict(mdlS,trainX,testX);
[svmTrain, svmTest] = getConfMat(trainY,svmPredYTrain,testY,svmPredYTest);

disp('----------------Naive Bayes----------------')
mdlNB = fitcnb(trainX,trainY');
[nbPredYTrain,nbPredYTest] = getPredict(mdlNB,trainX,testX);
[nbTrain, nbTest] = getConfMat(trainY,nbPredYTrain,testY,nbPredYTest);

disp('----------------Discriminant Analysis----------------')
mdlDA = fitcdiscr(trainX,trainY');
[daPredYTrain,daPredYTest] = getPredict(mdlDA,trainX,testX);
[daTrain, daTest] = getConfMat(trainY,daPredYTrain,testY,daPredYTest);

disp('----------------Nearest Neighbors----------------')
mdlNEI = fitcknn(trainX,trainY');
[neiPredYTrain,neiPredYTest] = getPredict(mdlNEI,trainX,testX);
[neiTrain, neiTest] = getConfMat(trainY,neiPredYTrain,testY,neiPredYTest);


%----------------Neural Network----------------
%net = fitnet(10);
%net = train(net,trainX',trainY);
%[netTrain, netTest] = mdl2ConfMat(mdlNET,trainX,trainY,testX,testY);

figure;
x = {'TreeBagger','SVM','Naive Bayes','Discriminant Analysis','Nearest Neighbors'};
y = [[treeTrain, treeTest];[svmTrain, svmTest];[nbTrain, nbTest];[daTrain, daTest];[neiTrain, neiTest]];
h = bar(y);
legend(h,{'Train error','Test error'});
set(gca,'xticklabel',x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRUE POSITIVE  FALSE POSITIVE %
% FALSE NEGATIVE TRUE NEGATIVE  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


