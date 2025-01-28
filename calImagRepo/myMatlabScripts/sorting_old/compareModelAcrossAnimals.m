load('C:\Users\David\Desktop\myData.mat')
load('C:\Users\David\Desktop\allValidCellMax.mat')
%load('C:\Users\David\Desktop\cellsPerBatch.mat')

dataFields = fieldnames(myData);
compactData = zeros([size(myData.(dataFields{1}),1),length(dataFields)]);
for k = 1:length(dataFields)
    tempData = myData.(dataFields{k});
    tempData(isnan(tempData)) = max(tempData);
    myData.(dataFields{k}) = tempData;
    compactData(:,k) = myData.(dataFields{k}); %No normalitzation
    %compactData(:,k) = (myData.(dataFields{k}) - min(myData.(dataFields{k})))./(max(myData.(dataFields{k})) - min(myData.(dataFields{k})));
end

trainXAll = compactData;
trainY = allValidCellMax;
clear myData allValidCellMax compactData

load('C:\Users\David\Desktop\28\myData.mat')
load('C:\Users\David\Desktop\28\allValidCellMax.mat')

compactData = zeros([size(myData.(dataFields{1}),1),length(dataFields)]);
for k = 1:length(dataFields)
    tempData = myData.(dataFields{k});
    tempData(isnan(tempData)) = max(tempData);
    myData.(dataFields{k}) = tempData;
    compactData(:,k) = myData.(dataFields{k}); %No normalitzation
    %compactData(:,k) = (myData.(dataFields{k}) - min(myData.(dataFields{k})))./(max(myData.(dataFields{k})) - min(myData.(dataFields{k})));
end

testXAll = compactData;
allValidCellMax(allValidCellMax==3) = 0;
testY = allValidCellMax;
clear myData allValidCellMax compactData
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:

featSelection = {1:numel(dataFields),[1,3,10],[11,12,16]};
permuteAcrossDays = true;
models = {'TreeBagger','SVM','Naive Bayes','Discriminant Analysis','Nearest Neighbors'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%
[trainScore,testScore,trainFpos,testFpos] = deal(zeros([length(featSelection),length(models)]));

for i = 1:length(featSelection)
 
    disp([num2str(100*(i/length(featSelection))),'%'])

    trainX = trainXAll(:,featSelection{i});
    testX = testXAll(:,featSelection{i});

    featNames{i} = dataFields(featSelection{i});
    
    for j = 1:length(models)

        model = models{j};

        switch model
            case 'TreeBagger'
              mdl = TreeBagger(60,trainX,trainY','Method','classification');
            case 'SVM' 
              mdl = fitcsvm(trainX,trainY','Standardize',true,'KernelFunction','linear');
            case 'Naive Bayes'
              mdl = fitcnb(trainX,trainY');
            case 'Discriminant Analysis'
              mdl = fitcdiscr(trainX,trainY');
            case 'Nearest Neighbors'
              mdl = fitcknn(trainX,trainY');
            otherwise
              error('no model specified')
        end 
        

        [predYTrain,predYTest] = getPredict(mdl,trainX,testX);
        trainFpos(i,j) = sum(trainY == 0 & predYTrain == 1)/length(trainY);
        testFpos(i,j) = sum(testY == 0 & predYTest == 1)/length(testY);     
        trainScore(i,j) = confusion(trainY,predYTrain);
        testScore(i,j) = confusion(testY,predYTest);
        clear mdl
    end 
end
%%
dispNames = [];
for i = 1:length(featNames)
    tempDispNames = featNames{i};
    tempDispNamesConcat = [];
    for k = 1:length(tempDispNames)
        tempDispNamesConcat = [tempDispNamesConcat,tempDispNames{k},' '];
    end
    dispNames = [dispNames,{tempDispNamesConcat}];
end

figure;
if length(featSelection) > 1
    h = bar(testScore);
    hold on
    bar(testFpos);    
    legend(h,models);
    set(gca,'xticklabel',dispNames)
    fix_xticklabels();

else
    bar(testScore);
    set(gca,'xticklabel',models)   
end
ylabel('Error')