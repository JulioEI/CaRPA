%Train:
[myData,allValidCellMax] = removeNansAnd3FromData(myData,allValidCellMax);
%Compact the data into a observationsXfeatures matrix
dataFields = fieldnames(myData);
cellN = size(myData.(dataFields{1}),1);

compactData = zeros([size(myData.(dataFields{1}),1),length(dataFields)]);

for k = 1:length(dataFields)
    
    %tempData = myData.(dataFields{k});
    %tempData(isnan(tempData) | isinf(tempData)) = max(tempData(~isinf(tempData)));
    %myData.(dataFields{k}) = tempData;
    
    compactData(:,k) = myData.(dataFields{k}); %No normalitzation
    %compactData(:,k) = (myData.(dataFields{k}) - min(myData.(dataFields{k})))./(max(myData.(dataFields{k})) - min(myData.(dataFields{k})));
end

%%
%Plot logistic regression boundary on the 3 best features.
normalizedData = compactData;%(compactData - min(min(compactData)))./(max(max((compactData))) - min(min(compactData)));
B = mnrfit(normalizedData,categorical(allValidCellMax));
[~,idx] = sort(abs(B(2:end)));
% idx = idx(end:-1:1);
sortedDF = dataFields(idx);
for i = 1:length(dataFields)
    myNormalizedData.(dataFields{i}) = normalizedData(:,i);
end
feature3DScatter(sortedDF(end-2:end),myNormalizedData,allValidCellMax,true)
hold on
B2 = mnrfit(normalizedData(:,idx(end-2:end)),categorical(allValidCellMax));
%[X,Y] = meshgrid(linspace(min(normalizedData(:,idx(end-2))),max(normalizedData(:,idx(end-2))),100),linspace(min(normalizedData(:,idx(end-1))),max(normalizedData(:,idx(end-1))),100));
%Z = (B2(1)+B2(2).*X+B2(3).*Y)./(-B2(4));
%mesh(X,Y,Z)

%%
%Get the confusion matrix of those features applyed to the training data.
wrongCells = zeros([1,size(normalizedData,1)]);
predLog = NaN([1,size(normalizedData,1)]);
figure;
for k = 1:size(normalizedData,1)
    val = B2(1) + B2(2)*normalizedData(k,idx(end-2)) + B2(3)*normalizedData(k,idx(end-1)) + B2(4)*normalizedData(k,idx(end));
    if ~isnan(val)
        if val<0;
            predicted = 1;
            predLog(k) = 1;
        else
            predicted = 0;
            predLog(k) = 0;
        end
        if allValidCellMax(k) == predicted
            plot3(normalizedData(k,idx(end-2)),normalizedData(k,idx(end-1)),normalizedData(k,idx(end)),'g.')
        else
            wrongCells(k) = k;
            plot3(normalizedData(k,idx(end-2)),normalizedData(k,idx(end-1)),normalizedData(k,idx(end)),'r.')
        end
    end
    hold on
end

targets = allValidCellMax;%(~isnan(predLog));
outputs = predLog;%(~isnan(predLog));

[confCoef,confMat,confInd] = confusion(targets,outputs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRUE POSITIVE  FALSE POSITIVE %
% FALSE NEGATIVE TRUE NEGATIVE  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%Sort the features by their p value
Y = categorical(allValidCellMax);
X = compactData;
holdoutCVP = cvpartition(Y,'holdout',1000);
dataTrain = X(holdoutCVP.training,:);
grpTrain = Y(holdoutCVP.training);
dataTrainG1 = dataTrain(grp2idx(grpTrain)==1,:);
dataTrainG2 = dataTrain(grp2idx(grpTrain)==2,:);
[h,p,ci,stat] = ttest2(dataTrainG1,dataTrainG2,'Vartype','unequal');
%ecdf(p);
%xlabel('P value');
%ylabel('CDF value')
[~,featureIdxSortbyP] = sort(p,2); % sort the features

feature3DScatter(dataFields(featureIdxSortbyP(1:3)),myData,allValidCellMax,true)

%%
%Sort features using random forest
idx = validCellMax;%kmeans(compactData,2);

Mdl = TreeBagger(60,compactData,idx','OOBPrediction','On',...
    'Method','classification');
% figure;
% oobErrorBaggedEnsemble = oobError(Mdl);
% plot(oobErrorBaggedEnsemble)
% xlabel 'Number of grown trees';
% ylabel 'Out-of-bag classification error';

%SVMModel = fitcsvm(compactData,idx','Standardize',true,'KernelFunction','linear');
%%
% Lambda = logspace(-6,-0.5,11);
% CVMdl = fitclinear(compactData,idx','KFold',5,...
%     'Learner','logistic','Solver','sparsa','Regularization','lasso',...
%     'Lambda',Lambda,'GradientTolerance',1e-8);
% ce = kfoldLoss(CVMdl);
% 
% Mdl = fitclinear(compactData,idx',...
%     'Learner','logistic','Solver','sparsa','Regularization','lasso',...
%     'Lambda',Lambda,'GradientTolerance',1e-8);
% numNZCoeff = sum(Mdl.Beta~=0);
% MdlFinal = selectModels(Mdl,9);

%%
%CVSVMModel = crossval(SVMModel);
%classLoss = kfoldLoss(CVSVMModel)

% myLineW = 1;
% mdlPredict = Mdl.predict(compactData);
% mdlPredict = str2double({mdlPredict{:}});
% svmPredict = SVMModel.predict(compactData);
% 
% mdlEval = mdlPredict==idx;
% svmEval = svmPredict'==idx;
% 
% cMdlEval = cell([1,size(idx)]);
% cMdlEval(mdlEval == 1) = {'g'};
% cMdlEval(mdlEval == 0) = {'r'};
% 
% cSvmEval = cell([1,size(idx)]);
% cSvmEval(svmEval == 1) = {'g'};
% cSvmEval(svmEval == 0) = {'r'};
% 
% for k = 1:cellN
%     plot([0,1],[k,k],cMdlEval{k},'linewidth',myLineW)
%     hold on
%     plot([1.5,2.5],[k,k],cSvmEval{k},'linewidth',myLineW)
% end

%Predict:

load 'D:\processing\david\inputSignals.mat'
load 'D:\processing\david\testpeaksArray.mat'
load 'D:\processing\david\inputImages.mat'
load 'D:\processing\david\inputMovie.mat'

cellImages = permute(inputImages,[2,3,1]);

newData = getCellInfo(inputSignals,testpeaksArray,cellImages,inputMovie);
compactNewData = zeros([size(newData.(dataFields{1}),1),length(dataFields)]);

for k = 1:length(dataFields)
    compactNewData(:,k) = newData.(dataFields{k});
end
mdlPredict = (MdlFinal.predict(compactNewData))';
%mdlPredict = Mdl.predict(newData);
%mdlPredict = str2double({mdlPredict{:}});
validCellMax = mdlPredict;
save('D:\processing\david\Mouse2029\20150316-icx\2015_03_16_p000_mouse2029_NULL000_emAnalysisSorted.mat','validCellMax')


%%
% 
% selectedDF = [2,4,6];
% for k = 1:cellN
%     x = myData.(dataFields{selectedDF(1)});
%     y = myData.(dataFields{selectedDF(2)});
%     if length(selectedDF)==2;plot(x(k),y(k),['.',colVect{k}]);end
%     if length(selectedDF)>2; z = myData.(dataFields{selectedDF(3)});plot3(x(k),y(k),z(k),['.',colVect{k}]);end
%     hold on
%     xlabel(dataFields{selectedDF(1)})
%     ylabel(dataFields{selectedDF(2)})
%     if length(selectedDF)>2; zlabel(dataFields{selectedDF(3)});end
% end

%%
kmeansIdx = kmeans(compactData,2);
kmeansIdx(kmeansIdx==2)=0;
figure;
subplot(1,2,1)
idxIs1 = (kmeansIdx == 1);
plot3(compactData(idxIs1,idx(1)),compactData(idxIs1,idx(2)),compactData(idxIs1,idx(3)),'.','color',[0,0.4470,0.7410])
hold on
plot3(compactData(~idxIs1,idx(1)),compactData(~idxIs1,idx(2)),compactData(~idxIs1,idx(3)),'.','color',[0.8500,0.3250,0.0980])
xlabel(dataFields{idx(1)})
ylabel(dataFields{idx(2)})
zlabel(dataFields{idx(3)})
legend('Class1','Class2')

subplot(1,2,2)
idxIsSame = (kmeansIdx == allValidCellMax');
plot3(compactData(idxIsSame,idx(1)),compactData(idxIsSame,idx(2)),compactData(idxIsSame,idx(3)),'g.')
hold on
plot3(compactData(~idxIsSame,idx(1)),compactData(~idxIsSame,idx(2)),compactData(~idxIsSame,idx(3)),'r.')
xlabel(dataFields{idx(1)})
ylabel(dataFields{idx(2)})
zlabel(dataFields{idx(3)})
legend('Same class','Different class')

xlabel(dataFields{idx(1)})
ylabel(dataFields{idx(2)})
zlabel(dataFields{idx(3)})