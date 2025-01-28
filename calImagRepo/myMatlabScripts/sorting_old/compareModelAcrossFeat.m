%Train:

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:

ratioTrainTest = 0.3;
averageSize = 1000;
featSelection = {1:numel(fields(myData))};%,[6,9,12,14,17,19,21],[1,2,3,4,5,7,8,10,11,13,15,16,18,20,22]};
permuteAcrossDays = false;
models = {'TreeBagger','SVM','Naive Bayes','Discriminant Analysis','Nearest Neighbors'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%
[trainScoreAvg,testScoreAvg, trainScoreStd, testScoreStd, ...
    trainFposAvg,testFposAvg, trainFposStd, testFposStd] = deal(zeros([length(featSelection),length(models)]));

featNames = cell([1,length(featSelection)]);
allValidCellMax(allValidCellMax==3)=0;    
for i = 1:length(featSelection)
    
    disp([num2str(100*(i/length(featSelection))),'%'])

    dataFields = fieldnames(myData);
    dataFields = dataFields(featSelection{i});

    featNames{i} = dataFields;
        
    compactData = zeros([size(myData.(dataFields{1}),1),length(dataFields)]);

    for k = 1:length(dataFields)
%         tempData = myData.(dataFields{k});
%         tempData(isnan(tempData)) = max(tempData);
%         myData.(dataFields{k}) = tempData;
        compactData(:,k) = myData.(dataFields{k}); %No normalitzation
        %compactData(:,k) = (myData.(dataFields{k}) - min(myData.(dataFields{k})))./(max(myData.(dataFields{k})) - min(myData.(dataFields{k})));
    end

    [trainScore,testScore, trainFpos, testFpos] = deal(zeros([averageSize,length(models)]));
    for k = 1:averageSize   

        if permuteAcrossDays
            batchIndx = [1+cumsum(cellsPerBatch)'-cellsPerBatch',cumsum(cellsPerBatch)'];
            randbatchIndx = batchIndx(randperm(size(batchIndx,1)),:);
            elemTrain = ceil(ratioTrainTest*size(randbatchIndx,1));
            trainIndx = randbatchIndx(1:elemTrain,:);
            testIndx = randbatchIndx(elemTrain+1:end,:); 
            trainX = [];
            trainY = [];
            testX = [];
            testY = [];
            for l = trainIndx'
                trainX = [trainX;compactData(l(1):l(2),:)];
                trainY = [trainY,allValidCellMax(l(1):l(2))];
            end
            for l = testIndx'
                testX = [testX;compactData(l(1):l(2),:)];
                testY = [testY,allValidCellMax(l(1):l(2))];
            end
        else
            randIndx = randperm(size(compactData,1));
            elemTrain = ceil(ratioTrainTest*size(randIndx,2));
            %trainX = compactData(randIndx(1:elemTrain),:);
            trainY = allValidCellMax(randIndx(1:elemTrain));
            testX = compactData(randIndx(elemTrain+1:end),:); 
            testY = allValidCellMax(randIndx(elemTrain+1:end));    
        end
 
        %------------------------%
        %Trimming data with basic rules
        allDataFields = fields(myData);
        for l = 1:length(allDataFields)
            trainXData.(allDataFields{l}) = myData.(allDataFields{l})(randIndx(1:elemTrain));
            testXData.(allDataFields{l}) = myData.(allDataFields{l})(randIndx(elemTrain+1:end));               
        end
        [predTrainBasic,trainXData] = rulePrd(trainXData);
        [predTestBasic] = rulePrd(testXData);
        
        trainX = zeros([size(trainXData.(dataFields{1}),1),length(dataFields)]);
        for l = 1:length(dataFields)
            trainX(:,l) = trainXData.(dataFields{l});
        end
        trainY = trainY(predTrainBasic);
        clear trainXData testXData
        %------------------------%
            
        
        for j = 1:length(models)
    
            model = models{j};
            
%             mask = find(sum(isnan(trainX),2)==0);
%             trainX = trainX(mask,:);
%             trainY = trainY(mask);
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
            predYTest = predYTest & predTestBasic';
            trainFpos(k,j) = sum(trainY == 0 & predYTrain == 1)/length(trainY);
            testFpos(k,j) = sum(testY == 0 & predYTest == 1)/length(testY);     
            trainScore(k,j) = sum(trainY~=predYTrain)/length(predYTrain);
            testScore(k,j) = sum(testY~=predYTest)/length(predYTest);
            clear mdl
        end
    end
    clear compactData

    trainScoreAvg(i,:) = mean(trainScore);
    testScoreAvg(i,:) = mean(testScore);
    trainFposAvg(i,:) = mean(trainFpos);
    testFposAvg(i,:) = mean(testFpos);
    trainScoreStd(i,:) = std(trainScore);
    testScoreStd(i,:) = std(testScore);
    trainFposStd(i,:) = std(trainFpos);
    testFposStd(i,:) = std(testFpos);    

end

dispNames = [];
for i = 1:length(featNames)
    tempDispNames = featNames{i};
    tempDispNamesConcat = [];
    for k = 1:length(tempDispNames)
        tempDispNamesConcat = [tempDispNamesConcat,tempDispNames{k},' '];
    end
    dispNames = [dispNames,{tempDispNamesConcat}];
end
%%
figure;
if length(featSelection) > 1
    h = barwitherr(testScoreStd,testScoreAvg);
    hold on
    bar(testFposAvg);

    %Manually print errors for the false positive    
    nCols = length(models);
    hBar = h;
    values = testFposAvg;
    xOrder = 1:size(values,1);
    hErrorbar = zeros(1,nCols);
    for col = 1:nCols
        % Extract the x location data needed for the errorbar plots:
        x = bsxfun(@plus, hBar(col).XData, [hBar(col).XOffset]');
        % Use the mean x values to call the standard errorbar fn; the
        % errorbars will now be centred on each bar; these are in ascending
        % order so use xOrder to ensure y values and errors are too:
        hErrorbar(col) = errorbar(mean(x,1), values(xOrder,col), testFposStd(xOrder,col), testFposStd(xOrder, col), '.k');
        set(hErrorbar(col), 'marker', 'none')
    end
    
    legend(h,models);
    set(gca,'xticklabel',dispNames)
    fix_xticklabels();

else
    h = barwitherr(testScoreStd,testScoreAvg);
    hold on
    bar(testFposAvg);

    %Manually print errors for the false positive    
    nCols = length(models);
    hBar = h;
    values = testFposAvg;
    xOrder = 1:size(values,1);
    hErrorbar = zeros(1,nCols);
    for col = 1:nCols
        % Extract the x location data needed for the errorbar plots:
        x = bsxfun(@plus, hBar.XData(col), [hBar.XOffset(1)]');
        % Use the mean x values to call the standard errorbar fn; the
        % errorbars will now be centred on each bar; these are in ascending
        % order so use xOrder to ensure y values and errors are too:
        hErrorbar(col) = errorbar(mean(x,2), values(xOrder,col), testFposStd(xOrder,col), testFposStd(xOrder, col), '.k');
        set(hErrorbar(col), 'marker', 'none')
    end
    set(gca,'xticklabel',models)   
end
title({['Error on unknown test partition averaging over ', num2str(averageSize),' random partions with '],[ num2str(100*ratioTrainTest),'% data for training and ',num2str(100*(1-ratioTrainTest)),'% for testing']})
ylabel('Error')


% figure;
% if length(featSelection) > 1
%     h = barwitherr(testFposStd,testFposAvg);
%     legend(h,models);
%     set(gca,'xticklabel',dispNames)
%     fix_xticklabels();
% else
%     barwitherr(testFposStd,testFposAvg);
%     set(gca,'xticklabel',models)   
% end
% title({['Error on unknown test partition averaging over ', num2str(averageSize),' random partions with '],[ num2str(100*ratioTrainTest),'% data for training and ',num2str(100*(1-ratioTrainTest)),'% for testing']})
% ylabel('False positives')


% for k = 1:length(models)
%     [y,myStd] = deal(zeros([size(testAvg,2),2]));
%     for i = 1:size(testAvg,1)
%         y(i,:) = [testAvg(i,1),trainAvg(i,1)];
%         myStd(i,:) = [testStd(i,k),trainStd(i,k)];
%     end
% 
%     figure;
%     title([models{k},' averaging over ', num2str(averageSize),' random partions')
%     h = barwitherr(myStd,y);
%     legend(h,{['Test error (',num2str(100*(1-ratioTrainTest)),%')'],['Train error (',num2str(100*ratioTrainTest),%')']});
%     set(gca,'xticklabel',dispNames)
%     fix_xticklabels();
% end


