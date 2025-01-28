function [predictedYTrain,predictedYTest] = getPredict(mdl,trainX,testX)
    
    predictedYTrain = mdl.predict(trainX)';
    
    if iscell(predictedYTrain)
        predictedYTrain = cellfun(@str2num,predictedYTrain);
    end
    
    predictedYTest = zeros([1,size(testX,1)]);
    mask = find(sum(isnan(testX),2)==0);
    predictedYTestCrop = mdl.predict(testX(mask,:))';

    if iscell(predictedYTestCrop)
        predictedYTestCrop = cellfun(@str2num,predictedYTestCrop);
    end
    
    predictedYTest(mask) = predictedYTestCrop;
end

