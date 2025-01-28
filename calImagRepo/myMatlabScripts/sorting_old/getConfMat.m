function [confCoefTrain,confCoefTest,confMatTrain,confMatTest,confIndTrain,confIndTest] = getConfMat(trainY,predictedYTrain,testY,predictedYTest)
    [confCoefTrain,confMatTrain,confIndTrain] = confusion(trainY,predictedYTrain);
    disp(['TRAIN ERROR: ',num2str(confCoefTrain)]);disp(' ')
    disp(confMatTrain)

    [confCoefTest,confMatTest,confIndTest] = confusion(testY,predictedYTest);
    disp(['TEST ERROR: ',num2str(confCoefTest)]);disp(' ')
    disp(confMatTest);disp(' ')
end

