function errorPerFold = quickBayesDecoder(x,y,kFolds)
    %Tests with one fold and trains with the others, k times
    
    N = length(y);
    
    foldSize = floor(N/kFolds);
    Nfold = foldSize*kFolds;
    
    errorPerFold = zeros([1,kFolds]);
    
    parfor k = 1:kFolds
        testIdx = (1+(k-1)*foldSize):(k*foldSize);
        testX = x(testIdx,:);
        testY = y(testIdx);
        
        trainIdx = setdiff(1:Nfold,testIdx);
        trainX = x(trainIdx,:);
        trainY = y(trainIdx,:);
        
        mdl = fitcnb(trainX,trainY);
        
        predYTest = mdl.predict(testX);
        errorPerFold(k) = mean(sqrt((predYTest-testY).^2));
    end
end

