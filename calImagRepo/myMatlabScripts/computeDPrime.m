function [dVect,dPrimeDist,other] = computeDPrime(x,r,param,shuffleTrials)
    [x,r] = dataAnalysis.curateXandR(x,r,param); %Bin, remove slow frames
    x = x(:,1);
    if shuffleTrials
        r = decoderAnalysis.shuffleKeepingTrials1D(r,x);
    end
    binN = param.binN(1);
    %%
    disp('Separating bins into matcices')
    myMat = cell([1,binN]);
    for bin = 1:max(x)
        indWithBin = find(x==bin);
        myMat{bin} = r(indWithBin,:);
    end
    %%
    disp('Sorting data into train and test for each kfold and pair of classes')
    kFolds = 10;
    numPairs = (binN*(binN+1))/2;
    pairVect = zeros([2,numPairs]);
    [testData.X,testData.Y,trainData.X,trainData.Y] = deal(cell([kFolds,numPairs]));
    idx = 0;
    for i = 1:binN
        for j = i:binN 
            idx = 1 + idx;
            pairVect(:,idx) = [i,j];

            y = [myMat{i};myMat{j}];
            x = [ones(size(myMat{i},1),1);-1*ones(size(myMat{j},1),1)];

            %Permute so trials are mixed
            randIdx = randperm(length(x));
            y = y(randIdx,:);
            x = x(randIdx);

            %Compute kfolds
            N = length(x);
            foldSize = floor(N/kFolds);
            Nfold = foldSize*kFolds;
            for fold = 1:kFolds
                testIdx = (1+(fold-1)*foldSize):(fold*foldSize);
                testData.X{fold,idx} = x(testIdx);
                testData.Y{fold,idx} = y(testIdx,:);

                trainIdx = setdiff(1:Nfold,testIdx);
                trainData.X{fold,idx} = x(trainIdx);
                trainData.Y{fold,idx} = y(trainIdx,:);
            end
        end
    end
    %%
    disp('Training a svm and fit prob for each pair and kfold')
    [models,transforms] = deal(cell([kFolds,numPairs]));
    parfor idx = 1:numPairs
        for fold = 1:kFolds
            CVSVMModel= fitcsvm(trainData.Y{fold,idx},trainData.X{fold,idx});
            [models{fold,idx},transforms{fold,idx}] = fitSVMPosterior(CVSVMModel);
        end
    end
    mySvm.model = models;
    mySvm.transform = transforms;
    %%
    disp('Testing for each pair and finding d''')
    dPrime = zeros([kFolds,numPairs]);
    parfor idx = 1:numPairs
        for fold = 1:kFolds
            [pred,predProb] = mySvm.model{fold,idx}.predict(testData.Y{fold,idx});
            %probs_correct = (testData.X{fold,idx} == 1).*predProb(:,2) + (testData.X{fold,idx} == -1).*predProb(:,1);
            %dPrime(fold,idx) = mean(probs_correct);
            plusIdx = find(pred == 1);
            negIdx = find(pred == -1);
%             
%             %%%%%WITH DISTANCES%%%%%%%
%             distance = mySvm.model{fold,idx}.margin(testData.Y{fold,idx},testData.X{fold,idx});
% %             %Not working
% %             testY = testData.Y{fold,idx};
% %             bias = mySvm.model{fold,idx}.Bias;
% %             sv = mySvm.model{fold,idx}.SupportVectors;
% %             alphaHat = mySvm.model{fold,idx}.Alpha;
% %             wT = (alphaHat'*sv) + bias;
% %             distance = (wT*testY')./norm(wT);
% %             plusPredProb = distance(plusIdx);
% %             negPredProb =  distance(negIdx);
% %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
%             %%%%%WITH PROBABILITIES%%%%%
            plusPredProb = predProb(plusIdx,2)-0.5;%plusPredProb = 2*(predProb(plusIdx,2)-0.5); %Scale from [0.5,1] to [0,1]
            negPredProb = 1-predProb(negIdx,2)-0.5;%negPredProb = 2*((1-predProb(negIdx,2))-0.5); %Scale from [0,0.5] to [0,1]
            
            
            errorVec = sign((testData.X{fold,idx} - pred));
            plusErrorVec = errorVec(plusIdx);
            negErrorVec = errorVec(negIdx);
            plusPredProb(find(plusErrorVec)) = -plusPredProb(find(plusErrorVec)); %Flip incorrect
            negPredProb(find(negErrorVec)) = -negPredProb(find(negErrorVec));%Flip incorrect
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
            %Compute d'
            if isempty(plusPredProb) && isempty(negPredProb)
                error('There is no test data for this fold and pair')
            elseif isempty(plusPredProb)
                disp(['All test data in pair i:', num2str(pairVect(1,idx)), ' j:', num2str(pairVect(2,idx)),' and fold ', num2str(fold), ' was predicted to lie in a single direction'])                       
                dPrime(fold,idx) = mean(negPredProb);
            elseif isempty(negPredProb)
                disp(['All test data in pair i:', num2str(pairVect(1,idx)), ' j:', num2str(pairVect(2,idx)),' and fold ', num2str(fold), ' was predicted to lie in single direction'])                    
                dPrime(fold,idx) = mean(plusPredProb);
            else
                dPrime(fold,idx) = mean(plusPredProb) + mean(negPredProb); %d' goes from [-2,2]
            end
        end
    end
    %%
    % %Create pairwise matrix (USELESS)
    dPrimeIdx = zeros(binN);
    upperTriangleIndices = triu(true(binN));
    dPrimeIdx(upperTriangleIndices) = 1:numPairs;
    dPrimeIdx = dPrimeIdx + dPrimeIdx.' - diag(diag(dPrimeIdx));
    
    dPrimeMat = zeros(binN);
    for i = 1:binN
        for j = 1:binN
            dPrimeMat(i,j) =mean(dPrime(:,dPrimeIdx(i,j)));
        end
    end
    %Show results
    figure;imagesc(dPrimeMat);colorbar;xlabel('Class i');ylabel('Class j')
    %%
    disp('Averaging results from each distance pair')
    dPrimePerDist = cell([1,binN]);
    for idx = 1:numPairs
        dPrimePerDist{diff(pairVect(:,idx))+1} = [dPrimePerDist{diff(pairVect(:,idx))+1},mean(dPrime(:,idx))];
    end
    dPrimeDist.mean = cellfun(@mean,dPrimePerDist);
    dPrimeDist.ste = cellfun(@(x) std(x)./sqrt(numel(x)),dPrimePerDist);
    dVect = 0:(binN-1);
    %%
    %Save other computed data if needed
    if nargout > 2
        other.dPrime = dPrime;
        other.mySvm = mySvm;
        other.testData = testData;
        other.trainData = trainData;
        other.myMat = myMat;
        other.param = param;
        other.numPairs = numPairs;
        other.pairVect = pairVect;
        other.kFolds = kFolds;
    end
end

