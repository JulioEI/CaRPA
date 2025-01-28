function [bootstrapError,errorsPerBin] = spatialBootstrap(x,y,n,varargin)
    bootstrapError = zeros([1,n]);
    errorsPerBin = zeros([n,max(y)]);
    %Subdivide data in n chunks
    xBoot = zeros([floor(length(y)/n),size(x,2),n]);
    yBoot = zeros([floor(length(y)/n),n]);
    for j = 1:n
        xBoot(:,:,j) = x(1+(j-1)*floor(length(y)/n):j*floor(length(y)/n),:);
        yBoot(:,j) = y(1+(j-1)*floor(length(y)/n):j*floor(length(y)/n));
    end
    %Do the bootstraping
    for j = 1:n
        trainX = reshape(permute(xBoot(:,:,setdiff(1:n,j)),[2,1,3]),size(x,2),[])';
        testX = xBoot(:,:,j);
        trainY = reshape(yBoot(:,setdiff(1:n,j)),[1,numel(yBoot(:,setdiff(1:n,j)))]);
        testY = yBoot(:,j);
        
%         binN = max(trainY(:));
%         timeFramesInBin = cell([1,binN]);
%         for k = 1:binN
%             timeFramesInBin{k} = find(trainY==k);
%         end
%         shuffledX = trainX;
%         for i = 1:size(testX,2)
%             for k = 1:binN
%                 originalX = trainX(timeFramesInBin{k},i);
%                 shuffledX(timeFramesInBin{k},i) = originalX(randperm(length(originalX)));
%             end
%         end
%         trainX = shuffledX;

        if isempty(varargin);model = 'Naive Bayes';else;model = varargin{1};end
        switch model
            case 'TreeBagger'
              mdl = TreeBagger(60,trainX,trainY,'Method','classification');
            case 'SVM' 
              mdl = fitcsvm(trainX,trainY,'Standardize',true,'KernelFunction','linear');
            case 'Naive Bayes'
              mdl = fitcnb(trainX,trainY);
            case 'Discriminant Analysis'
              mdl = fitcdiscr(trainX,trainY,'discrimType','pseudoLinear');
            case 'Nearest Neighbors'
              mdl = fitcknn(trainX,trainY);
            otherwise
              error('no model specified')
        end
     
%         binN = max(testY(:));
%         timeFramesInBin = cell([1,binN]);
%         for k = 1:binN
%             timeFramesInBin{k} = find(testY==k);
%         end
%         shuffledX = testX;
%         for i = 1:size(testX,2)
%             for k = 1:binN
%                 originalX = testX(timeFramesInBin{k},i);
%                 shuffledX(timeFramesInBin{k},i) = originalX(randperm(length(originalX)));
%             end
%         end
%         testX = shuffledX;
        
        predYTest = mdl.predict(testX);
        if iscell(predYTest)
            predYTest = cellfun(@str2num,predYTest);
        end
        se = sqrt((predYTest-testY).^2);
        for k = 1:max(y)
            errorsPerBin(j,k) = sum(testY~=predYTest);
        end
        bootstrapError(j) = mean(se);
    end
end

