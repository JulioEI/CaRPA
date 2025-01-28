function [predXTest, mixturePoissLogLikeAllCells] = poissonMixtureDecoder(trainR,trainX,testR)

[~,~,~,mixturePoissLogLikeAllCells] = PoissonMixtureSpkCountLikelihood(trainX,trainR,max(trainX));

predXTest = zeros([size(testR,1),1]);
for t = 1:size(testR,1)
    [~,maxIdx] = max(mixturePoissLogLikeAllCells(testR(t,:)'));
    predXTest(t) = maxIdx;
end

