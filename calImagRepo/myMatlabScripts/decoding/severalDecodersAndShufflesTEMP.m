% binN = max(x(:));
% timeFramesInBin = cell([1,binN]);
% for k = 1:binN
%     timeFramesInBin{k} = find(x==k);
% end
% shuffledR = r;
% for i = 1:size(r,2)
%     for k = 1:binN
%         originalX = r(timeFramesInBin{k},i);
%         shuffledR(timeFramesInBin{k},i) = originalX(randperm(length(originalX)));
%     end
% end
%%
decoders = {'Naive Bayes','TreeBagger','Nearest Neighbors'};
    
allShuffledMat = zeros([10,3]);
allShuffledMat(:,1) = spatialBootstrap(shuffledR,x,10,decoders{1});disp(['done ',decoders{1}])
allShuffledMat(:,2) = spatialBootstrap(shuffledR,x,10,decoders{2});disp(['done ',decoders{2}])
allShuffledMat(:,3) = spatialBootstrap(shuffledR,x,10,decoders{3});disp(['done ',decoders{3}])

originalMat = zeros([10,3]);
originalMat(:,1) = spatialBootstrap(r,x,10,decoders{1});disp(['done ',decoders{1}])
originalMat(:,2) = spatialBootstrap(r,x,10,decoders{2});disp(['done ',decoders{2}])
originalMat(:,3) = spatialBootstrap(r,x,10,decoders{3});disp(['done ',decoders{3}])

%%
testShuffledMat = zeros([10,3]);
testShuffledMat(:,1) = spatialBootstrap(r,x,10,decoders{1});disp(['done ',decoders{1}])
testShuffledMat(:,2) = spatialBootstrap(r,x,10,decoders{2});disp(['done ',decoders{2}])
testShuffledMat(:,3) = spatialBootstrap(r,x,10,decoders{3});disp(['done ',decoders{3}])
%%
trainShuffledMat = zeros([10,3]);
trainShuffledMat(:,1) = spatialBootstrap(r,x,10,decoders{1});disp(['done ',decoders{1}])
trainShuffledMat(:,2) = spatialBootstrap(r,x,10,decoders{2});disp(['done ',decoders{2}])
trainShuffledMat(:,3) = spatialBootstrap(r,x,10,decoders{3});disp(['done ',decoders{3}])
%%
figure;
subplot(2,2,1);
boxplot(originalMat,decoders);title('No shuffle');axis square;grid on;ylabel('error (bins)')
subplot(2,2,2);
boxplot(allShuffledMat,decoders);title('Train + test shuffled');axis square;grid on;ylabel('error (bins)')
subplot(2,2,3);
boxplot(testShuffledMat,decoders);title('Only test shuffled');axis square;grid on;ylabel('error (bins)')
subplot(2,2,4);
boxplot(trainShuffledMat,decoders);title('Only train shuffled');axis square;grid on;ylabel('error (bins)')
linkaxes;set(gcf,'color','w')
annotation('textbox',[0,0,1,1],'String','10fold decoder 20bin error with dfof traces','FitBoxToText','on');
%%
figure;
myText = {'No shuffle','Train + test','Only test','Only train'};
for k = 1:3
    subplot(3,1,k);
    boxplot([originalMat(:,k),allShuffledMat(:,k),testShuffledMat(:,k),trainShuffledMat(:,k)],myText);title(decoders{k});grid on;ylabel('error (bins)')
end
linkaxes;set(gcf,'color','w')
annotation('textbox',[0,0,1,1],'String','10fold decoder 20bin error with dfof traces','FitBoxToText','on');