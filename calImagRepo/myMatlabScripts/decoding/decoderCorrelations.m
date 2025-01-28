moviePath = 'C:\Users\csc\Desktop\caImagg\Behavior\Linear-track\Mouse-2028\Behavioral-videos\Mouse2028_20150303_103504.avi';
mousePath = 'E:\Processing\Mouse2028 - 500itEM\20150303-icx\';
% [position,~] = getMouseTrajectory(moviePath);
% load([mousePath,'upFilteredTraces.mat'])
load([mousePath,'upScaledPhi.mat'])
load([mousePath,'validCells.mat'])
scaledPhi = upScaledPhi(logical(validCells),:);

normalizedSP = scaledPhi./repmat(max(scaledPhi')',[1,size(scaledPhi,2)]);
x = normalizedSP';

normalizedPosX = (position(:,1)-min(position(:,1)))/(max(position(:,1))-min(position(:,1)));
binN = 20;
y = zeros([1,length(normalizedPosX)]);
for k = 1:binN
    if k ~= binN
        idx = find(normalizedPosX>=(k-1)/binN&normalizedPosX<(k/binN));
    else
        idx = find(normalizedPosX>=(k-1)/binN&normalizedPosX<=(k/binN));
    end
    y(idx) = k;
end


%Find timeframes in the spatial bins
timeFramesInBin = cell([1,binN]);
for k = 1:binN
    timeFramesInBin{k} = find(y==k);
end
%Shufle the timeframes for each neuron
% shuffledX = x;
% for j = 1:size(x,2)
%     for k = 1:binN
%         originalX = x(timeFramesInBin{k},j);
%         shuffledX(timeFramesInBin{k},j) = originalX(randperm(length(originalX)));
%     end
% end

%%%%%%%%%%%%%%%
binLength = range(position(:,1))/binN;
velocity = abs(diff(position(:,1)));
hist = histogram(velocity(velocity>1.3));
pks = findpeaks(hist.Values);
acuracyLow = find(pks(1)==hist.Values)/binLength;
acuracyHigh = find(pks(end)==hist.Values)/binLength;
%%%%%%%%%%%%%%%%

n = 10;
cellDivisions = 50;
nSamplings = 10;
nShuffles = 10;
cellDivVect = round(linspace(1,size(x,2),cellDivisions));
errorBootMean = zeros([nSamplings,cellDivisions]);
errorBootMeanShuffled = zeros([nSamplings,cellDivisions,nShuffles]);
for i = 1:nSamplings
    disp(['Computing sampling ',num2str(i),' of ',num2str(nSamplings)])
    selectCells = randperm(size(x,2));
    for j = 1:cellDivisions
        errorBootMean(i,j) = mean(spatialBootstrap(x(:,selectCells(1:cellDivVect(j))),y,n));
        for k = 1:nShuffles
            %Shufle the timeframes for each neuron
            shuffledX = x;
            for q = 1:size(x,2)
                for r = 1:binN
                    originalX = x(timeFramesInBin{r},q);
                    shuffledX(timeFramesInBin{r},q) = originalX(randperm(length(originalX)));
                end
            end
            errorBootMeanShuffled(i,j,k) = mean(spatialBootstrap(shuffledX(:,selectCells(1:cellDivVect(j))),y,n));
        end
    end
end
errorBootMeanSamp = mean(errorBootMean,1);
errorShuffled = squeeze(mean(errorBootMeanShuffled,1));
errorShuffledMean = mean(errorShuffled,2);
errorShuffledStd = std(errorShuffled,[],2);
figure;plot(cellDivVect(end:-1:1),errorBootMeanSamp(end:-1:1));
hold on;
shadedErrorBar(cellDivVect(end:-1:1),errorShuffledMean(end:-1:1),errorShuffledStd(end:-1:1),'r')
xlabel('cellNumber')
ylabel('mean rsquared error')
title([num2str(cellDivisions),' datapoints ',num2str(n), 'partitions bootstrapping, averaged ',num2str(nSamplings),' times']);
legend('Original','Frame shuffled')

plot(get(gca,'xlim'), [acuracyLow acuracyLow],'k--')
plot(get(gca,'xlim'), [acuracyHigh acuracyHigh],'k--')


