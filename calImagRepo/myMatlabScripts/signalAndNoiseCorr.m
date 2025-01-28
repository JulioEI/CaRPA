load('E:\Processing\Mouse2028 - 500itEM\20150303-icx\upScaledPhi.mat')
load('E:\Processing\Mouse2028 - 500itEM\20150303-icx\validCells.mat')
moviePath = 'C:\Users\csc\Desktop\caImagg\Behavior\Linear-track\Mouse-2028\Behavioral-videos\Mouse2028_20150303_103504.avi';

%%%%%%%%%%%%%%%%
% [position,~] = getMouseTrajectory(moviePath);

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
%%%%%%%%%%%%%%

scaledPhi = upScaledPhi(logical(validCells),:);

x = zscore(scaledPhi)';


%Remove slow bins
velY = abs(diff(position(:,1)));
idxFast = find(velY > 4);
x = x(idxFast,:);
y = y(idxFast);
%

timeFramesInBin = cell([1,binN]);
for k = 1:binN
    timeFramesInBin{k} = find(y==k);
end

%Compute signal correlations
meanSignalPerBin = zeros(binN,size(x,2));
for k = 1:binN
    meanSignalPerBin(k,:) = mean(x(timeFramesInBin{k},:));
end
signalCorr = corr(meanSignalPerBin);
signalCorrV = signalCorr(~tril(ones(size(signalCorr))));

%Compute noise correlations
noiseCorr = zeros([size(x,2),size(x,2),binN]);
for k = 1:binN
    noiseCorr(:,:,k) = corr(x(timeFramesInBin{k},:));
end
noiseCorrV = zeros(binN,length(signalCorrV));
for k = 1:binN
    noiseCorrTemp = noiseCorr(:,:,k);
    noiseCorrV(k,:) = noiseCorrTemp(~tril(ones(size(signalCorr))));
end
% figure;plot(mean(noiseCorrV,2));hold on;plot(mean(noiseCorrV,2)+std(noiseCorrV,[],2),'--r');plot(mean(noiseCorrV,2)-std(noiseCorrV,[],2),'--r')
figure;plot(noiseCorrV,'b.')