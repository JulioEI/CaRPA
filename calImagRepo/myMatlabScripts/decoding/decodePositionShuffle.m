moviePath = 'C:\Users\csc\Desktop\caImagg\Behavior\Linear-track\Mouse-2028\Behavioral-videos\Mouse2028_20150303_103504.avi';
mousePath = 'E:\Processing\Mouse2028 - 500itEM\20150303-icx\';
% [position,~] = getMouseTrajectory(moviePath);
% load([mousePath,'upFilteredTraces.mat'])
load([mousePath,'upScaledPhi.mat'])
load([mousePath,'validCells.mat'])
scaledPhi = upScaledPhi(logical(validCells),:);

%Divide trials
mov = VideoReader(moviePath);
sess = 1;
period = 1/20;
if size(position,1) ~= mov.numberOfFrames
    error('Movie and calcium have movie differ in frames')
end
realWidth = 1;
metersPerPixel = realWidth/mov.Width;
Trajectory{sess} = [period*(0:(size(position,1)-1))',position*metersPerPixel];
Trajectory{sess}(:,2) = Trajectory{sess}(:,2) - (max(Trajectory{sess}(:,2)) + min(Trajectory{sess}(:,2)))/2; 

%Divide trials%%%%%%%%%%%%%%%%%
absTraX = abs(Trajectory{sess}(:,2));
% figure;plot(absTraX)
absTraXDiff = abs(diff(absTraX));
absTraX(absTraXDiff > median(absTraXDiff)) = nan;
% hold on;plot(absTraX)
traXNaN = find(~isnan(absTraX));
traXNaNDiff = diff(traXNaN);
traXFillIdx = find(traXNaNDiff > 1 & traXNaNDiff < 50); 
traXFill = [traXNaN(traXFillIdx),traXNaN(traXFillIdx)+traXNaNDiff(traXFillIdx)];
for k = 1:size(traXFill,1)
    absTraX(traXFill(k,1):traXFill(k,2)) = abs(Trajectory{sess}(traXFill(k,1):traXFill(k,2),2));
end
absTraX(absTraX<max(absTraX)*.5) = nan;
% hold on;plot(absTraX)
idffIdx = find(~isnan(absTraX));
idffIdxDiff = diff(idffIdx);
trialBoundaries = idffIdx(find(idffIdxDiff>1))';

posAtBoundaries = position(trialBoundaries,1);
dxBetwenBoundaries = abs(diff(posAtBoundaries));
mergeIdx = find(dxBetwenBoundaries<(max(dxBetwenBoundaries)/5));
trialBoundaries(mergeIdx) = [];
trialBoundaries = [trialBoundaries,length(position(:,1))];

% figure;
% plot(position);
% hold on
% plot([trialBoundaries;trialBoundaries],get(gca,'ylim'),'color',[.8,.8,.8])
% hold off

nTrials = length(trialBoundaries)-1;
%%%%%%%%%%%%%%%%%

%Get Y (position)%%%%%%%%%%%%%%%%%

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

yTrial = cell([nTrials,1]);
for k = 1:nTrials
    yTrial{k} = y(trialBoundaries(k):trialBoundaries(k+1));
end

%%%%%%%%%%%%%%%%%

%Get X (f traces)%%%%%%%%%%%%%%%%%
normalizedSP = scaledPhi./repmat(max(scaledPhi')',[1,size(scaledPhi,2)]);
x = normalizedSP';
% validCellsIdx = find(validCells);
% normalizedUT = upFilteredTraces./repmat(max(upFilteredTraces')',[1,size(upFilteredTraces,2)]);
% spikeDeconv = zeros(size(scaledPhi'));
% for k = 1:length(validCellsIdx)
%     [c, s, options] = deconvolveCa(normalizedUT(validCellsIdx(k),:)','exp2','foopsi');
%     spikeDeconv(:,k) = c;
% end
% x = spikeDeconv;

xTrial = cell([nTrials,1]);
for k = 1:nTrials
    xTrial{k} = x(trialBoundaries(k):trialBoundaries(k+1),:);
end
%%%%%%%%%%%%%%%%%
%%
%Assess mean decoder acuracy without shuffling
newX = [];
newY = [];
for k = 1:nTrials
    newX = [newX;xTrial{k}];
    newY = [newY,yTrial{k}];
end
noShuffleScore = mean(spatialBootstrap(newX,newY,10));
%%
%Assess mean decoder acuracy shuffling
nShuffles = 10;
shuffleScore = zeros([nShuffles,1]);
for j = 1:nShuffles
    newX = [];
    newY = [];
    for k = randperm(nTrials)
        newX = [newX;xTrial{k}];
        newY = [newY,yTrial{k}];
    end
    shuffleScore(j) = mean(spatialBootstrap(newX,newY,10));    
end

%%
%Find timeframes in the spatial bins
binN = 20;
timeFramesInBin = cell([1,binN]);
for k = 1:binN
    timeFramesInBin{k} = find(y==k);
end
%Shufle the timeframes for each neuron
shuffledX = x;
for j = 1:size(x,2)
    for k = 1:binN
        originalX = x(timeFramesInBin{k},j);
        shuffledX(timeFramesInBin{k},j) = originalX(randperm(length(originalX)));
    end
end

