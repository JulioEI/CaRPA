moviePath = 'C:\Users\csc\Desktop\caImagg\Behavior\Linear-track\Mouse-2028\Behavioral-videos\Mouse2028_20150305_132108.avi';
% mousePath = 'E:\Processing\newTurboreg\Mouse2028_500\20150303-icx';
mousePath = 'E:\Processing\Mouse2028 - 500itEM\20150305-icx';
% [position,~] = getMouseTrajectory(moviePath);
% 
% %%%upscale em%%
% hinf = hdf5info([mousePath,'\2015_03_05_p000_mouse2028_NULL000_turboreg_crop_dfof_1.h5']);
% inputMovie20hz = hdf5read(hinf.GroupHierarchy.Datasets);
% robustMovie = inputMovie20hz;
% robustMovie(isnan(inputMovie20hz))=min(inputMovie20hz(:));
% [upScaledPhi, upFilteredTraces] = recalcPhiAndDetectEvents(robustMovie,emAnalysisOutput.cellImages,emAnalysisOutput.dsCellTraces,emAnalysisOutput.CELLMaxoptions,'runEventDetection',0);
% upScaledPhi = upScaledPhi*10^5;
% save([mousePath,'\upFilteredTraces'],'upFilteredTraces')
% save([mousePath,'\upScaledPhi'],'upScaledPhi')
% % %%%%%%%%%%%%%%%
% %%%predict em%%
% cInf = cellInfo(emAnalysisOutput.CELLMaxoptions.movieFilename,permute(emAnalysisOutput.cellImages,[3,1,2]),emAnalysisOutput.scaledProbability,[]);
% cInf.predictAll('showProgress',1,'tresholds',{'getOverlap','>=0.6','getCellShapePropieties.Eccentricity','<=0.9','getCellShapePropieties.Solidity','>0.8'},'eventPercent',.2);
% validCellMax = logical(cInf.validCells);
% save([mousePath,'\2015_03_07_p000_mouse2028_NULL000_emAnalysisSorted'],'validCellMax')
%%%%%%%%%%%%%%%

%Upscale ICA
% [upScaledPhi] = extractFTfromMovie('E:\Processing\Mouse2028 - 500itEM\20150303-icx\2015_03_03_p000_mouse2028_NULL000_turboreg_crop_dfof_1.h5',pcaicaAnalysisOutput.IcaFilters);
% % Predict ICA
% cInf = cellInfo('E:\Processing\Mouse2028 - 500itEM\20150303-icx\2015_03_03_p000_mouse2028_NULL000_turboreg_crop_dfof_downsampleTime_1.h5',permute(pcaicaAnalysisOutput.IcaFilters,[3,1,2]),pcaicaAnalysisOutput.IcaTraces,[]);
% cInf.predictAll('showProgress',1,'tresholds',{'getOverlap','>=0.6','getCellShapePropieties.Eccentricity','<=0.9','getCellShapePropieties.Solidity','>0.8'},'eventPercent',.2);
% validCells = logical(cInf.validCells);
%%%%%%%

scaledPhi = upScaledPhi(logical(validCells),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mov = VideoReader(moviePath);

period = 1/20;

normalizedSP = scaledPhi./repmat(max(scaledPhi')',[1,size(scaledPhi,2)]);
x = normalizedSP';
% 
signalPeakIdx = [];
validCellsIdx = find(validCells);
normalizedUT = upFilteredTraces./repmat(max(upFilteredTraces')',[1,size(upFilteredTraces,2)]);
spikeDeconv = zeros(size(scaledPhi'));
for k = 1:length(validCellsIdx)
    [c, s, options] = deconvolveCa(normalizedUT(validCellsIdx(k),:)','exp2','foopsi');
    spikeDeconv(:,k) = c;
    signalPeakIdx{k} = find(s);
end

myXAxis = 0.05*(0:1:size(scaledPhi,2)-1);
figure;plot(myXAxis,normalizedUT(validCellsIdx(k),:));hold on;
plot(myXAxis,spikeDeconv(:,k));plot(myXAxis,scaledPhi(k,:)');
line([0.05*signalPeakIdx{k} 0.05*signalPeakIdx{k}],get(gca,'YLim'),'Color',[.8 .8 .8]);legend('Raw trace','EM probability','Denoised trace','UnderlyingSpikeTrain')
xlabel('dfof');xlabel('sec')
% 
% x = spikeDeconv;
% 
% robustInputSignals = double(scaledPhi);
% robustInputSignals(isnan(robustInputSignals)) = 0;
% [~, signalPeakIdx] = computeSignalPeaks(robustInputSignals,'numStdsForThresh',2);
% spikeBinned = zeros(size(scaledPhi));
% winPad = -30:30;
% for k = 1:size(spikeBinned,1)
%     winIdx = bsxfun(@plus,signalPeakIdx{k},winPad);
%     winIdx = reshape(winIdx,[1,numel(winIdx)]);
%     winIdx(winIdx<1)=[];
%     winIdx(winIdx>size(scaledPhi,2))=[];
%     spikeBinned(k,winIdx) = 1;
% end
% % figure;spy(spikeBinned)
% % axis square
% x = spikeBinned';

% integrationWin = 50;
% meanFireVect = zeros(size(spikeBinned));
% for k = 1:size(spikeBinned,2)
%     window = (k-integrationWin/2):(k+integrationWin/2);
%     window(window<=0) = [];
%     window(window>size(spikeBinned,2)) = [];
%     meanFireVect(:,k) = mean(spikeBinned(:,window),2);
% end
% 
% x = meanFireVect';


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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%ELIMINATE SLOW BINS%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

velY = abs(diff(position(:,1)));
idxFast = find(velY > 0);
x = x(idxFast,:);
y = y(idxFast);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%PLOT THE TRAJECTORIES%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Train with 90%, show results on the other 10% movie
trainTestRatio = 0.9;
trainX = x(1:ceil(length(y)*trainTestRatio),:);
testX = x(ceil(length(y)*trainTestRatio)+1:end,:);
trainY = y(1:ceil(length(y)*trainTestRatio));
testY = y(ceil(length(y)*trainTestRatio)+1:end);

%mdl = TreeBagger(60,trainX,trainY,'Method','classification');
mdl = fitcnb(trainX,trainY);
predYTest = mdl.predict(testX)';
if iscell(predYTest)
    predYTest = cellfun(@str2num,predYTest);
end

%Test shuffle (temp)
timeFramesInBin = cell([1,binN]);
for k = 1:binN
    timeFramesInBin{k} = find(testY==k);
end
testXS = testX;
for j = 1:size(testX,2)
    for k = 1:binN
        originalX = testX(timeFramesInBin{k},j);
        testXS(timeFramesInBin{k},j) = originalX(randperm(length(originalX)));
    end    
end
predYTestShuffled = mdl.predict(testXS)';

%Same but shuffling
timeFramesInBin = cell([1,binN]);
for k = 1:binN
    timeFramesInBin{k} = find(y==k);
end

shuffledX = x;
for j = 1:size(x,2)
    for k = 1:binN
        originalX = x(timeFramesInBin{k},j);
        shuffledX(timeFramesInBin{k},j) = originalX(randperm(length(originalX)));
    end
end
trainXShuffled = shuffledX(1:ceil(length(y)*trainTestRatio),:);
testXShuffled = shuffledX(ceil(length(y)*trainTestRatio)+1:end,:);
mdlShuffled = fitcnb(trainXShuffled,trainY);
predYTestShuffled = mdl.predict(testXShuffled)';
if iscell(predYTest)
    predYTestShuffled = cellfun(@str2num,predYTestShuffled);
end

% figure;histogram(predError);title('Error histogram');
figure;plot(testY);hold on;plot(predYTest);plot(predYTestShuffled);legend('GT','Decoded','Decoded shuffled')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%VISUALIZE THE PREDICTIONS%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Save
% 
videoName = 'vidONLYPred';
frames = 100;
frameRate = 20;
frames = min([frames,mov.NumberOfFrames]);
v = VideoWriter(videoName);
v.FrameRate = frameRate;
v.Quality = 100;
%Smooth trajectories
b = 0.2*[1,1,1,1,1];%filter(b,1,testY)
% gtTR = position(ceil(length(y)*trainTestRatio)+1:end,1);
gtTR = ((testY-0.7)/binN)*range(position(:,1))+min(position(:,1));
% (testY-1)*range(position(:,1))/binN + min(position(:,1));
gtS = filtfilt(b,1,gtTR);

predTR = ((predYTest-0.7)/binN)*range(position(:,1))+min(position(:,1));
%(predYTest-1)*range(position(:,1))/binN + min(position(:,1));
predS = filtfilt(b,1,predTR);
predShufTR = ((predYTestShuffled-0.7)/binN)*range(position(:,1))+min(position(:,1)); 
%(predYTestShuffled-1)*range(position(:,1))/binN + min(position(:,1));
predSShuffled = filtfilt(b,1,predShufTR);

open(v)
h = figure('units','normalized','outerposition',[0 0 1 1]);
myCMap = lines;
t = 0;
for k = idxFast(ceil(length(y)*trainTestRatio)+1:end)'%idxFast(idxFast > ceil(length(y)*trainTestRatio)+1)'
    t = t + 1;
    imshow(mat2gray(mov.read(k)));
    caxis([0,1]);
    hold on;
%     plot(gtS(t),mov.height/2-5,'.','markersize',40)
%     plot(predS(t),mov.height/2,'.','markersize',40)
    plot(predSShuffled(t),mov.height/2+5,'.','markersize',40)
    hold off
%     legend('Ground truth','Prediction','Trial/response shuffle','Location','northoutside','Orientation','horizontal')%,'Prediction','Prediction Shuffled')
%     pause(1/30)
    v.writeVideo(getframe(h));
end
close(v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
figure;
mov = VideoReader(moviePath);
myCMap = lines;
for k = ceil(length(y)*trainTestRatio)+1:mov.numberOfFrames
    imshow(mat2gray(mov.read(k)));
    caxis([0,1]);
    hold on;
    for i = 1:binN
    leftXVert = (i-1)*range(position(:,1))/binN + min(position(:,1));
    rightXVert = i*range(position(:,1))/binN + min(position(:,1));
    xVert = [leftXVert,rightXVert,rightXVert,leftXVert];
    yVert = [0,0,mov.Height/2,mov.Height/2];
    if predYTest(k-ceil(length(y)*trainTestRatio)) == i
        patch('XData',xVert,'YData',yVert,'FaceColor',myCMap(i,:),'EdgeColor','none','FaceAlpha',0.6);
    else
        patch('XData',xVert,'YData',yVert,'FaceColor','none','EdgeColor',myCMap(i,:),'EdgeAlpha',0.2); 
    end
    end
    
    for i = 1:binN
    leftXVert = (i-1)*range(position(:,1))/binN + min(position(:,1));
    rightXVert = i*range(position(:,1))/binN + min(position(:,1));
    xVert = [leftXVert,rightXVert,rightXVert,leftXVert];
    yVert = [mov.Height/2,mov.Height/2,mov.Height,mov.Height];
%     if predYTest2(k-ceil(length(y)*trainTestRatio)) == i%y(k) == i
%         patch('XData',xVert,'YData',yVert,'FaceColor',myCMap(i,:),'EdgeColor','none','FaceAlpha',0.6);
%     else
%         patch('XData',xVert,'YData',yVert,'FaceColor','none','EdgeColor',myCMap(i,:),'EdgeAlpha',0.2); 
%     end
    end
    
    hold off
    drawnow;
end

%%
%Bootstrap with length(y)/n data points, in uniform chunks
n = 10;
bootstrapError = zeros([1,n]);
%Subdivide data in n chunks
xBoot = zeros([floor(length(y)/n),size(x,2),n]);
yBoot = zeros([floor(length(y)/n),n]);
for k = 1:n
    xBoot(:,:,k) = x(1+(k-1)*floor(length(y)/n):k*floor(length(y)/n),:);
    yBoot(:,k) = y(1+(k-1)*floor(length(y)/n):k*floor(length(y)/n));
end
%Do the bootstraping
for k = 1:n
    trainX = reshape(permute(xBoot(:,:,setdiff(1:n,k)),[2,1,3]),size(x,2),[])';
    testX = xBoot(:,:,k);
    trainY = reshape(yBoot(:,setdiff(1:n,k)),[1,numel(yBoot(:,setdiff(1:n,k)))]);
    testY = yBoot(:,k); 
    mdl = fitcnb(trainX,trainY);
    predYTest = mdl.predict(testX);
    bootstrapError(k) = mean(sqrt((predYTest-testY).^2));
end
figure;bar(bootstrapError);title(['Bootstraping in ',num2str(n),' chunks']);ylabel('rsquared error');xlabel('blocks')

%%
%Look at how error changes with cell number
n = 10;
cellDivisions = 100;
selectCells = randperm(size(x,2));
cellDivVect = round(linspace(1,size(x,2),cellDivisions));
errorBootMean = zeros([1,cellDivisions]);
parfor j = 1:cellDivisions
    x2 = x(:,selectCells(1:cellDivVect(j)));
    bootstrapError = zeros([1,n]);
    %Subdivide data in n chunks
    xBoot = zeros([floor(length(y)/n),size(x2,2),n]);
    yBoot = zeros([floor(length(y)/n),n]);
    for k = 1:n
        xBoot(:,:,k) = x2(1+(k-1)*floor(length(y)/n):k*floor(length(y)/n),:);
        yBoot(:,k) = y(1+(k-1)*floor(length(y)/n):k*floor(length(y)/n));
    end
    %Do the bootstraping
    for k = 1:n
        trainX = reshape(permute(xBoot(:,:,setdiff(1:n,k)),[2,1,3]),size(x2,2),[])';
        testX = xBoot(:,:,k);
        trainY = reshape(yBoot(:,setdiff(1:n,k)),[1,numel(yBoot(:,setdiff(1:n,k)))]);
        testY = yBoot(:,k); 
        mdl = fitcnb(trainX,trainY);
        predYTest = mdl.predict(testX);
        bootstrapError(k) = mean(sqrt((predYTest-testY).^2));
    end
    errorBootMean(j) = mean(bootstrapError);
end
figure;plot(cellDivVect(end:-1:1),errorBootMean(end:-1:1))
xlabel('cellNumber')
ylabel('mean rsquared error')
%%
%Bootstrap with 90% data n times, taking random values
n = 100;
trainTestRatio = 0.9;
bootstrapError = perFrameBootstrap(x,y,n,trainTestRatio);
figure;histogram(bootstrapError);title(['Bootstraping ',num2str(n),' times taking random frames']);xlabel('rsquared error');
%%
%Error decay with cell number, taking random frames
n = 20;
trainTestRatio = 0.9;
cellDivisions = 50;
nSamplings = 10;
cellDivVect = round(linspace(1,size(x,2),cellDivisions));
errorBootMean = zeros([nSamplings,cellDivisions]);

for i = 1:nSamplings
    disp(['Computing sampling ',num2str(i),' of ',num2str(nSamplings)])
    selectCells = randperm(size(x,2));
    parfor j = 1:cellDivisions
        errorBootMean(i,j) = mean(perFrameBootstrap(x(:,selectCells(1:cellDivVect(j))),y,n));
    end
end
errorBootMeanSamp = 5*mean(errorBootMeanPerFrame,1);
figure;plot(cellDivVect(end:-1:1),errorBootMeanSamp(end:-1:1))
xlabel('cellNumber')
ylabel('mean rsquared error')
axis square
grid on
title([num2str(cellDivisions),' datapoints ',num2str(n), 'times bootstrapping of random frames, averaged ',num2str(nSamplings),' times']);
figure;plot(cellDivVect(end:-1:1),log(errorBootMeanSamp(end:-1:1)))
xlabel('cellNumber')
ylabel('Log of mean rsquared error')
axis square
grid on
title([num2str(cellDivisions),' datapoints ',num2str(n), 'times bootstrapping of random frames, averaged ',num2str(nSamplings),' times']);

%%
%Error decay with cell number, per division

n = 10;
cellDivisions = 50;
nSamplings = 10;
cellDivVect = round(linspace(1,size(x,2),cellDivisions));
errorBootMean = zeros([nSamplings,cellDivisions]);

for i = 1:nSamplings
    disp(['Computing sampling ',num2str(i),' of ',num2str(nSamplings)])
    selectCells = randperm(size(x,2));
    parfor j = 1:cellDivisions
        errorBootMean(i,j) = mean(spatialBootstrap(shuffledX(:,selectCells(1:cellDivVect(j))),y,n));
    end
end
errorBootMeanSamp = mean(errorBootMean,1);
figure;plot(cellDivVect(end:-1:1),errorBootMeanSamp(end:-1:1))
xlabel('cellNumber')
ylabel('mean rsquared error')
title([num2str(cellDivisions),' datapoints ',num2str(n), 'partitions bootstrapping, averaged ',num2str(nSamplings),' times']);
figure;plot(cellDivVect(end:-1:1),log(errorBootMeanSamp(end:-1:1)))
xlabel('cellNumber')
ylabel('Log of mean rsquared error')
title([num2str(cellDivisions),' datapoints ',num2str(n), 'partitions bootstrapping, averaged ',num2str(nSamplings),' times']);