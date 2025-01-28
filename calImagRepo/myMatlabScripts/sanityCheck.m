%Load animals
%clear all
home
animals = {
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2022'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2011'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2021'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2024'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2026'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2019'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2012'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2028'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2025'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2023'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2010'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2029'};
sessions = {
    '20150326'
    '20150207'
    '20150326'
    '20150311'
    '20150228'
    '20150228'
    '20150118'
    '20150228'
    '20150228'
    '20150326'
    '20150128'
    '20150311'};
pxPerCm = cell([1,length(animals)]);
load pxPerCm
data = [];
%%
%Create decoder analysis object for each animal
mouseNames = cell([1,length(animals)]);
for k = 1:length(animals)
    [~,folderName,~] = fileparts(animals{k});
    regexpOut = regexp(folderName, 'ouse\D*(\d+)','tokens');
    mouseNames{k} = regexpOut{1}{1};
    data.(['m',mouseNames{k}]).decA = decoderAnalysis(animals{k});
    if isempty(pxPerCm{k})
        [aviFile,aviPath] = uigetfile([animals{k},filesep,'*.avi'],'Select a movie to find px scaling');
        mov = VideoReader([aviPath,filesep,aviFile]);
        frame100 = mov.read(100);
        myFig = figure;imagesc(frame100);title('Delimitate the track')
        rect = getrect(myFig);
        trackLengthPx = rect(3);
        trackWidthPx = rect(4);
        answer = inputdlg({'Length in cm','Width in cm'});
        close(myFig);
        pxPerCm{k} =  [trackLengthPx/str2num(answer{1}),trackWidthPx/str2num(answer{2})];
    end
    data.(['m',mouseNames{k}]).decA.pxPerCm = pxPerCm{k};
end

%%
% Set the parameters
param.dtCamera = 0.05; %period of miniscope (seconds/frame)
param.dT = 0.05; %Integration parameter (seconds)
param.binN = [20,1]; %Number of bins in X Y
param.minVel = [4, 0]; %Frames under this speed will be discarded (cm/s)
param.pxPerCm = data.m2024.decA.pxPerCm; %How many px a centimeter is (px/cm)
param.traceField = 'rawProb'; %Default field to extract traces
%%
sessions = data.m2024.decA.chooseSessions('all');
%Get traces and positions
output = data.m2024.decA.loadOutput(sessions(1));
x = output.position(:,1);
r = output.(param.traceField);
%%
figure;
subplot(2,3,1);
plot(x,'-o');title('Position');grid on;
subplot(2,3,4);
plot(r(:,1:10),'-o');title('10 Responses');grid on;
[xCur,rCur] = dataAnalysis.curateXandR(x,r,param);
subplot(2,3,2);
plot(xCur,'-o');title('Curated position');grid on;
subplot(2,3,5);
plot(rCur(:,1:10),'-o');title('10 Curated responses');grid on;
%
%Randomly permute
randPermIdx = randperm(length(xCur));%1:length(xCur);%
xPerm = xCur(randPermIdx);
rPerm = rCur(randPermIdx,:);
subplot(2,3,3);
plot(xPerm(:),'-o');title('Position permuted');grid on;
subplot(2,3,6);
plot(rPerm(:,1:10),'-o');title('10 Responses permuted');grid on;
%%

% %Show PF taking half a sessoion for permuted and not permuted
%Get activity per bin
PFPerm = cell([max(xPerm),size(rPerm,2)]);
PFCur = cell([max(xPerm),size(rPerm,2)]);
for t = 1:(round(length(xPerm)/2))
    for celli = 1:100%size(rPerm,2)
        PFPerm{xPerm(t),celli} = [PFPerm{xPerm(t),celli},rPerm(t,celli)];
        PFCur{xCur(t),celli} = [PFCur{xCur(t),celli},rCur(t,celli)];
    end
end
meanPFPerm = cellfun(@mean,PFPerm);
stePFPerm = cellfun(@std,PFPerm)./sqrt(cellfun(@length,PFPerm));
meanPFCur = cellfun(@mean,PFCur);
stePFCur = cellfun(@std,PFCur)./sqrt(cellfun(@length,PFCur));
%%
figure;
cellsToShow = 40;
subplotN = numSubplots(cellsToShow);
for celli = 1:cellsToShow%size(rPerm,2)
    subplot_er(subplotN(1),subplotN(2),celli);
    errorbar(meanPFPerm(:,celli),stePFPerm(:,celli),'lineWidth',1.3);
    hold on;
    errorbar(meanPFCur(:,celli),stePFCur(:,celli),'lineWidth',1.3);
    axis square;
    %set(gca,'xtick',[])
    %set(gca,'ytick',[])
    box off
    ylim([0,max(max(meanPFPerm(:,1:cellsToShow)))]);
    grid on;
    if celli == 1
        legend('Permuted','Raw');
    end
end
%annotation('textbox',[0,0,1,1],'string',['PF taking half a session (',param.traceField,')'],'fontSize',20,'fontWeight','bold')

set(gcf,'color','white')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
export_fig asdf.jpg -m3
%save2pdf(['C:\Users\csc\Desktop\caImagg\Graphics\PFhalfsession'])

%%
%Make folds
kFolds = 10;
N = length(xPerm);
foldSize = floor(N/kFolds);
Nfold = foldSize*kFolds;
testIdx = zeros([kFolds,foldSize]);
trainIdx = zeros([kFolds,Nfold-foldSize]);
for k = 1:kFolds
    testIdx(k,:) = (1+(k-1)*foldSize):(k*foldSize);
    trainIdx(k,:) = setdiff(1:Nfold,testIdx(k,:));
end
%%
%Run the svm
errorPerFold = zeros([1,kFolds]);
xVar = xPerm;
rVar = rPerm;
tic;
for k = 1:kFolds
    disp(['Fold ', num2str(k),'/',num2str(kFolds)]);
    testX = xVar(testIdx(k,:));
    testR = rVar(testIdx(k,:),:);
    
    trainX = xVar(trainIdx(k,:));
    trainR = rVar(trainIdx(k,:),:);
    
    %Permute training?
    %randPermIdx = randperm(length(trainX));
    %trainX = trainX(randPermIdx);
    %trainR = trainR(randPermIdx,:);
    
    mdl = fitcecoc(trainR,trainX);
    errorPerFold(k) = mean(sqrt((mdl.predict(testR)-testX').^2));
    
end
toc
disp(errorPerFold)
%%
%Compare closest 2 folds to train with all the other folds (avoid extremes)
separateFolds = cell([1,kFolds]);
for k = 1:kFolds
    separateFolds{k} = (1+(k-1)*foldSize):(k*foldSize);
end

errorPerClosestFold = nan([1,kFolds-2]);
errorPerClosestFoldSte = nan([1,kFolds-2]);
errorPerOtherFold = nan([1,kFolds-2]);
errorPerOtherFoldSte = nan([1,kFolds-2]);
xVar = xPerm;
rVar = rPerm;
tic;
for k = 2:(kFolds-1)
    disp(['Fold ', num2str(k),'/',num2str(kFolds)]);
    testX = xVar(cat(2,separateFolds{k}));
    testR = rVar(cat(2,separateFolds{k}),:);
    
    closestIdx = [k+1,k-1];
    otherIdx = setdiff(1:kFolds,[k,closestIdx]);
    otherIdx = randsample(otherIdx,2);
    
    trainXClosest = xVar(cat(2,separateFolds{closestIdx}));
    trainRClosest = rVar(cat(2,separateFolds{closestIdx}),:);
    
    trainXOther = xVar(cat(2,separateFolds{otherIdx}));
    trainROther = rVar(cat(2,separateFolds{otherIdx}),:);
        
    mdlClosest = fitcecoc(trainRClosest,trainXClosest);
    mdlOther = fitcecoc(trainROther,trainXOther);
    
    errorPerClosestFold(k-1) = mean(sqrt((mdlClosest.predict(testR)-testX').^2));  
    errorPerClosestFoldSte(k-1) = std(sqrt((mdlClosest.predict(testR)-testX').^2))./sqrt(length(testR));
    
    errorPerOtherFold(k-1) = mean(sqrt((mdlOther.predict(testR)-testX').^2));  
    errorPerOtherFoldSte(k-1) = std(sqrt((mdlOther.predict(testR)-testX').^2))./sqrt(length(testR));
end
toc

% figure;
% boxplot([errorPerClosestFold;errorPerOtherFold]',{'Closest','Other'})
% ylabel('Error')
% set(gcf,'color','white')
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
% title(['Decoding by kfold proximity non perm (',param.traceField,')'],'fontSize',20,'fontWeight','bold')
% ylim([0,3.5])
%save2pdf(['C:\Users\csc\Desktop\caImagg\Graphics\closeKDecodingPerm'])
%%
%%
%Look at the slope of the linear fits of the mean over the kfolds
rCurFolds = zeros([kFolds,size(r,2)]);
rCurFoldsSte = zeros([kFolds,size(r,2)]);
for k = 1:kFolds
    rCurFolds(k,:) = mean(rCur(testIdx(k,:),:));
    rCurFoldsSte(k,:) = std(rCur(testIdx(k,:),:))/sqrt(length(testIdx(k,:)));
end
slopeRCur = zeros([1,size(r,2)]);
model = [];
for celli = 1:size(r,2)
    model(celli).mdl = fitlm(1:kFolds,rCurFolds(:,celli)); 
    slopeRCur(celli) = model(celli).mdl.Coefficients.Estimate(2);
end
figure;
histogram(slopeRCur);grid on;axis square
title('Slope of the linear fits of the mean over the kfolds')
xlabel('Slope')
set(gcf,'color','white')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
save2pdf(['C:\Users\csc\Desktop\caImagg\Graphics\histogramLinearFits'])
%%
figure;
cellsToShow = 40;
subplotN = numSubplots(cellsToShow);
for celli = 1:cellsToShow%size(rPerm,2)
    subplot_er(subplotN(1),subplotN(2),celli);
    errorbar(rCurFolds(:,celli),rCurFoldsSte(:,celli),'lineWidth',2);
    hold on;
    plot(model(celli).mdl.predict([1:length(rCurFolds(:,celli))]'),'lineWidth',2);
    axis square;
    grid on;
    if celli == 1
        legend('Original','Linear fit');
    end
end
annotation('textbox',[0,0,1,1],'string',['Mean value over 10 folds (',param.traceField,')'],'fontSize',20,'fontWeight','bold')
set(gcf,'color','white')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
save2pdf(['C:\Users\csc\Desktop\caImagg\Graphics\meanValOver10Folds'])
%%
%%
%Run svms traiing with data progressively further from test data
N = length(xPerm);
randIdx = randperm(N);
testIdx = randIdx(1:round(N*0.1));
trainIdx = setdiff(randIdx,testIdx);
%Separate the data by proximity in 10 bunchs
trainDist = zeros(size(trainIdx));
trainDistIdx = zeros(size(trainIdx));
for k = 1:length(trainIdx)
    [trainDist(k),trainDistIdx(k)] = min(abs(testIdx-trainIdx(k)));
end
%Sort for distance
nBins = 5;
[~,sortIdx] = sort(trainDist,'ascend');
sortTrainIdx = trainIdx(sortIdx);
distEdges = floor(linspace(0,length(trainDist),nBins+1));
trainIdxCell = cell([1,nBins]);
for k = 2:length(distEdges)
    trainIdxCell{k-1} = sortTrainIdx(((distEdges(k-1))+1):distEdges(k));
end
%%
%Check that its sorted by distance
distCheck = cat(2,trainIdxCell{:});
if length(distCheck)~= length(unique(distCheck))
    error('Non unique val')
end
distCheckMat = cell([1,length(trainIdxCell)]);
for j = 1:length(trainIdxCell)
    for k = 1:length(trainIdxCell{j})
        distCheckMat{j}(k) = min(abs(testIdx-trainIdxCell{j}(k)));
    end
end
figure;hold on;
plotOffset = 0;
for j = 1:length(trainIdxCell)
    plot(plotOffset+(1:length(trainIdxCell{j})),distCheckMat{j});
    plotOffset = plotOffset + length(trainIdxCell{j});
end
%%
%Run the svm
errorPerBin = zeros([1,nBins]);
errorPerBinSte = zeros([1,nBins]);
tic;
for k = 1:nBins
    disp(['Fold ', num2str(k),'/',num2str(nBins)]);
    testX = xCur(testIdx);
    testR = rCur(testIdx,:);
    
    trainX = xCur(trainIdxCell{k});
    trainR = rCur(trainIdxCell{k},:);
        
    mdl = fitcecoc(trainR,trainX);
    errorPerBin(k) = mean(sqrt((mdl.predict(testR)-testX').^2));
    errorPerBinSte(k) = std(sqrt((mdl.predict(testR)-testX').^2))./sqrt(length(testR));
end
toc
disp(errorPerBin)
%%
figure;
subplot(1,2,1)
hold on;
plotOffset = 0;
for j = 1:length(trainIdxCell)
    plot(plotOffset+(1:length(trainIdxCell{j})),distCheckMat{j},'linewidth',2.5);
    plotOffset = plotOffset + length(trainIdxCell{j});
end
xlabel('frames')
ylabel('distance to closest testing data point')
grid on;

subplot(1,2,2)
errorbar(cellfun(@mean,distCheckMat),errorPerBin,errorPerBinSte,'-o','linewidth',2.5);
xlabel('Mean distance');
ylabel('Mean error');axis square;grid on
xticks(cellfun(@mean,distCheckMat))

set(gcf,'color','white')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
annotation('textbox',[0,0,1,1],'string',['Decoding by proximity (',param.traceField,')'],'fontSize',20,'fontWeight','bold')
save2pdf(['C:\Users\csc\Desktop\caImagg\Graphics\proximityDecoding'])
%%
