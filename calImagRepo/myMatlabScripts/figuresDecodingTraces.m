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

for dT = 0.25%[0.05,0.15,0.2,0.25,0.5]

    disp(dT)


    % Set the parameters
    param.dtCamera = 0.05; %period of miniscope (seconds/frame)
    param.dT = dT; %Integration parameter (seconds)
    param.binN = [20,1]; %Number of bins in X Y
    param.minVel = [4, 0]; %Frames under this speed will be discarded (cm/s)
    param.pxPerCm = data.m2019.decA.pxPerCm; %How many px a centimeter is (px/cm)
    param.traceField = 'rawProb'; %Default field to extract traces
    %%
    sessions = data.m2026.decA.chooseSessions('all');
    %Get traces and positions
    output = data.m2019.decA.loadOutput(sessions(1));
    x = output.position(:,1);
    r = output.(param.traceField);
    %%
    [xCur,rCur] = dataAnalysis.curateXandR(x,r,param);
    %%
    %Make folds
    kFolds = 10;
    N = length(xCur);
    foldSize = floor(N/kFolds);
    Nfold = foldSize*kFolds;
    testIdx = zeros([kFolds,foldSize]);
    trainIdx = zeros([kFolds,Nfold-foldSize]);
    for k = 1:kFolds
        testIdx(k,:) = (1+(k-1)*foldSize):(k*foldSize);
        trainIdx(k,:) = setdiff(1:Nfold,testIdx(k,:));
    end
    %%
    %Choose the best fold
    errorRaw = zeros([1,kFolds]);
    errorShuffle = zeros([1,kFolds]);
    for k = 1:kFolds
        disp(['Fold ', num2str(k),'/',num2str(kFolds)]);
        testX = xCur(testIdx(k,:));
        testR = rCur(testIdx(k,:),:);

        trainX = xCur(trainIdx(k,:));
        trainR = rCur(trainIdx(k,:),:);

        mdlRaw = fitcecoc(trainR,trainX);
        errorRaw(k) = mean(sqrt((mdlRaw.predict(testR)-testX).^2))*param.pxPerCm(1);

        trainRShuffle = decoderAnalysis.shuffleKeepingTrials1D(trainR,trainX);
        testRShuffle = decoderAnalysis.shuffleKeepingTrials1D(testR,testX);
        mdlShuff = fitcecoc(trainRShuffle,trainX);
        errorShuffle(k) = mean(sqrt((mdlShuff.predict(testRShuffle)-testX).^2))*param.pxPerCm(1);
    end
    %Take the bigest difference
    [~,fold] = max(errorRaw - errorShuffle);
    %fold = 6
    %%
    %
    %segNum = 10;
    N = [1,10,50,150,300];%round(linspace(1,size(r,2),segNum));
    predRawSmooth = cell([1,length(N)]);
    predShuffleSmooth = cell([1,length(N)]);
    predGTSmooth = cell([1,length(N)]);
    for k = 1:length(N)
        disp([num2str(k),'/',num2str(length(N))])

        xVar = xCur;
        rVar = rCur(:,1:N(k));

        testX = xVar(testIdx(fold,:));
        testR = rVar(testIdx(fold,:),:);

        trainX = xVar(trainIdx(fold,:));
        trainR = rVar(trainIdx(fold,:),:);

        mdlRaw = fitcecoc(trainR,trainX);
        predictRaw = mdlRaw.predict(testR);

        trainRShuffle = decoderAnalysis.shuffleKeepingTrials1D(trainR,trainX);
        testRShuffle = decoderAnalysis.shuffleKeepingTrials1D(testR,testX);
        mdlShuff = fitcecoc(trainRShuffle,trainX);
        predictShuffle = mdlShuff.predict(testRShuffle);

        windowSize = 5;
        b = (1/windowSize)*ones(1,windowSize);
        a = 1;

        predRawSmooth{k} = filtfilt(b,a,predictRaw);
        predShuffleSmooth{k} = filtfilt(b,a,predictShuffle);
        predGTSmooth{k} = filtfilt(b,a,testX);
    end
    %%
    figure;
    seg2show = 1:length(N);
    subplotN = numSubplots(length(seg2show));
    for k = 1:length(seg2show)
        subplot_er(subplotN(1),subplotN(2),k);
        plot(predGTSmooth{seg2show(k)}*param.pxPerCm(1),'linewidth',3);
        hold on;
        plot(predRawSmooth{seg2show(k)}*param.pxPerCm(1),'linewidth',2);
        plot(predShuffleSmooth{seg2show(k)}*param.pxPerCm(1),'linewidth',2);
        if k == 1
            legend('Position','Prediction','Prediction shuffled')
        end
        grid on;
        title([num2str(N(seg2show(k))),' neurons'])
        xlabel('time (frames)')
        ylabel('position (cm)')
    end
    %print(['C:\Users\csc\Desktop\caImagg\Graphics\pairPredictions',num2str(1000*param.dT)],'-dsvg')

    figure;
    seg2show = 1:length(N);
    myCMap = lines;

    subplot_er(2,1,1);
    plot(predGTSmooth{1}*param.pxPerCm(1),'linewidth',3);
    hold on;
    for k = 1:length(seg2show)
        toLeg(k) = plot(predRawSmooth{seg2show(k)}*param.pxPerCm(1),'linewidth',2,'color',(k/length(seg2show))*myCMap(2,:));
        %toLeg(k) = plot(predShuffleSmooth{seg2show(k)}*param.pxPerCm(1),'linewidth',2,'color',(k/length(seg2show))*myCMap(3,:));
        grid on;
        xlabel('time (frames)')
        ylabel('position (cm)')
    end
    legend(toLeg,cellfun(@(x) [num2str(x),' cells'],num2cell(N(seg2show)),'UniformOutput',0))
    title('raw')

    subplot_er(2,1,2);
    plot(predGTSmooth{1}*param.pxPerCm(1),'linewidth',3);
    hold on;
    for k = 1:length(seg2show)
        %toLeg(k) = plot(predRawSmooth{seg2show(k)}*param.pxPerCm(1),'linewidth',2,'color',(k/length(seg2show))*myCMap(2,:));
        toLeg(k) = plot(predShuffleSmooth{seg2show(k)}*param.pxPerCm(1),'linewidth',2,'color',(k/length(seg2show))*myCMap(3,:));
        grid on;
        xlabel('time (frames)')
        ylabel('position (cm)')
    end
    legend(toLeg,cellfun(@(x) [num2str(x),' cells'],num2cell(N(seg2show)),'UniformOutput',0))
    title('shuffle')

    %print(['C:\Users\csc\Desktop\caImagg\Graphics\superposedPredictions',num2str(1000*param.dT)],'-dsvg')
    
end