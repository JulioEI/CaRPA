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
animals = fieldnames(data);
%%

for dT = [0.05,0.15,0.2,0.25,0.5]

    disp(dT)

    for animal = 1:length(animals)
        disp([num2str(animal),'/',num2str(length(animals))])
        % Set the parameters
        param.dtCamera = 0.05; %period of miniscope (seconds/frame)
        param.dT = dT;%0.05; %Integration parameter (seconds)
        param.binN = [20,1]; %Number of bins in X Y
        param.minVel = [4, 0]; %Frames under this speed will be discarded (cm/s)
        param.pxPerCm = data.(animals{animal}).decA.pxPerCm; %How many px a centimeter is (px/cm)
        param.traceField = 'rawProb'; %Default field to extract traces
        %%
        sessions = data.(animals{animal}).decA.chooseSessions('all');
        %Get traces and positions
        output = data.(animals{animal}).decA.loadOutput(sessions(1));
        x = output.position;
        r = output.(param.traceField);
        %%
        %Curate
        [xCur,rCur,cmPerBin] = dataAnalysis.curateXandR(x,r,param);
        xCur = xCur(:,1);

        %Shuffle
        xShuf = xCur;
        rShuf = decoderAnalysis.shuffleKeepingTrials1D(rCur,xCur);

        %Trial shuffle
        %[rTrialShuff,xTrialShuff] = decoderAnalysis.shuffleSegments(rCur,xCur,1);

        %Randomly permute
        randPermIdx = randperm(length(xCur));%1:length(xCur);%
        xPerm = xCur(randPermIdx);
        rPerm = rCur(randPermIdx,:);
        %%
        [correlations.(animals{animal}).raw.mean,correlations.(animals{animal}).raw.ste,correlations.(animals{animal}).raw.x,correlations.(animals{animal}).raw.mat] = computeCorrelations(xCur,rCur);
        [correlations.(animals{animal}).shuffle.mean,correlations.(animals{animal}).shuffle.ste,correlations.(animals{animal}).shuffle.x,correlations.(animals{animal}).shuffle.mat] = computeCorrelations(xShuf,rShuf);

        [correlationsPairs.(animals{animal}).raw.signalCorrelation,correlationsPairs.(animals{animal}).raw.trialByTrialCorrelation] = computeCorrelationPairs(xCur,rCur);
        [correlationsPairs.(animals{animal}).shuffle.signalCorrelation,correlationsPairs.(animals{animal}).shuffle.trialByTrialCorrelation] = computeCorrelationPairs(xShuf,rShuf);
           %%
    end
    %%
    figure;
    nSub = numSubplots(length(animals));
    for animal = 1:length(animals)
        ax = subplot_er(nSub(1),nSub(2),animal);
        hold on;
        myCmap = lines;

        scRaw = correlationsPairs.(animals{animal}).raw.signalCorrelation;
        tcRaw = correlationsPairs.(animals{animal}).raw.trialByTrialCorrelation;
        scShuffle = correlationsPairs.(animals{animal}).shuffle.signalCorrelation;
        tcShuffle = correlationsPairs.(animals{animal}).shuffle.trialByTrialCorrelation;
        plot(scRaw,tcRaw,'.','MarkerSize',1,'color',myCmap(1,:));
        plot(scShuffle,tcShuffle,'.','MarkerSize',1,'color',myCmap(2,:));

        mdlRaw = fitlm(scRaw,tcRaw,'RobustOpts','on');
        mdlShuffle = fitlm(scShuffle,tcShuffle,'RobustOpts','on');

        plot((min(scRaw):0.01:max(scRaw))',mdlRaw.predict((min(scRaw):0.01:max(scRaw))'),'color',myCmap(1,:),'lineWidth',2);
        plot((min(scShuffle):0.01:max(scShuffle))',mdlShuffle.predict((min(scShuffle):0.01:max(scShuffle))'),'color',myCmap(2,:),'lineWidth',2);



        ylabel('Noise correlation')
        if animal == 1
            legend('Raw','Shuffle');
        end
        grid on;
        title((animals{animal}))
        xlabel('Signal correlation')
        ylim([-1,1])
    end
    linkaxes;
    set(gcf,'color','white')
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    annotation('textbox',[0,0,1,1],'String',[param.traceField,num2str(1000*param.dT),'msec'],'FitBoxToText','off','FontSize',20,'FontWeight','bold','LineStyle','none');
    print(['C:\Users\csc\Desktop\caImagg\Graphics\pairwiseNoiseCorrelations',num2str(1000*param.dT),'msec'],'-dsvg')
    %export_fig asdf.jpg -m3

    %%
    figure;
    nSub = numSubplots(length(animals));
    for animal = 1:length(animals)
        ax = subplot_er(nSub(1),nSub(2),animal);

        myCmap = lines;
        fNames = fields(correlations.(animals{animal}));

        %subplot(4,1,[1,2,3])
        for f = 1:length(fNames)
            h = shadedErrorBar(correlations.(animals{animal}).(fNames{f}).x,correlations.(animals{animal}).(fNames{f}).mean,1.96*correlations.(animals{animal}).(fNames{f}).ste,'lineprops',{'-o','linewidth',2,'color',myCmap(f,:)},'patchSaturation',0.2);h.mainLine.DisplayName = (fNames{f});
        end
        ylabel('Noise correlation')
    %     myXTicks = cellfun(@num2str,num2cell([-1,correlations.(animals{animal}).(fNames{1}).x]),'UniformOutput',0);
    %     myXTicks2 = cell([1,length(myXTicks)-1]);
    %     for k = 1:length(myXTicks2)
    %         myXTicks2{k} = [myXTicks{k}, ' to ', myXTicks{k+1}];
    %     end
        xticks(correlations.(animals{animal}).(fNames{1}).x)
        myLabels = linspace(-1,1,length(correlations.(animals{animal}).(fNames{1}).x));
        xticklabels(round(myLabels*10)/10)
        %xticklabels(myXTicks2)
        if animal == 1
            legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'),'Location','best');
        end
        grid on;
        title((animals{animal}))
        %axis square;
        % subplot(4,1,4)
        % barData = [];
        % for f = 1:length(fNames)
        %     barData = [barData,cellfun(@length,correlations.(fNames{f}).mat.PFcorrPerPair)/sum(cellfun(@length,correlations.(fNames{f}).mat.PFcorrPerPair))];
        % end
        % h = bar(barData);
        % for k = 1:length(h)
        %     h(k).FaceColor = myCmap(k,:);
        % end
        % xticklabels(myXTicks2)
        % ylabel('% of pairs in that bin')
        xlabel('Signal correlation')
        % %legend(fNames{:})
        % grid on;
    end
    linkaxes;
    set(gcf,'color','white')
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    annotation('textbox',[0,0,1,1],'String',[param.traceField,num2str(1000*param.dT),'msec'],'FitBoxToText','off','FontSize',20,'FontWeight','bold','LineStyle','none');
    print(['C:\Users\csc\Desktop\caImagg\Graphics\meanNoiseCorrelations',num2str(1000*param.dT),'msec'],'-dsvg')
    %export_fig asdf.jpg -m3
    
    
    
    
end