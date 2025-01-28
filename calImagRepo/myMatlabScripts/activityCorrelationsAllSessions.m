%Load animals
%clear all
home
animals = {
    'C:\Users\csc\Desktop\caImagg\AllData\Mouse2022'
    'C:\Users\csc\Desktop\caImagg\AllData\Mouse2011'
    'C:\Users\csc\Desktop\caImagg\AllData\Mouse2021'
    'C:\Users\csc\Desktop\caImagg\AllData\Mouse2024'
    'C:\Users\csc\Desktop\caImagg\AllData\Mouse2026'
    'C:\Users\csc\Desktop\caImagg\AllData\Mouse2019'
    'C:\Users\csc\Desktop\caImagg\AllData\Mouse2012'
    'C:\Users\csc\Desktop\caImagg\AllData\Mouse2028'
    'C:\Users\csc\Desktop\caImagg\AllData\Mouse2025'
    'C:\Users\csc\Desktop\caImagg\AllData\Mouse2023'
    'C:\Users\csc\Desktop\caImagg\AllData\Mouse2010'
    'C:\Users\csc\Desktop\caImagg\AllData\Mouse2029'};
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
param.traceField = 'spikeML'; %Default field to extract traces
    
x = cell([length(animals),1]);
r = cell([length(animals),1]);
disp('Collecting data...')
for animal = 1:length(animals)
    disp([num2str(animal),'/',num2str(length(animals))])
    sessions = data.(animals{animal}).decA.chooseSessions('all');
    x{animal} = cell([length(sessions),1]);
    r{animal} = cell([length(sessions),1]);
    %Get traces and positions
    for sess = 1:length(sessions)
        output = data.(animals{animal}).decA.loadOutput(sessions(sess));
        x{animal}{sess} = output.position;
        r{animal}{sess} = output.(param.traceField);
    end
end
%%
% %Check traces
% figure;
% for animal = 1:length(animals)
%     subplotN = numSubplots(length(x{animal}));
%     sessions = data.(animals{animal}).decA.chooseSessions('all');
%     for sess = 1:length(x{animal})
%         subplot_er(subplotN(1),subplotN(2),sess);
%         plot(x{animal}{sess})
%         title(sessions(sess).folderName(12:(end-13)))
%     end
%     annotation('textbox',[0,0,1,1],'String',animals{animal},'FitBoxToText','off','FontSize',20,'FontWeight','bold','LineStyle','none');
%     pause;
%     clf;
% end

for dT = 0.5%[0.05,0.15,0.2,0.25,0.5]
    
    % Set the parameters
    param.dtCamera = 0.05; %period of miniscope (seconds/frame)
    param.dT = dT;%0.05; %Integration parameter (seconds)
    param.binN = [20,1]; %Number of bins in X Y
    param.minVel = [4, 0]; %Frames under this speed will be discarded (cm/s)
    
    disp(['Finding noise correlations for dt ',num2str(param.dT),' ...'])
    for animal = 1:length(animals)
        disp([num2str(animal),'/',num2str(length(animals))])
        param.pxPerCm = data.(animals{animal}).decA.pxPerCm; %How many px a centimeter is (px/cm)
        sessions = data.(animals{animal}).decA.chooseSessions('all');
        for sess = 1:length(sessions)
            
            %Curate
            [xCur,rCur,cmPerBin] = dataAnalysis.curateXandR(x{animal}{sess},r{animal}{sess},param);
            xCur = xCur(:,1);

            %Shuffle
            xShuf = xCur;
            rShuf = decoderAnalysis.shuffleKeepingTrials1D(rCur,xCur);
            
            %Find correlations
            [correlations.(['dT',num2str(round(1000*dT))]).(animals{animal})(sess).raw.mean,correlations.(['dT',num2str(round(1000*dT))]).(animals{animal})(sess).raw.ste,correlations.(['dT',num2str(round(1000*dT))]).(animals{animal})(sess).raw.x,correlations.(['dT',num2str(round(1000*dT))]).(animals{animal})(sess).raw.mat] = computeCorrelations(xCur,rCur);
            [correlations.(['dT',num2str(round(1000*dT))]).(animals{animal})(sess).shuffle.mean,correlations.(['dT',num2str(round(1000*dT))]).(animals{animal})(sess).shuffle.ste,correlations.(['dT',num2str(round(1000*dT))]).(animals{animal})(sess).shuffle.x,correlations.(['dT',num2str(round(1000*dT))]).(animals{animal})(sess).shuffle.mat] = computeCorrelations(xShuf,rShuf);

            [correlationsPairs.(['dT',num2str(round(1000*dT))]).(animals{animal})(sess).raw.signalCorrelation,correlationsPairs.(['dT',num2str(round(1000*dT))]).(animals{animal})(sess).raw.trialByTrialCorrelation] = computeCorrelationPairs(xCur,rCur);
            [correlationsPairs.(['dT',num2str(round(1000*dT))]).(animals{animal})(sess).shuffle.signalCorrelation,correlationsPairs.(['dT',num2str(round(1000*dT))]).(animals{animal})(sess).shuffle.trialByTrialCorrelation] = computeCorrelationPairs(xShuf,rShuf);
               
        end
    end
end
    %%
%     figure;
%     nSub = numSubplots(length(animals));
%     for animal = 1:length(animals)
%         ax = subplot_er(nSub(1),nSub(2),animal);
%         hold on;
%         myCmap = lines;
% 
%         scRaw = correlationsPairs.(animals{animal}).raw.signalCorrelation;
%         tcRaw = correlationsPairs.(animals{animal}).raw.trialByTrialCorrelation;
%         scShuffle = correlationsPairs.(animals{animal}).shuffle.signalCorrelation;
%         tcShuffle = correlationsPairs.(animals{animal}).shuffle.trialByTrialCorrelation;
%         plot(scRaw,tcRaw,'.','MarkerSize',1,'color',myCmap(1,:));
%         plot(scShuffle,tcShuffle,'.','MarkerSize',1,'color',myCmap(2,:));
% 
%         mdlRaw = fitlm(scRaw,tcRaw,'RobustOpts','on');
%         mdlShuffle = fitlm(scShuffle,tcShuffle,'RobustOpts','on');
% 
%         plot((min(scRaw):0.01:max(scRaw))',mdlRaw.predict((min(scRaw):0.01:max(scRaw))'),'color',myCmap(1,:),'lineWidth',2);
%         plot((min(scShuffle):0.01:max(scShuffle))',mdlShuffle.predict((min(scShuffle):0.01:max(scShuffle))'),'color',myCmap(2,:),'lineWidth',2);
% 
% 
% 
%         ylabel('Noise correlation')
%         if animal == 1
%             legend('Raw','Shuffle');
%         end
%         grid on;
%         title((animals{animal}))
%         xlabel('Signal correlation')
%         ylim([-1,1])
%     end
%     linkaxes;
%     set(gcf,'color','white')
%     set(gcf,'units','normalized','outerposition',[0 0 1 1])
%     annotation('textbox',[0,0,1,1],'String',[param.traceField,num2str(1000*param.dT),'msec'],'FitBoxToText','off','FontSize',20,'FontWeight','bold','LineStyle','none');
%     print(['C:\Users\csc\Desktop\caImagg\Graphics\pairwiseNoiseCorrelations',num2str(1000*param.dT),'msec'],'-dsvg')
%     %export_fig asdf.jpg -m3
% 
%     %%
%myDT = 'dT50';  
myDT = {'dT500'};%{'dT50','dT150','dT200','dT250'};
dtnum = .5%[0.05,0.15,0.2,0.25];
for k = 1:length(myDT)
    for animal = 1:length(animals)
        sessions = data.(animals{animal}).decA.chooseSessions('all');
        figure;
        nSub = numSubplots(length(sessions));
        for sess = 1:length(sessions)
            ax = subplot_er(nSub(1),nSub(2),sess);

            myCmap = lines;
            fNames = fields(correlations.(myDT{k}).(animals{animal})(1));

            for f = 1:length(fNames)
                h = shadedErrorBar(correlations.(myDT{k}).(animals{animal})(sess).(fNames{f}).x,correlations.(myDT{k}).(animals{animal})(sess).(fNames{f}).mean,1.96*correlations.(myDT{k}).(animals{animal})(sess).(fNames{f}).ste,'lineprops',{'-o','linewidth',2,'color',myCmap(f,:)},'patchSaturation',0.2);h.mainLine.DisplayName = (fNames{f});
            end
            ylabel('Noise correlation')
            xticks(correlations.(myDT{k}).(animals{animal})(sess).(fNames{1}).x)
            myLabels = linspace(-1,1,length(correlations.(myDT{k}).(animals{animal})(sess).(fNames{1}).x));
            xticklabels(round(myLabels*10)/10)
            if sess == 1
                legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'),'Location','best');
            end
            grid on;
            xlabel('Signal correlation')
            title(sessions(sess).folderName(12:(end-13)))
        end
        linkaxes;
        annotation('textbox',[0,0,1,1],'String',[animals{animal},'-',num2str(1000*dtnum(k)),'ms'],'FitBoxToText','off','FontSize',20,'FontWeight','bold','LineStyle','none');
        set(gcf,'color','white')
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        save2pdf(['C:\Users\csc\Desktop\caImagg\Graphics\corrAcrossSess',animals{animal},'_',num2str(1000*dtnum(k)),'ms'])
        %pause;
        %clf;
    end
end
close all
joinPDFs('C:\Users\csc\Desktop\caImagg\Graphics\corrAcrossSess.pdf','C:\Users\csc\Desktop\caImagg\Graphics')

%annotation('textbox',[0,0,1,1],'String',[param.traceField,num2str(1000*param.dT),'msec'],'FitBoxToText','off','FontSize',20,'FontWeight','bold','LineStyle','none');
%print(['C:\Users\csc\Desktop\caImagg\Graphics\meanNoiseCorrelations',num2str(1000*param.dT),'msec'],'-dsvg')
%export_fig asdf.jpg -m3

%%
%myDT = 'dT50';  
myDT = {'dT50','dT150','dT200','dT250'};
dtnum = [0.05,0.15,0.2,0.25];
for k = 1:length(myDT)
    figure;
    nSub = numSubplots(length(animals));
    for animal = 1:length(animals)
        
        sessions = data.(animals{animal}).decA.chooseSessions('all');
        ax = subplot_er(nSub(1),nSub(2),animal);
        
        for sess = 1:length(sessions)
            

            myCmap = lines;
            fNames = fields(correlations.(myDT{k}).(animals{animal})(1));

            for f = 1:length(fNames)
                h = shadedErrorBar(correlations.(myDT{k}).(animals{animal})(1).(fNames{f}).x,correlations.(myDT{k}).(animals{animal})(sess).(fNames{f}).mean,1.96*correlations.(myDT{k}).(animals{animal})(sess).(fNames{f}).ste,'lineprops',{'-o','linewidth',2,'color',myCmap(f,:)},'patchSaturation',0.2);h.mainLine.DisplayName = (fNames{f});
            end
            ylabel('Noise correlation')
            xticks(correlations.(myDT{k}).(animals{animal})(sess).(fNames{1}).x)
            myLabels = linspace(-1,1,length(correlations.(myDT{k}).(animals{animal})(sess).(fNames{1}).x));
            xticklabels(round(myLabels*10)/10)
            if sess == 1
                legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'),'Location','best');
            end
            grid on;
            xlabel('Signal correlation')
            title(sessions(sess).folderName(12:(end-13)))
        end
        linkaxes;
        annotation('textbox',[0,0,1,1],'String',[animals{animal},'-',num2str(1000*dtnum(k)),'ms'],'FitBoxToText','off','FontSize',20,'FontWeight','bold','LineStyle','none');
        set(gcf,'color','white')
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        %save2pdf(['C:\Users\csc\Desktop\caImagg\Graphics\corrAcrossSess',animals{animal},'_',num2str(1000*dtnum(k)),'ms'])
        %pause;
        %clf;
    end
end
close all
%joinPDFs('C:\Users\csc\Desktop\caImagg\Graphics\corrAcrossSess.pdf','C:\Users\csc\Desktop\caImagg\Graphics')

%annotation('textbox',[0,0,1,1],'String',[param.traceField,num2str(1000*param.dT),'msec'],'FitBoxToText','off','FontSize',20,'FontWeight','bold','LineStyle','none');
%print(['C:\Users\csc\Desktop\caImagg\Graphics\meanNoiseCorrelations',num2str(1000*param.dT),'msec'],'-dsvg')
%export_fig asdf.jpg -m3

%%
