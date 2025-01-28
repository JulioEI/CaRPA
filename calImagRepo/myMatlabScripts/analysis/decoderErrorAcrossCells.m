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
% sessions = {
%     '20150326'
%     '20150207'
%     '20150326'
%     '20150311'
%     '20150228'
%     '20150228'
%     '20150118'
%     '20150228'
%     '20150228'
%     '20150326'
%     '20150128'
%     '20150311'};
pxPerCm = cell([1,length(animals)]);
load pxPerCm
data = [];
animalVec = 1:length(animals);%[4,6];
%%
%Create decoder analysis object for each animal
mouseNames = cell([1,length(animals)]);
for k = animalVec%1:length(animals)
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
timeWindow = .5;%[0.45,0.5,0.6,0.8];%%0.25;%[.05,.5];
responseType = {'spikeML'};%{'rawProb'};%spikeDeconv','spikeML','rawProb'};
directionality = 0;
decoder = {'SVM'};%,'NB','NBnormal'};
%%
ctTotal = length(animals)*length(timeWindow)*length(responseType)*length(directionality)*length(decoder);
for k = animalVec%1:length(animals)
sessionVect = 'all';%sessionVect = {7,{1}};%{3,{1},round(3+(length(data.(['m',mouseNames{k}]).decA.dayStruct)-3)/2),{1},length(data.(['m',mouseNames{k}]).decA.dayStruct),{1}};
ct = 1;
disp(repmat('%',[3,100]))
disp(['Animal ',num2str(k),'/',num2str(length(animals))])
    for tw = timeWindow
        for rt = responseType
            for dir = directionality
                for deco = decoder
                    data.(['m',mouseNames{k}]).decCells(ct).timeWindow = tw;
                    data.(['m',mouseNames{k}]).decCells(ct).responseType = rt{1};
                    data.(['m',mouseNames{k}]).decCells(ct).directionality = dir;
                    data.(['m',mouseNames{k}]).decCells(ct).decoder = deco{1};
                    data.(['m',mouseNames{k}]).decCells(ct).sessions = sessionVect;
                    data.(['m',mouseNames{k}]).decCells(ct).vec = ...
                        data.(['m',mouseNames{k}]).decA.decoderAcrossCellsAdaptative('dT',tw,'rField',rt{1},'dir',dir,'decoder',deco{1},'sessions',sessionVect);
                    ct = ct + 1;
                end
            end
        end
    end
end
%%
% for k = 1:length(animals)
%     myCmPerBin = data.(['m',mouseNames{k}]).decA.getCmPerBin('dT',tw,'rField',rt{1},'dir',dir,'decoder',deco{1},'sessions',sessionVect);
%     for c = 1:length(myCmPerBin)
%         data.(['m',mouseNames{k}]).decCells.vec(c).cmPerBin = myCmPerBin{c};
%     end
% end

%%
joinDir = 0;
if ~joinDir;figureRatio = 3/4;end
for k = animalVec%1:length(animals)
    
    animalData = data.(['m',mouseNames{k}]).decCells;
    %[xSess,~] = data.(['m',mouseNames{k}]).decA.getXR('sessions',data.(['m',mouseNames{k}]).decCells(1).sessions,'rField','rawProb');%Assuming all conf have the same sessions
    if joinDir
        figIDTemp = 2*((1:length(animalData)/2) - 1)+1;
        figID = [figIDTemp,figIDTemp+1];
        figIDK = 0.5;
    else
        figID = 1:length(animalData);
    end
    for conf = figID
        if joinDir
            figure(100*(k-1)+ceil(figIDK));
            figIDK = figIDK + 0.5;
        else
            figure(100*(k-1)+conf);clf;
        end
        
        sessionN = length(animalData(conf).vec);
        for sessIdx = 1:sessionN
            
            if ~joinDir               
                subplotRow = round(sqrt(sessionN/figureRatio));
                subplotCol = round(figureRatio*subplotRow);
                if subplotCol*subplotRow < sessionN;subplotRow = subplotRow + 1;end
                subplot(subplotRow,subplotCol,sessIdx)
            end

            toPlot.cellVec = animalData(conf).vec(sessIdx).cellVect;
            cmPerBin = animalData(conf).vec(sessIdx).cmPerBin(1);
            if animalData(conf).directionality
                
                if joinDir
                    subplot(2,sessionN,sessIdx)
                end
                
                toPlot.rawMeanValL = nanmean(nanmean(animalData(conf).vec(sessIdx).rawData.L,3));
                toPlot.rawMeanValR = nanmean(nanmean(animalData(conf).vec(sessIdx).rawData.R,3));
                toPlot.rawSteValL = 1.96*nanstd(nanmean(animalData(conf).vec(sessIdx).rawData.L,3))/sqrt(size(nanmean(animalData(conf).vec(sessIdx).rawData.L,3),1));
                toPlot.rawSteValR = 1.96*nanstd(nanmean(animalData(conf).vec(sessIdx).rawData.R,3))/sqrt(size(nanmean(animalData(conf).vec(sessIdx).rawData.R,3),1));
                toPlot.shuffleMeanValL = nanmean(nanmean(animalData(conf).vec(sessIdx).shuffledData.L,3));
                toPlot.shuffleMeanValR = nanmean(nanmean(animalData(conf).vec(sessIdx).shuffledData.R,3));
                toPlot.shuffleSteValL = 1.96*nanstd(nanmean(animalData(conf).vec(sessIdx).shuffledData.L,3))/sqrt(size(nanmean(animalData(conf).vec(sessIdx).shuffledData.L,3),1));
                toPlot.shuffleSteValR = 1.96*nanstd(nanmean(animalData(conf).vec(sessIdx).shuffledData.R,3))/sqrt(size(nanmean(animalData(conf).vec(sessIdx).shuffledData.R,3),1));  

                h = shadedErrorBar(toPlot.cellVec,cmPerBin*toPlot.rawMeanValL,cmPerBin*toPlot.rawSteValL,'lineprops',{'-o','color',[178, 83, 0]/255,'linewidth',1},'patchSaturation',0.2);h.mainLine.DisplayName = 'Raw L';
                h = shadedErrorBar(toPlot.cellVec,cmPerBin*toPlot.rawMeanValR,cmPerBin*toPlot.rawSteValR,'lineprops',{'-o','color',[255, 167, 90]/255,'linewidth',1},'patchSaturation',0.2);h.mainLine.DisplayName = 'Raw R';
                h = shadedErrorBar(toPlot.cellVec,cmPerBin*toPlot.shuffleMeanValL,cmPerBin*toPlot.shuffleSteValL,'lineprops',{'-o','color',[70, 160, 178]/255,'linewidth',1},'patchSaturation',0.2);h.mainLine.DisplayName = 'Shuffle L';
                h = shadedErrorBar(toPlot.cellVec,cmPerBin*toPlot.shuffleMeanValR,cmPerBin*toPlot.shuffleSteValR,'lineprops',{'-o','color',[0, 211, 255]/255,'linewidth',1},'patchSaturation',0.2);h.mainLine.DisplayName = 'Shuffle R';
            else

                if joinDir
                    subplot(2,sessionN,sessIdx + sessionN)
                end
                                
                toPlot.rawMeanVal = nanmean(nanmean(animalData(conf).vec(sessIdx).rawData,3),1);
                toPlot.rawSteVal = 1.96*nanstd(nanmean(animalData(conf).vec(sessIdx).rawData,3),[],1)/sqrt(size(nanmean(animalData(conf).vec(sessIdx).rawData,3),1));
                toPlot.shuffleMeanVal = nanmean(nanmean(animalData(conf).vec(sessIdx).shuffledData,3),1);
                toPlot.shuffleSteVal = 1.96*nanstd(nanmean(animalData(conf).vec(sessIdx).shuffledData,3),[],1)/sqrt(size(nanmean(animalData(conf).vec(sessIdx).shuffledData,3),1));
                %toPlot.permutedMeanVal = nanmean(nanmean(animalData(conf).vec(sessIdx).permData,3),1);
                %toPlot.permutedSteVal = 1.96*nanstd(nanmean(animalData(conf).vec(sessIdx).permData,3),[],1)/sqrt(size(nanmean(animalData(conf).vec(sessIdx).shuffledData,3),1));

                h = shadedErrorBar((toPlot.cellVec),(cmPerBin*toPlot.rawMeanVal),(cmPerBin*toPlot.rawSteVal),'lineprops',{'.-','color',[178, 83, 0]/255,'linewidth',1},'patchSaturation',0.2);h.mainLine.DisplayName = 'Raw';
                %h = shadedErrorBar((toPlot.cellVec),(cmPerBin*toPlot.permutedMeanVal),(cmPerBin*toPlot.permutedSteVal),'lineprops',{'.-','color',[0, 211, 255]/255,'linewidth',1},'patchSaturation',0.2);h.mainLine.DisplayName = 'Permuted';
                h = shadedErrorBar((toPlot.cellVec),(cmPerBin*toPlot.shuffleMeanVal),(cmPerBin*toPlot.shuffleSteVal),'lineprops',{'.-','color',[70, 160, 178]/255,'linewidth',1},'patchSaturation',0.2);h.mainLine.DisplayName = 'Shuffle';
            end
            clear toPlot
            xlabel('Number of cells')
            ylabel('AbsoluteError(cm)')
            legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'),'Location','best');
            
            
            %Plot velocity limits
%             [xSess,~] = data.(['m',mouseNames{k}]).decA.getXR;%Assuming all conf have the same sessions
%             for sessIdx = 1:length(xSess)
            %if ~joinDir || ceil(mod(figIDK,1))
%                 x = xSess{sessIdx}./data.(['m',mouseNames{k}]).decA.pxPerCm;
%                 xVel = dataAnalysis.removeFramesUnderMinVel(x,x,data.(['m',mouseNames{k}]).decA.minVel,data.(['m',mouseNames{k}]).decA.dtCamera);
%                 [~,xI] = dataAnalysis.integrateInTime(xVel,xVel,data.(['m',mouseNames{k}]).decCells(conf).timeWindow,data.(['m',mouseNames{k}]).decA.dtCamera);
%                 v = abs(diff(xI)/data.(['m',mouseNames{k}]).decCells(conf).timeWindow);
%                 velQuantiles = [quantile(v(:,2),0.05),quantile(v(:,1),0.95)]*data.(['m',mouseNames{k}]).decCells(conf).timeWindow;
%                 hold on;
%                 plot(animalData(conf).vec(sessIdx).cellVect,velQuantiles(1)*ones(size(animalData(conf).vec(sessIdx).cellVect)),'k--');
%                 plot(animalData(conf).vec(sessIdx).cellVect,velQuantiles(2)*ones(size(animalData(conf).vec(sessIdx).cellVect)),'k--');
%                 hold off;
            %end
%             histogram(v(:,1),100);
%             pause;
%             end
           
            %ylim([0,cmPerBin*10]);
            grid on;%title(myTitle);
        end
        figText = {['m',mouseNames{k}],animalData(conf).decoder,animalData(conf).responseType,[num2str(animalData(conf).timeWindow),' sec']};
        if ~joinDir || ceil(mod(figIDK,1))
            figText = ['m',mouseNames{k},' ',animalData(conf).decoder,' ',animalData(conf).responseType,' ',num2str(1000*animalData(conf).timeWindow),'msec'];
            annotation('textbox',[0.3,0,1,1],'String',figText,'FitBoxToText','on','FontSize',12);
            save2pdf(['C:\Users\csc\Desktop\caImagg\Graphics\rawProbSeparate_','m',mouseNames{k},animalData(conf).decoder,'_',animalData(conf).responseType,'_',num2str(1000*animalData(conf).timeWindow)])
        end
    end
end
%%
joinPDFs('C:\Users\csc\Desktop\caImagg\Graphics\rawProbSeparate.pdf','C:\Users\csc\Desktop\caImagg\Graphics')

%%
%%
%TTEST
% numComparisons = 10;
% [pValShuf,pValRaw] = deal(zeros([1,numComparisons+1]));
% meanData = nanmean(animalData(conf).vec(sessIdx).shuffledData,3);
% for t = 0:numComparisons
%     [~,pValShuf(t+1)] = ttest(meanData(:,end-t),meanData(:,end-t-1));
% end
% hShuf= pValShuf < 0.05;
% table(hShuf',pValShuf','VariableNames',{'SHUFFLE','pVal'})
% meanData = nanmean(animalData(conf).vec(sessIdx).rawData,3);
% for t = 0:numComparisons
%     [~,pValRaw(t+1)] = ttest(meanData(:,end-t),meanData(:,end-t-1));
% end
% hRaw = pValRaw < 0.05;
% table(hRaw',pValRaw','VariableNames',{'RAW','pVal'})