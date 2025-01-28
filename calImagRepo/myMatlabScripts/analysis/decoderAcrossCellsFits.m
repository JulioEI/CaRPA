%clear all;
%close all;
load 'rawProbSeparate'

param.dtCamera = 0.05; %period of miniscope (seconds/frame)
param.dT = 0.5; %Integration parameter (seconds)
param.binN = [20,1]; %Number of bins in X Y
param.minVel = [4, 0]; %Frames under this speed will be discarded (cm/s)
param.pxPerCm = data.m2024.decA.pxPerCm; %How many px a centimeter is (px/cm)
param.traceField = 'spikeML'; %Default field to extract traces

animals = fieldnames(data);
%animals(7) = [];
%%
dataTypes = {'raw','shuffle'};
yMax = -inf;
sessIdx = 1;
%Run the fits and get the coefficients
for k = 1:length(dataTypes)
    dataType = dataTypes{k};
    
    coefficients.(dataType) = nan([3,length(animals)]);
    meanYVal.(dataType) = nan([1,length(animals)]);
    
    dataInfo.meanSpikeCount.(dataType) = nan([1,length(animals)]);
    dataInfo.fastFramesCount.(dataType) = nan([1,length(animals)]);
    dataInfo.driftVal.(dataType) = nan([1,length(animals)]);
    for animal = 1:length(animals)
        animalData = data.(animals{animal}).decCells;
        cmPerBin = animalData(conf).vec(sessIdx).cmPerBin(1);
        x = animalData(conf).vec(sessIdx).cellVect;
        
        [pos,r] = data.(animals{animal}).decA.getXR('sessions','all','rField','spikeML','dT',0.5);
        v = abs(diff(pos{1}(:,1)./data.(animals{animal}).decA.pxPerCm(1))/data.(animals{animal}).decA.dtCamera);
        [posCut,rCut] = dataAnalysis.removeFramesUnderMinVel(pos{1},r{1},v,data.(animals{animal}).decA.minVel);
        dataInfo.fastFramesCount.(dataType)(animal) = length(posCut);
        dataInfo.meanSpikeCount.(dataType)(animal) = mean(mean(rCut));
        
        slopeRCur = zeros([1,size(r,2)]);
        for celli = 1:size(r,2)
            mdl = fitlm(posCut,rCut(:,celli)); 
            slopeRCur(celli) = mdl.Coefficients.Estimate(2);
        end
        dataInfo.driftVal.(dataType)(animal) = mean(abs(slopeRCur));
        
        switch dataType
            case 'raw'
                y = cmPerBin*mean(mean(animalData(conf).vec(sessIdx).rawData,3));
            case 'shuffle'
                y = cmPerBin*mean(mean(animalData(conf).vec(sessIdx).shuffledData,3));
        end
        meanYVal.(dataType)(animal) = mean(y);
        x = x(6:end);
        y = y(6:end);
        if max(y) > yMax
            yMax = max(y);
        end
        %Polynomial fit 2 order
        options.FunValCheck = 'off';
        if k == 1
            pol2Fn = @(b,x) (b(1)-b(2))*exp(-b(3)*x) + b(2);%@(b,x) b(1)*(x.^-b(2)) + b(3);
            b0Pol2 = [8,1,0];%[30,1,0.2];%[1,1,0];
            coefficients.(dataType)(:,animal) = nlinfit(x,y,pol2Fn,b0Pol2',options);
        else
            b1 = coefficients.(dataTypes{1})(1,animal);
            pol2Fn = @(b,x,b1) (b1-b(1))*exp(-b(2)*x) + b(1);%@(b,x) b(1)*(x.^-b(2)) + b(3);
            b0Pol2 = [1,0];%[30,1,0.2];%[1,1,0];
            coefficients.(dataType)(:,animal) = [b1;nlinfit(x,y,@(b,x) pol2Fn(b,x,b1),b0Pol2',options)];
        end
    end
end
pol2Fn = @(b,x) (b(1)-b(2))*exp(-b(3)*x) + b(2);%@(b,x) b(1)*(x.^-b(2)) + b(3);
%%
% %Look at the mean correlations
correlationVal = zeros([1,length(animals)]);
corrData = load('correlationDataMLSPIKES');
dataInfo.rawMinusShuffleCorr.raw = nan([1,length(animals)]);
dataInfo.rawMinusShuffleCorr.shuffle = nan([1,length(animals)]);
dataInfo.sumCorr.raw = nan([1,length(animals)]);
dataInfo.sumCorr.shuffle = nan([1,length(animals)]);
for animal = 1:length(animals)
    disp([num2str(animal),'/',num2str(length(animals))]);
%     [pos,r] = data.(animals{animal}).decA.getXR('sessions','all','rField','spikeML','dT',0.5);
%     %Curate
%     [xCur,rCur] = dataAnalysis.curateXandR(pos{1},r{1},param);
%     xCur = xCur(:,1);
%     %Shuffle
%     xShuf = xCur;
%     rShuf = decoderAnalysis.shuffleKeepingTrials1D(rCur,xCur);
%     %Trial shuffle
%     %[rTrialShuff,xTrialShuff] = decoderAnalysis.shuffleSegments(rCur,xCur,1);
%     %Randomly permute
%     randPermIdx = randperm(length(xCur));%1:length(xCur);%
%     xPerm = xCur(randPermIdx);
%     rPerm = rCur(randPermIdx,:);
%     [correlations.raw.mean,correlations.raw.ste,correlations.raw.x,correlations.raw.mat] = computeCorrelations(xCur,rCur);
%     [correlations.shuffle.mean,correlations.shuffle.ste,correlations.shuffle.x,correlations.shuffle.mat] = computeCorrelations(xShuf,rShuf);
    rawCorr = (corrData.correlations.(animals{animal}).raw.mean);
    shuffleCorr = (corrData.correlations.(animals{animal}).shuffle.mean);

    dataInfo.rawMinusShuffleCorr.raw(animal) = nansum(rawCorr - shuffleCorr);
    dataInfo.rawMinusShuffleCorr.shuffle(animal) = nansum(rawCorr - shuffleCorr);
    dataInfo.sumCorr.raw(animal) = nansum(rawCorr);
    dataInfo.sumCorr.shuffle(animal) = nansum(shuffleCorr);
    
    %correlationVal(animal) = nansum(abs(corrData.correlations.(animals{animal}).raw.mean - corrData.correlations.(animals{animal}).shuffle.mean));
    
end
% = correlationVal;
%dataInfo.rawMinusShuffleCorr.shuffle = correlationVal;
%%
%Look at correlation between spike mean and mean dec value
figure;
dataInfoF = fields(dataInfo);
nSub = numSubplots(length(dataInfoF));
for k = 1:length(dataInfoF)
    subplot_er(nSub(1),nSub(2),k);
    transparent(1) = plot(dataInfo.(dataInfoF{k}).raw,meanYVal.raw,'color','white');
    text(dataInfo.(dataInfoF{k}).raw, meanYVal.raw, cellfun(@(x) x(4:end),animals,'UniformOutput',0), 'HorizontalAlignment','center', 'VerticalAlignment','baseline','Color','blue','FontWeight','bold')
    hold on
    transparent(2) = plot(dataInfo.(dataInfoF{k}).shuffle,meanYVal.shuffle,'color','white');
    text(dataInfo.(dataInfoF{k}).shuffle, meanYVal.shuffle, cellfun(@(x) x(4:end),animals,'UniformOutput',0), 'HorizontalAlignment','center', 'VerticalAlignment','baseline','Color','red','FontWeight','bold')
    
    for i = 1:length(transparent)
        transparent(i).Color(4) = 0;
    end
    
    mdlRaw = fitlm(dataInfo.(dataInfoF{k}).raw,meanYVal.raw,'RobustOpts','on');
    mdlShuffle = fitlm(dataInfo.(dataInfoF{k}).shuffle,meanYVal.shuffle,'RobustOpts','on');
    
    semiTransparent(1) = plot(dataInfo.(dataInfoF{k}).raw',mdlRaw.predict(dataInfo.(dataInfoF{k}).raw'),'b','lineWidth',2);
    semiTransparent(2) = plot(dataInfo.(dataInfoF{k}).shuffle',mdlShuffle.predict(dataInfo.(dataInfoF{k}).shuffle'),'r','lineWidth',2);
     
    for i = 1:length(semiTransparent)
        semiTransparent(i).Color(4) = .2;
    end   
    
    dummy(1) = plot(nan,nan,'.b','MarkerSize',20);
    dummy(2) = plot(nan,nan,'.r','MarkerSize',20);
    legend(dummy,['raw',' r^2 ',num2str(round(100*mdlRaw.Rsquared.Adjusted)/100)], ['shuffle',' r^2 ',num2str(round(100*mdlShuffle.Rsquared.Adjusted)/100)])
    xlabel(dataInfoF{k})
    ylabel('Mean decoding error')
    grid on;
    axis square;
end
set(gcf,'color','white')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
%export_fig asdf.jpg -m3
%%
%Get the axis
decAxis = [0,yMax];
coeffAxis = [min(0,min([min(coefficients.raw,[],2),min(coefficients.shuffle,[],2)],[],2)),max([max(coefficients.raw,[],2),max(coefficients.shuffle,[],2)],[],2)];%[min([min(coefficients.raw,[],2),min(coefficients.shuffle,[],2)],[],2),max([max(coefficients.raw,[],2),max(coefficients.shuffle,[],2)],[],2)];

%Look at the coefficients
figure;
myCMFull = jet;
myCM = myCMFull(round(linspace(1,length(myCMFull)-length(myCMFull)/length(animals),length(animals))),:);
nPoints = 1000;
ax = [];
for k = 1:length(dataTypes)
    ax(k) = subplot_er(2,(size(coefficients.raw,1)+1),1+(k-1)*(size(coefficients.raw,1)+1));
    hold on;
    dataType = dataTypes{k};
    toLeg = [];
    for animal = 1:length(animals)
        animalData = data.(animals{animal}).decCells;
        cmPerBin = animalData(conf).vec(sessIdx).cmPerBin(1);
        x = animalData(conf).vec(sessIdx).cellVect;
        switch dataType
            case 'raw'
                y = cmPerBin*mean(mean(animalData(conf).vec(sessIdx).rawData,3));
            case 'shuffle'
                y = cmPerBin*mean(mean(animalData(conf).vec(sessIdx).shuffledData,3));
        end
        toLeg = [toLeg,plot(nan,nan,'.','LineWidth',2,'MarkerSize',12,'color',myCM(animal,:))];
        plot(x,(y),'.-','LineWidth',2,'MarkerSize',12,'color',myCM(animal,:))
        xTest = linspace(x(6),x(end)*3,3*nPoints);
        plot(xTest,(pol2Fn(coefficients.(dataType)(:,animal)',xTest)),'LineWidth',2,'color',0.7*myCM(animal,:))
    end
    grid on;
    %axis square;
    %set(gca, 'YScale', 'log')
    xlim([0,700])
    title(dataType)
    xlabel('Cell num')
    ylabel('Log of absolute decoding error')
    
    lgnd = legend(toLeg, cellfun(@(x) x(4:end),animals,'UniformOutput',0));
    set(lgnd,'color','none');
    set(lgnd,'Box','off');
end

linkaxes;
for k = 1:size(coefficients.raw,1)
    for j = 1:length(dataTypes)
        subplot_er(2,(size(coefficients.raw,1)+1),(k+1)+(j-1)*(size(coefficients.raw,1)+1));
        hold on
        barInfo = coefficients.(dataTypes{j})(k,:)';
        for i = 1:size(barInfo,1)
            myBar = bar(i,barInfo(i,1));
            set(myBar,'FaceColor',myCM(i,:));
        end
        xlabel('Animal (20xx)')
        xticks(1:length(animals))
        xticklabels(cellfun(@(x) x(4:end),animals,'UniformOutput',0))
        ylabel('MLE')
        ylim(coeffAxis(k,:));
        switch k
            case 1
                title(['a',' (',dataTypes{j},')'])
            case 2
                title(['b',' (',dataTypes{j},')'])
            case 3
                title(['c',' (',dataTypes{j},')'])
        end
        axis square;
        grid on;
    end
end
%funName = 'ax^{-b}+c';
funname = '(a-b)*e^{-cx} + b';
annotation('textbox',[0,0,1,1],'String',funname,'FitBoxToText','off','FontSize',20,'FontWeight','bold','LineStyle','none');
set(gcf,'color','white')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
%export_fig asdf.jpg -m3
%save2pdf(['C:\Users\csc\Desktop\caImagg\Graphics\powerFits'])
%%
%Look at the coefficients
figure;
ax = [];
for k = 1:length(dataTypes)
    ax(k) = subplot_er(2,2,k)%(size(coefficients.raw,1)+1),1+(k-1)*(size(coefficients.raw,1)+1))%subplot_er(2,(size(coefficients.raw,1)+1),1+(k-1)*(size(coefficients.raw,1)+1));
    hold on;
    dataType = dataTypes{k};
    toLeg = [];
    for animal = 1:length(animals)
        animalData = data.(animals{animal}).decCells;
        cmPerBin = animalData(conf).vec(sessIdx).cmPerBin(1);
        x = animalData(conf).vec(sessIdx).cellVect;
        switch dataType
            case 'raw'
                y = cmPerBin*mean(mean(animalData(conf).vec(sessIdx).rawData,3));
            case 'shuffle'
                y = cmPerBin*mean(mean(animalData(conf).vec(sessIdx).shuffledData,3));
        end
        toLeg = [toLeg,plot(nan,nan,'.','LineWidth',2,'MarkerSize',12,'color',myCM(animal,:))];
        plot(x,(y),'.-','LineWidth',2,'MarkerSize',12,'color',myCM(animal,:))
        xTest = linspace(x(6),x(end)*3,3*nPoints);
        plot(xTest,(pol2Fn(coefficients.(dataType)(:,animal)',xTest)),'LineWidth',2,'color',0.7*myCM(animal,:))
    end
    grid on;
    %axis square;
    %set(gca, 'YScale', 'log')
    ylim([0,15])
    xlim([0,700])
    title(dataType)
    xlabel('Cell num')
    ylabel('Log of absolute decoding error')
    
    lgnd = legend(toLeg, cellfun(@(x) x(4:end),animals,'UniformOutput',0));
    set(lgnd,'color','none');
    set(lgnd,'Box','off');
end

linkaxes;


for k = 2:3
    subplot_er(2,2,k+1);
    coeffRaw = coefficients.raw(k,:)';
    coeffShuffle = coefficients.shuffle(k,:)';
    offset = 0;
    for coeff = 1:length(coeffRaw)
        plot([0,1],coeff*offset + [coeffRaw(coeff),coeffShuffle(coeff)],'.-','linewidth',2,'markersize',30,'color',myCM(coeff,:))
        hold on;
    end
    grid on
    [~,pval] = ttest(coeffRaw - coeffShuffle);
    switch k
        case 2
            title(['b ',num2str(pval)])
        case 3
            title(['c ',num2str(pval)])
    end
    ylim([0,max(max(coefficients.raw(k,:)))])
    xlim([-0.5,1.5])
end

funname = '(a-b)*e^{-cx} + b';
annotation('textbox',[0,0,1,1],'String',funname,'FitBoxToText','off','FontSize',20,'FontWeight','bold','LineStyle','none');
set(gcf,'color','white')
set(gcf,'units','normalized','outerposition',[0 0 1 1])

%%
%Show all cells
figure;
nPoints = 1000;

hold on;
toLeg = [];
lastDiffVal = zeros([1,length(animals)]);
for animal = 1:length(animals)
    
    animalData = data.(animals{animal}).decCells;
    x = animalData(conf).vec(sessIdx).cellVect;
    xTest = linspace(x(1),x(end)*3,3*nPoints);
    
    rawMeanVal = nanmean(nanmean(animalData(conf).vec(sessIdx).rawData,3),1);
    shuffleMeanVal = nanmean(nanmean(animalData(conf).vec(sessIdx).shuffledData,3),1);
    cmPerBin = animalData(conf).vec(sessIdx).cmPerBin(1);
    
    %plot(x,cmPerBin*(rawMeanVal-shuffleMeanVal),'LineWidth',2,'color',myCM(animal,:))
    
    y = (pol2Fn(coefficients.raw(:,animal)',xTest)-pol2Fn(coefficients.shuffle(:,animal)',xTest));
    lastDiffVal(animal) = y(end);
    plot(xTest,cmPerBin*y,'-','LineWidth',2,'MarkerSize',12,'color',0.7*myCM(animal,:))
    
    toLeg = [toLeg,plot(nan,nan,'.','LineWidth',2,'MarkerSize',12,'color',myCM(animal,:))];
    %h = shadedErrorBar((toPlot.cellVec),(cmPerBin*toPlot.rawMeanVal),(cmPerBin*toPlot.rawSteVal),'lineprops',{'-','color',[178, 83, 0]/255,'linewidth',1},'patchSaturation',0.2);
    %h = shadedErrorBar((toPlot.cellVec),(cmPerBin*toPlot.shuffleMeanVal),(cmPerBin*toPlot.shuffleSteVal),'lineprops',{'-','color',[70, 160, 178]/255,'linewidth',1},'patchSaturation',0.2);
end
lgnd = legend(toLeg, cellfun(@(x) x(4:end),animals,'UniformOutput',0));
set(lgnd,'color','none');
set(lgnd,'Box','off');
set(lgnd,'color','none');
set(lgnd,'Box','off');
grid on;
axis square;
xlabel('Cell num')
ylabel('Error raw - error shuffle')
set(gcf,'color','white')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
%export_fig asdf.jpg -m3
%%

figure;
hold on
transparent(1) = plot(dataInfo.rawMinusShuffleCorr.raw,lastDiffVal,'color','white');
text(dataInfo.rawMinusShuffleCorr.raw, lastDiffVal, cellfun(@(x) x(4:end),animals,'UniformOutput',0), 'HorizontalAlignment','center', 'VerticalAlignment','baseline','Color','blue','FontWeight','bold')
for i = 1:length(transparent)
    transparent(i).Color(4) = 0;
end
mdl = fitlm(dataInfo.rawMinusShuffleCorr.raw,lastDiffVal,'RobustOpts','on');
semiTransparent(1) = plot(dataInfo.rawMinusShuffleCorr.raw',mdl.predict(dataInfo.rawMinusShuffleCorr.raw'),'b','lineWidth',2);
for i = 1:length(semiTransparent)
    semiTransparent(i).Color(4) = .2;
end
dummy(1) = plot(nan,nan,'.b','MarkerSize',20);
legend(dummy,['r^2' ,num2str(round(100*mdl.Rsquared.Adjusted)/100)])
xlabel('Raw correlations - shuffle correlations')
ylabel('Raw error - shuffle error')
grid on;
axis square;
set(gcf,'color','white')
set(gcf,'units','normalized','outerposition',[0 0 1 1])

%export_fig asdf.jpg -m3

%%
%Show means 
rawY.data = cell([1,length(animals)]);
shuffleY.data = cell([1,length(animals)]);
rawY.fit = cell([1,length(animals)]);
shuffleY.fit = cell([1,length(animals)]);

for k = 1:length(dataTypes)
    dataType = dataTypes{k};
    for animal = 1:length(animals)        
        animalData = data.(animals{animal}).decCells;
        cmPerBin = animalData(conf).vec(sessIdx).cmPerBin(1);
        x = animalData(conf).vec(sessIdx).cellVect;
        switch dataType
            case 'raw'
                y = cmPerBin*mean(mean(animalData(conf).vec(sessIdx).rawData,3));
            case 'shuffle'
                y = cmPerBin*mean(mean(animalData(conf).vec(sessIdx).shuffledData,3));
        end
        xTest = linspace(x(6),x(end)*10,10*nPoints);
        switch dataType
            case 'raw'
                rawY.data{animal} = y;
                rawX.data{animal} = x;
                rawY.fit{animal} = pol2Fn(coefficients.(dataType)(:,animal)',xTest);
                rawX.fit{animal} = xTest;
            case 'shuffle'
                shuffleY.data{animal} = y;
                shuffleX.data{animal} = x;
                shuffleY.fit{animal} = pol2Fn(coefficients.(dataType)(:,animal)',xTest);  
                shuffleX.fit{animal} = xTest;                  
        end
    end
end
%%
figure;
typeData = {'data','fit'};
binWidth = [[0:10:40],[50:50:200],270,500];
for k = 1:length(typeData)
    subplot(1,2,k)
    [x,y,ste,N] = meanForDifferentX(rawX.(typeData{k}),rawY.(typeData{k}),binWidth);
    disp(N)
    h = shadedErrorBar(x,y,ste,'lineprops',{'.-','color',[178, 83, 0]/255,'linewidth',1},'patchSaturation',0.2);h.mainLine.DisplayName = 'Raw';
    [x,y,ste,N] = meanForDifferentX(shuffleX.(typeData{k}),shuffleY.(typeData{k}),binWidth);
    h = shadedErrorBar(x,y,ste,'lineprops',{'.-','color',[70, 160, 178]/255,'linewidth',1},'patchSaturation',0.2);h.mainLine.DisplayName = 'Shuffle';
    title((typeData{k}))
    xlim([0,500])
    ylim([0,15])
    grid on;
end
xlabel('Number of cells')
ylabel('AbsoluteError(cm)')
legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'),'Location','best');

