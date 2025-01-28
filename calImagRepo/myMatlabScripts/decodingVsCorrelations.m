%load('myStuff')

load spikeCorrelations
load quickDataSpikes
myCorr = correlations.dT500;
myData = data;
%clearvars -except myCorr myData
%%
animals = fields(myCorr);

%Find the sum of noise correlations for each animal and session
noiseCorr = cell([1,length(animals)]);
for animal = 1:length(animals)
    for sess = 1:length(myCorr.(animals{animal}))
        noiseCorr{animal} = [noiseCorr{animal},sum(abs(myCorr.(animals{animal})(sess).raw.mean))];
    end
end

%%
%Find the mean decoding error
decErr.raw = cell([1,length(animals)]);
decErr.shuffle = cell([1,length(animals)]);
decErr.rawMShuffle = cell([1,length(animals)]);
goodSess = cell([1,length(animals)]);

figure;hold on
myCMFull = hsv;
myCM = myCMFull(round(linspace(1,length(myCMFull)-length(myCMFull)/length(animals),length(animals))),:);
for animal = 1:length(animals)
    vecs = myData.(animals{animal}).decCells.vec;
    animalData = myData.(animals{animal}).decCells;
    for sess = 1:length(vecs)
        cmPerBin = animalData.vec(sess).cmPerBin(1);
        
        decErr.raw{animal} = [decErr.raw{animal},mean(cmPerBin*mean(vecs(sess).rawData,3))];
        decErr.shuffle{animal} = [decErr.shuffle{animal},mean(cmPerBin*mean(vecs(sess).shuffledData,3))];
        decErr.rawMShuffle{animal} = [decErr.rawMShuffle{animal},sum(cmPerBin*mean(vecs(sess).rawData,3)-cmPerBin*mean(vecs(sess).shuffledData,3))];
        
        %plot(mean(vecs(sess).shuffledData,3),'color',myCM(animal,:));
        
        tmp = mean(vecs(sess).shuffledData,3);
        if (tmp(end)-tmp(1))./tmp(1) > -0.3
            plot(mean(vecs(sess).shuffledData,3),'--r');
            %plot(mean(vecs(sess).rawData,3),'-.r');
        else
            goodSess{animal} = [goodSess{animal},sess];
            plot(mean(vecs(sess).shuffledData,3),'--b');
            %plot(mean(vecs(sess).rawData,3),'-.b');
        end
    end
end
%%
%Remove bad sessions
for animal = 1:length(animals)
    decErr.raw{animal} = decErr.raw{animal}(goodSess{animal});
    decErr.shuffle{animal} = decErr.shuffle{animal}(goodSess{animal});
    decErr.rawMShuffle{animal} = decErr.rawMShuffle{animal}(goodSess{animal});
    noiseCorr{animal} = noiseCorr{animal}(goodSess{animal});
end
%%

allDecErr.raw = cat(2,decErr.raw{:});
allDecErr.shuffle = cat(2,decErr.shuffle{:});
allNoiseCorr = cat(2,noiseCorr{:});

myCMFull = hsv;
myCM = myCMFull(round(linspace(1,length(myCMFull)-length(myCMFull)/length(animals),length(animals))),:);

figure;hold on
for animal = 1:length(animals)
    %Look at correlation between spike mean and mean dec value
    transparent(1) = plot(noiseCorr{animal},decErr.raw{animal},'x','color',myCM(animal,:));
    transparent(2) = plot(noiseCorr{animal},decErr.shuffle{animal},'d','color',myCM(animal,:));

%     for i = 1:length(transparent)
%         transparent(i).Color(4) = 0;
%     end
end

mdlRaw = fitlm(allNoiseCorr,allDecErr.raw,'RobustOpts','on');
mdlShuffle = fitlm(allNoiseCorr,allDecErr.shuffle,'RobustOpts','on');

% xTest = linspace(min(allNoiseCorr),max(allNoiseCorr),1000);
% semiTransparent(1) = plot(xTest',mdlRaw.predict(xTest'),'k','lineWidth',2);
% semiTransparent(2) = plot(xTest'',mdlShuffle.predict(xTest'),'k','lineWidth',2);

% for i = 1:length(semiTransparent)
%     semiTransparent(i).Color(4) = .2;
% end   
 
dummy(1) = plot(nan,nan,'xk');
dummy(2) = plot(nan,nan,'ok');
legend(dummy,'raw','shuffle')
%legend(dummy,['raw',' r^2 ',num2str(round(100*mdlRaw.Rsquared.Adjusted)/100)], ['shuffle',' r^2 ',num2str(round(100*mdlShuffle.Rsquared.Adjusted)/100)])

xlabel('Noise corr')
ylabel('Mean decoding error')
grid on;
axis square;
set(gcf,'color','white')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
%export_fig asdf.jpg -m3

%%

myCMFull = hsv;
myCM = myCMFull(round(linspace(1,length(myCMFull)-length(myCMFull)/length(animals),length(animals))),:);

figure;
for animal = 1:length(animals)
    transparent(1) = plot(noiseCorr{animal},decErr.rawMShuffle{animal},'.','color',myCM(animal,:),'MarkerSize',20);
    hold on
end

xlabel('Noise corr')
ylabel('Mean raw - shuffle error')
grid on;
axis square;
set(gcf,'color','white')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
%export_fig asdf.jpg -m3

%%
%Look at how diff depends on several other things
y = decErr.rawMShuffle;
x = zeros([1,length(y)]);
clear dataInfo
dataInfo.fastFramesCount = cell([1,length(animals)]);
dataInfo.meanSpikeCount = cell([1,length(animals)]);
dataInfo.driftVal = cell([1,length(animals)]);
dataInfo.maxNeurons = cell([1,length(animals)]);

for animal = 1:length(animals)
    disp([num2str(animal),'/',num2str(length(animals))])
    vecs = myData.(animals{animal}).decCells.vec;
    animalData = myData.(animals{animal}).decCells;
    allSessions = {myData.(animals{animal}).decA.dayStruct.sessions};allSessions = cat(2,allSessions{:});
    for sess = 1:length(vecs)
        cmPerBin = animalData.vec(sess).cmPerBin(1);
        x = animalData.vec(sess).cellVect;
        %[pos,r] = myData.(animals{animal}).decA.getXR('sessions',sess,'rField','spikeML','dT',0.5);
        out = load([allSessions(sess).rootFolder,filesep,allSessions(sess).folderName,filesep,allSessions(sess).tracesEventsFileName]);
        pos = out.tracesEvents.position;
        r = out.tracesEvents.rawProb;
        v = abs(diff(pos(:,1)./myData.(animals{animal}).decA.pxPerCm(1))/myData.(animals{animal}).decA.dtCamera);
        [posCut,rCut] = dataAnalysis.removeFramesUnderMinVel(pos,r,v,myData.(animals{animal}).decA.minVel);

        dataInfo.fastFramesCount{animal} = [dataInfo.fastFramesCount{animal},length(posCut)];
        dataInfo.meanSpikeCount{animal} = [dataInfo.meanSpikeCount{animal},mean(mean(rCut))];

        slopeRCur = zeros([1,size(r,2)]);
        for celli = 1:size(r,2)
            mdl = fitlm(posCut,rCut(:,celli)); 
            slopeRCur(celli) = mdl.Coefficients.Estimate(2);
        end
        dataInfo.driftVal{animal} = [dataInfo.driftVal{animal},mean(abs(slopeRCur))];
        dataInfo.maxNeurons{animal} = [dataInfo.maxNeurons{animal},max(x)];
    end
end
%%
%Remove bad sessions
for animal = 1:length(animals)
    dataInfo.driftVal{animal} = dataInfo.driftVal{animal}(goodSess{animal});
    dataInfo.fastFramesCount{animal} = dataInfo.fastFramesCount{animal}(goodSess{animal});
    dataInfo.maxNeurons{animal} = dataInfo.maxNeurons{animal}(goodSess{animal}); 
    dataInfo.meanSpikeCount{animal} = dataInfo.meanSpikeCount{animal}(goodSess{animal});
end
%%
myCMFull = hsv;
myCM = myCMFull(round(linspace(1,length(myCMFull)-length(myCMFull)/length(animals),length(animals))),:);

figure;
dataField = fields(dataInfo);
for k = 1:length(dataField)
    subplot_er(1,length(dataField),k)
    for animal = 1:length(animals)
        plot(dataInfo.(dataField{k}){animal},decErr.rawMShuffle{animal},'.','color',myCM(animal,:),'MarkerSize',20);
        hold on
    end

    xlabel(dataField{k})
    ylabel('Mean raw - shuffle error')
    grid on;
    axis square;
end
set(gcf,'color','white')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
%export_fig asdf.jpg -m3

