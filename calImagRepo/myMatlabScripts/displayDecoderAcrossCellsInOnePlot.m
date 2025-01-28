joinDir = 0;
conf = 1;
sessIdx = 1;
animalVec = 1:length(animals);
figure;
nSub = numSubplots(length(animalVec));
for k = animalVec%1:length(animals)
    
    ax = subplot_er(nSub(1),nSub(2),k);
    
    animalData = data.(['m',mouseNames{k}]).decCells;   

    sessionN = length(animalData(conf).vec);

    toPlot.cellVec = animalData(conf).vec(sessIdx).cellVect;
    cmPerBin = animalData(conf).vec(sessIdx).cmPerBin(1);

    toPlot.rawMeanVal = nanmean(nanmean(animalData(conf).vec(sessIdx).rawData,3),1);
    toPlot.rawSteVal = 1.96*nanstd(nanmean(animalData(conf).vec(sessIdx).rawData,3),[],1)/sqrt(size(nanmean(animalData(conf).vec(sessIdx).rawData,3),1));
    toPlot.shuffleMeanVal = nanmean(nanmean(animalData(conf).vec(sessIdx).shuffledData,3),1);
    toPlot.shuffleSteVal = 1.96*nanstd(nanmean(animalData(conf).vec(sessIdx).shuffledData,3),[],1)/sqrt(size(nanmean(animalData(conf).vec(sessIdx).shuffledData,3),1));
    toPlot.permutedMeanVal = nanmean(nanmean(animalData(conf).vec(sessIdx).permData,3),1);
    toPlot.permutedSteVal = 1.96*nanstd(nanmean(animalData(conf).vec(sessIdx).permData,3),[],1)/sqrt(size(nanmean(animalData(conf).vec(sessIdx).shuffledData,3),1));

    h = shadedErrorBar((toPlot.cellVec),(cmPerBin*toPlot.rawMeanVal),(cmPerBin*toPlot.rawSteVal),'lineprops',{'.-','color',[178, 83, 0]/255,'linewidth',1},'patchSaturation',0.2);h.mainLine.DisplayName = 'Raw';
    h = shadedErrorBar((toPlot.cellVec),(cmPerBin*toPlot.permutedMeanVal),(cmPerBin*toPlot.permutedSteVal),'lineprops',{'.-','color',[211, 0, 255]/255,'linewidth',1},'patchSaturation',0.2);h.mainLine.DisplayName = 'Permuted';
    h = shadedErrorBar((toPlot.cellVec),(cmPerBin*toPlot.shuffleMeanVal),(cmPerBin*toPlot.shuffleSteVal),'lineprops',{'.-','color',[70, 160, 178]/255,'linewidth',1},'patchSaturation',0.2);h.mainLine.DisplayName = 'Shuffle';

    clear toPlot
    xlabel('Number of cells')
    ylabel('AbsoluteError(cm)')
    if k == 1
        legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'),'Location','best');
    end
    title(['m',mouseNames{k}])
    %linkaxes(ax,'y')
    ylim([0,70]);
    grid on;%title(myTitle);
end

%figText = {['m',mouseNames{k}],animalData(conf).decoder,animalData(conf).responseType,[num2str(animalData(conf).timeWindow),' sec']};
    
set(gcf,'color','white')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
export_fig asdf.jpg -m3
%exportFig(frame2im(getframe(gcf)),'newImage.jpg','jpg','Quality',100)
