doMean = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%Histograms%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataPerAnimalConcat = concatDataPerAnimal(dataPerAnimal,'doMean',doMean);
% dataPerAnimalConcat.m28ICA500.globalSNR =  dataPerAnimalConcat.m28EMit125.globalSNR+1;%+ (0.5*rand([1,length(dataPerAnimalConcat.m28EMit125.globalSNR)])-0.5);
% dataPerAnimalConcat.m28EMit50.expRatio(abs(dataPerAnimalConcat.m28EMit50.expRatio)>3) = [];
% dataPerAnimalConcat.m28EMit125.expRatio(abs(dataPerAnimalConcat.m28EMit125.expRatio)>3) = [];
%%
animalFields = {'m28EM'};%'m28EMit50','m28EMit125','m28EMit200','m28EMit500','m28ICA500'}';%fields(dataPerAnimalConcat);
data = dataPerAnimalConcat.(animalFields{1});
dataFields = fields(data);
binN = 20;
figure;
for k = 1:length(dataFields)
    minAxis = [];
    maxAxis = [];
    for animal = animalFields'
        
        data = dataPerAnimalConcat.(animal{1});      
        
        subplot(ceil(sqrt(length(dataFields))),ceil(sqrt(length(dataFields))),k)
        histData = data.(dataFields{k});
        
        histDataCrop = histData(~isnan(histData) & ~isinf(histData));
        histDataCrop = histDataCrop(abs(histDataCrop) < mean(histDataCrop)+3*std(histDataCrop));
        minAxisAnimal = min(histDataCrop(abs(histDataCrop) < mean(histDataCrop)+3*std(histDataCrop)));
        maxAxisAnimal = max(histDataCrop(abs(histDataCrop) < mean(histDataCrop)+3*std(histDataCrop)));
        minAxis = [minAxis,minAxisAnimal];
        maxAxis = [maxAxis,maxAxisAnimal];
        histogram(histDataCrop,minAxisAnimal:maxAxisAnimal/binN:maxAxisAnimal)
        hold on;
    end
    xlim([min(minAxis),max(maxAxis)])
    title(dataFields{k})
%     text(1.25,0,animal{1},'fontsize',14)
end
legend(animalFields)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Cell comparison%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellsPerAnimalCell = {};
animalFields = {'m28EMit50','m28EMit100','m28EMit125','m28EMit200','m28EMit500'}';
for animal = animalFields'
    dayData = dataPerAnimal.(animal{1});
    dayFields = fields(dayData);
    cellsPerDay = [];
    for day = dayFields'
        featData = dayData.(day{1});
        featFields = fields(featData); 
        cellsPerDay = [cellsPerDay,length(featData.(featFields{1}))];
    end
    cellsPerAnimalCell = [cellsPerAnimalCell;cellsPerDay];
end
daysInAnimal = cellfun(@length,cellsPerAnimalCell);
cellsPerAnimal = nan([length(cellsPerAnimalCell),max(daysInAnimal)]);
for k = 1:length(cellsPerAnimalCell)
    cellsPerAnimal(k,1:daysInAnimal(k)) = cellsPerAnimalCell{k};
end
figure;
ax=subplot(1,2,1);bar(cellsPerAnimal);xticklabels(animalFields);ax.XTickLabelRotation=45;title('total cells per day');grid on
ax=subplot(1,2,2);bar(nanmean(cellsPerAnimal,2));xticklabels(animalFields);ax.XTickLabelRotation=45;title('average cells per day');grid on
%%
%Remove non valid cells
animalFields = fields(predPerAnimalAll);
clear dataPerAnimalValid
for animalField = animalFields'
    dayFields = fields(dataPerAnimal.(animalField{1}));
    for dayField = dayFields'
        try
            validCells = logical(predPerAnimalAll.(animalField{1}).(dayField{1}));     
            featFields = fields(dataPerAnimal.(animalField{1}).(dayField{1}));
            for featField = featFields'
                dataTemp = dataPerAnimal.(animalField{1}).(dayField{1}).(featField{1});
                dataValid = dataTemp(validCells);
                dataPerAnimalValid.(animalField{1}).(dayField{1}).(featField{1}) = dataValid;
                clear dataTemp dataValid
            end
        catch
            disp(['predictions not found for ',animalField{1}, ' ',dayField{1},' discarding day...'])
        end   
    end
end
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%Scatters%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compare2DScatter({'imMovCorr','CellArea','CellSolidity','CellEccentricity','peakSNR','expRatio'},dataPerAnimalValid,true,true,'doMean',true,'hideOutliers',true)%,'animalFields',{'m28EMit50','m28EMit125','m28EMit200','m28EMit500','m28ICA500'})

%%
compare3DScatter({'peakSNR','CellSolidity','imMovCorr'},dataPerAnimalValid,'doMean',true,'hideOutliers',true)%,'animalFields',{'m28EMit50','m28EMit125','m28EMit200','m28EMit500','m28ICA500'})