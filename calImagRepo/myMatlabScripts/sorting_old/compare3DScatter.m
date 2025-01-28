function compare3DScatter(features,myData,varargin)

param.names = features;
param.doMean = false;
param.hideOutliers = false;
param.animalFields = fields(myData);

myCombsFeat = combnk(features,3);
myCombsNames = combnk(param.names,3);

numFeat = size(myCombsFeat,1);



if sum(strcmp(param.animalFields,myCombsFeat{1,1}))
    tempData.Default = myData;
    myData = tempData;
    param.animalFields = fields(myData);
end

if ~isempty(varargin)
    for k = 1:2:length(varargin)
        param.(varargin{k}) = varargin{k+1};
    end  
end

if ~sum(strcmp(fields(myData.(param.animalFields{1})),myCombsFeat{1,1}))
    try
        myData = concatDataPerAnimal(myData,'doMean',param.doMean);
%         animalFields = fields(myData);
    catch
        error([myCombsFeat{1,1},' is a nonexistent feature'])
    end
end
animalFields = param.animalFields;
figure;
myColormap = lines;
outlierSTD = 3;
for i = 1:numFeat
    subplot(ceil(sqrt(numFeat)),ceil(sqrt(numFeat)),i)
    hold on
    axis square
    grid on
    for k = 1:length(animalFields)
        data1 = myData.(animalFields{k}).(myCombsFeat{i,1});
        data2 = myData.(animalFields{k}).(myCombsFeat{i,2});
        data3 = myData.(animalFields{k}).(myCombsFeat{i,3});
        if param.hideOutliers
            outliers = (abs(data1)>abs(nanmean(data1))+outlierSTD*abs(nanstd(data1))) | (abs(data2)>abs(nanmean(data2))+outlierSTD*abs(nanstd(data2))) | (abs(data3)>abs(nanmean(data3))+outlierSTD*abs(nanstd(data3)));
            data1(outliers) = [];
            data2(outliers) = [];
            data3(outliers) = [];
        end
        plot3(data1,data2,data3,'.','Color',myColormap(k,:));
        xlabel(myCombsNames{i,1})
        ylabel(myCombsNames{i,2})
        zlabel(myCombsNames{i,3})
    end
    if length(animalFields) > 1
        h = findobj(gca,'Type','line');
        if length(h) > length(animalFields)
            h = h(1:2:end);
        end
            legend(h(end:-1:1),animalFields)
    end
end

