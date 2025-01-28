function compare2DScatter(features,myData,showData,showEllipse,varargin)

param.names = features;
param.doMean = false;
param.hideOutliers = false;
param.animalFields = fields(myData);

myCombsFeat = combnk(features,2);
myCombsNames = combnk(param.names,2);

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
% 
%         if strcmp(myCombsFeat{i,1},'globalSNR')
%             data1 = (data1-min(data1))/(max(data1)-min(data1));
%             disp('normalizing data')
%         end
%         if strcmp(myCombsFeat{i,2},'globalSNR')
%             data2 = (data2-min(data2))/(max(data2)-min(data2));
%             disp('normalizing data')
%         end
                
        if param.hideOutliers
            outliers = (abs(data1)>abs(nanmean(data1))+outlierSTD*abs(nanstd(data1))) | (abs(data2)>abs(nanmean(data2))+outlierSTD*abs(nanstd(data2)));
            data1(outliers) = [];
            data2(outliers) = [];
        end
        if showData  
            plot(data1,data2,'.','Color',myColormap(k,:));
        end
        xlabel(myCombsNames{i,1})
        ylabel(myCombsNames{i,2})

       if showEllipse
           [r_ellipse,X0,Y0] = fitEllipseToData([data1(~isnan(data1)&~isnan(data2)&~isinf(data1)&~isinf(data2)),data2(~isnan(data1)&~isnan(data2)&~isinf(data1)&~isinf(data2))],.9);
           plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-','Color',myColormap(k,:),'LineWidth',2);
       end
    end
end
if length(animalFields) > 1
    h = findobj(gca,'Type','line');
    if length(h) > length(animalFields)
        h = h(1:2:end);
    end
        legend(h(end:-1:1),animalFields)
end
