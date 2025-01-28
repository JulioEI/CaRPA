function feature2DScatter(features,myData,validCell,showData,showBoundary,showEllipse,varargin)
%This plots features and valid/nonvalid points in 3d. Last parameter is
%optional and specifies the label names. If not set, the variable names
%will be used instead.

names = features;
if ~isempty(varargin)
    if length(names) ~= length(features)
        warning('DIFERENT NUM OF FEATURES AND NAMES, IGNORING NAMES')
    else
        names = varargin{1};
    end
end

myCombsFeat = combnk(features,2);
myCombsNames = combnk(names,2);

numFeat = size(myCombsFeat,1);
figure;
for i = 1:numFeat
    subplot(ceil(sqrt(numFeat)),ceil(sqrt(numFeat)),i)
    hold on
    axis square
    grid on
    
    data1 = myData.(myCombsFeat{i,1});
    data2 = myData.(myCombsFeat{i,2});

    idx = validCell == 1;  
    if showData  
        plot(data1(idx),data2(idx),'.','color',[121/255,180/255,0],'markersize',1)
        hold on
        plot(data1(~idx),data2(~idx),'.','color',[174/255,0,45/255],'markersize',1)
    end
    xlabel(myCombsNames{i,1})
    ylabel(myCombsNames{i,2})
    
   if showBoundary   
       tresholdData1 = nanmean(data1) + 3*nanstd(data1);
       tresholdData2 = nanmean(data2) + 3*nanstd(data2);
       tresholdIdx = data1 < tresholdData1 & data2 < tresholdData2;
       data1Robust = data1(tresholdIdx);
       data2Robust = data2(tresholdIdx);
       validCellRobust = validCell(tresholdIdx);
       B2 = mnrfit([data1Robust,data2Robust],categorical(validCellRobust));
       X = linspace(min(data1Robust),max(data1Robust),50);
       Y = (B2(1)+B2(2).*X)./(-B2(3));
       plot(X,Y,'color',[84/255,220/255,196/255],'linewidth',2)
   end
    
   if showEllipse
       [r_ellipse,X0,Y0] = fitEllipseToData([data1(idx),data2(idx)],.95);
       plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-','Color',[121/255,180/255,0],'LineWidth',2);
       [r_ellipse,X0,Y0] = fitEllipseToData([data1(~idx),data2(~idx)],.95);
       plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-','Color',[174/255,0,45/255],'LineWidth',2);
   end
end

