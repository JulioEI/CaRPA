function feature3DScatter(features,myData,validCell,showBoundary,varargin)
%This plots features and valid/nonvalid points in 2d. Last parameter is
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

myCombsFeat = combnk(features,3);
myCombsNames = combnk(names,3);

numFeat = size(myCombsFeat,1);
figure;
for i = 1:numFeat
    subplot(ceil(sqrt(numFeat)),ceil(sqrt(numFeat)),i)
  
    data1 = myData.(myCombsFeat{i,1});
    data2 = myData.(myCombsFeat{i,2});
    data3 = myData.(myCombsFeat{i,3});
    
    idx = validCell == 1;    
    plot3(data1(idx),data2(idx),data3(idx),'.','color',[121/255,180/255,0])
    axis square
    grid on

    hold on
    plot3(data1(~idx),data2(~idx),data3(~idx),'.','color',[174/255,0,45/255])
    axis square
    grid on

    if showBoundary
       tresholdData1 = nanmean(data1) + 3*nanstd(data1);
       tresholdData2 = nanmean(data2) + 3*nanstd(data2);
       tresholdData3 = nanmean(data3) + 3*nanstd(data3);
       tresholdIdx = data1 < tresholdData1 & data2 < tresholdData2 & data3 < tresholdData3;
       data1Robust = data1(tresholdIdx);
       data2Robust = data2(tresholdIdx);
       data3Robust = data3(tresholdIdx);
       validCellRobust = validCell(tresholdIdx);         
       B2 = mnrfit([data1Robust,data2Robust,data3Robust],categorical(validCellRobust));
       [X,Y] = meshgrid(linspace(min(data1Robust),max(data1Robust),50),linspace(min(data2Robust),max(data2Robust),100));
       Z = (B2(1)+B2(2).*X+B2(3).*Y)./(-B2(4));
       mesh(X,Y,Z);colormap([84/255,220/255,196/255])
    end
    
    xlabel(myCombsNames{i,1})
    ylabel(myCombsNames{i,2})
    zlabel(myCombsNames{i,3})
end

end


