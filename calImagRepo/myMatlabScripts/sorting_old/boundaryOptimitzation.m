dataDir = 'C:\Users\David\Desktop\idibaps\caImagging\Data\';
animalNames = {'m2029','m2028'};

[myData,allValidCellMax] = getData(dataDir,animalNames);
[myData,allValidCellMax] = removeNansAnd3FromData(myData,allValidCellMax);

fieldNames = fields(myData);
fitFun = @(x) fitnessFun(x,myData,allValidCellMax);
LB = nan([1,2*length(fieldNames)]);
UB = nan([1,2*length(fieldNames)]);
x0 = nan([1,2*length(fieldNames)]);
for k = 1:2:length(fieldNames)*2
    dataField = myData.(fieldNames{ceil(k/2)});
    dataFieldCrop = dataField(abs(dataField) < mean(dataField) + 3*std(dataField));
    LB(k) = min(dataFieldCrop);
    UB(k) = max(dataFieldCrop);
    LB(k+1) = 0;
    UB(k+1) = range(dataFieldCrop);
    x0(k) = mean(dataFieldCrop);
    x0(k+1) = range(dataFieldCrop)/2;
end

varN = 2*length(fieldNames);

%[bound,fval] = ga(fitFun,varN,[],[],[],[],LB,UB);
%[bound,fval] = patternsearch(fitFun,x0,[],[],[],[],LB,UB);
[bound,fval] = particleswarm(fitFun,varN,LB,UB);


%%
centers = bound(1:2:length(bound));
ranges = bound(2:2:length(bound));
uLimits = centers + ranges;
lLimits = centers - ranges;

pred = rulePrdParametric(myData,lLimits,uLimits);
plotconfusion(allValidCellMax,pred')

swarmLimits.uLimits = uLimits;
swarmLimits.lLimits = lLimits;

%%

figure;
dataFields = fields(myData);
dataRange = [0,2;0,4;0,20;0,4;0,1;0.5,1;0,1.3;0,5;0,0.5;0,150;0,2;-0.1,1;0,20;-0.1,1;0,1;0,1;0.5,1;0,1.2;0,5;-0.1,0.6;0,150;0,0.5];
binN = 20;
for k = 1:length(dataFields)
    subplot(ceil(sqrt(length(dataFields))),ceil(sqrt(length(dataFields))),k)
    histData = myData.(dataFields{k});
    histDataCrop = histData(~isnan(histData));
    validCrop = allValidCellMax(~isnan(histData)); 
    %minAxis = min(histDataCrop(abs(histDataCrop) < mean(histDataCrop)+3*std(histDataCrop)));
    %maxAxis = max(histDataCrop(abs(histDataCrop) < mean(histDataCrop)+3*std(histDataCrop)));
    histogram(histDataCrop(~logical(validCrop)),dataRange(k,1):dataRange(k,2)/binN:dataRange(k,2),'FaceColor',[174/255,0,45/255]);
    hold on
    histogram(histDataCrop(logical(validCrop)),[dataRange(k,1):dataRange(k,2)/binN:dataRange(k,2)],'FaceColor',[121/255,180/255,0])
    plot([lLimits(k),lLimits(k)],get(gca,'ylim'),'b')
    plot([uLimits(k),uLimits(k)],get(gca,'ylim'),'r')
    hold off
    %xlim(dataRange(k,:))%([minAxis,maxAxis])
    title(dataFields{k})
end
text(.75,1.25,animalNames,'FontSize',14)
