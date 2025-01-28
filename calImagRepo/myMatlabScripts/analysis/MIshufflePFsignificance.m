[x,r] = PC.getXR('sessions',4,'rField','spikeDeconv');
x = x{1};
r = r{1};

x = x./PC.pxPerCm;
%Get absolute velocity
v = abs(diff(x)/PC.dtCamera);
%Get frames under min velocity
framesUnderMinVel = unique([find(v(:,1) < PC.minVel(1)); find(v(:,2) < PC.minVel(2))]);
%Remove frames under min velocity
x(framesUnderMinVel,:) = [];
r(framesUnderMinVel,:) = [];
            
myX = PCAnalysis.binPosition(x,[20,1]);
myX = myX(:,1);

miCellS = cell([1000,1]);%zeros([10000,size(r,2)]);
%ComputeMI
parfor k = 1:size(miCellS,1)
    tempMI = zeros([1,size(r,2)]);
    for celli = 1:size(r,2)
        response = r(:,celli);
        if sum(response) == 0
            tempMI(celli) = nan;
        else
            tempMI(celli) = mutualinfo(response(randperm(length(response))),myX);
%             tempMI(celli) = miTemp(response(randperm(length(response)))',myX',20);
        end
    end
    miCellS{k} = tempMI;
end
miCellSCat = cat(1,miCellS{:});

miCell = zeros([1,size(r,2)]);
for celli = 1:size(r,2)
%     response = r(:,celli);        
%     miCell(celli) = myMI(response(randperm(length(response)))',myX',20);
    response = r(:,celli);
    if sum(response) == 0
        tempMI(celli) = nan;
    else
        miCell(celli) = mutualinfo(response,myX);
%         miCell(celli) = miTemp(response',myX',20);
    end
end
    
meanMiS = mean(miCellSCat);
steMeanMiS = std(miCellSCat)./sqrt(size(miCellSCat,1));
%%
pvalVec = nan([1,size(r,2)]);
for celli = 1:size(r,2)
    pvalVec(celli) = sum(miCell(celli)<miCellSCat(:,celli))/length(miCellSCat(:,celli));
end
significativeVec = pvalVec < .05;

%%
isPC = zeros([1,size(r,2)]);
figure;
for celli = 1:size(r,2)
    clf;
    histogram(miCellSCat(:,celli));hold on;
    if significativeVec(celli)
        barColor = 'g';
    else
        barColor = 'r';
    end
    plot([miCell(celli),miCell(celli)],ylim,barColor,'lineWidth',2);
    title(['Cell: ',num2str(celli),' pVal ',num2str(pvalVec(celli))])
    pause;
end

