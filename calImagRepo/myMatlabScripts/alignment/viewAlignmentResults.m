% load('E:\Processing\newTurboreg\Mouse2028 - reg - linear\globalIDs.mat')%load('E:\Processing\Mouse2028 - 200itEM\globalIDs.mat')

load('E:\Processing\temp\extractionDecisionsm2028\day20150303 m2028\2015_03_03_p000_mouse2028_day20150303_emAnalysis.mat')
load('E:\Processing\temp\extractionDecisionsm2028\day20150303 m2028\2015_03_03_p000_mouse2028_day20150303_emAnalysisSorted.mat')
images1 = emAnalysisOutput.cellImages;
images1= images1(:,:,logical(validCellMax));

load('E:\Processing\temp\extractionDecisionsm2028\day20150305 m2028\2015_03_05_p000_mouse2028_day20150305_emAnalysis.mat')
load('E:\Processing\temp\extractionDecisionsm2028\day20150305 m2028\2015_03_05_p000_mouse2028_day20150305_emAnalysisSorted.mat')
images2 = emAnalysisOutput.cellImages;
images2= images2(:,:,logical(validCellMax));

%%
% 
% baseImg1 = mean(images1,3)*20;
% baseImg2 = mean(images2,3)*20;
% figure;
% for k = 1:length(idxFirstDay)
%     try
%         i1 = images1(:,:,globalIDs(idxFirstDay(k),1));
%         i2 = images2(:,:,globalIDs(idxFirstDay(k),2));
%         [cy1,cx1] = find(max((max(i1)))==i1);
%         [cy2,cx2] = find(max((max(i2)))==i2);
%         disp(sqrt((cx1-cx2)^2+(cy1-cy2)^2));
%         
%         subplot(2,2,[1,2]);imagesc(mean(cat(3,baseImg1,baseImg2),3)+images1(:,:,globalIDs(idxFirstDay(k),1))+images2(:,:,globalIDs(idxFirstDay(k),2)));axis square;hold on;plot(cx1,cy1,'r.');plot(cx2,cy2,'g.')
%         subplot(2,2,3);imagesc(baseImg1+images1(:,:,globalIDs(idxFirstDay(k),1)));axis square;hold on;plot(cx1,cy1,'r.');plot(cx2,cy2,'g.')
%         subplot(2,2,4);imagesc(baseImg2+images2(:,:,globalIDs(idxFirstDay(k),2)));axis square;hold on;plot(cx1,cy1,'r.');plot(cx2,cy2,'g.')
% 
%     catch
%         clf;
%     end
%         pause;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%
cellDist = [];
for k = 1:length(idxFirstDay)
    try
        i1 = images1(:,:,globalIDs(idxFirstDay(k),1));
        i2 = images2(:,:,globalIDs(idxFirstDay(k),2));
        [cx1,cy1] = find(max((max(i1)))==i1);
        [cx2,cy2] = find(max((max(i2)))==i2);
        cellDist = [cellDist,sqrt((cx1-cx2)^2+(cy1-cy2)^2)];
    end
end
figure;histogram(cellDist,15)
%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Remove duplicates
for k = 1:size(globalIDs,2)
    sessionK = globalIDs(:,k);
    sessionK = sessionK(sessionK~=0);
    repeatedGlobalIdx = find(hist(sessionK,unique(sessionK))>1);
    
    for repeatedIdx = repeatedGlobalIdx
        disp([num2str(k),' session, index ',num2str(repeatedIdx)])
%         globalIDs(find(globalIDs(:,k)==repeatedIdx,1,'last'),k) = 0;
    end
end
%%
idxFirstDay = find(globalIDs(:,1));

cellLog = zeros([1,size(globalIDs,2)]);
for k = 1:size(globalIDs,2)
    cellLog(k) = length(find(prod(globalIDs(:,1:k),2)));
end

cellLog2 = zeros([1,size(globalIDs,2)]);
for k = 1:size(globalIDs,2)
    cellLog2(k) = sum(sum(globalIDs(:,setdiff(1:size(globalIDs,2),k)),2) & globalIDs(:,k));
end

cellsInDay = zeros([1,size(globalIDs,2)]);
for k = 1:size(globalIDs,2)
    cellsInDay(k) = length(find(globalIDs(:,k)));
end

cellLog3 = zeros([1,size(globalIDs,2)]);
for k = 1:size(globalIDs,2)
    cellsInOtherSessions = sum(globalIDs(:,setdiff(1:size(globalIDs,2),k))~=0,2);
    cellsInOtherSessions = cellsInOtherSessions(cellsInOtherSessions&globalIDs(:,k));
    cellLog3(k) = sum(cellsInOtherSessions > floor(size(globalIDs,2)/2));
end

cellMat = zeros(size(globalIDs,2));
for j = 1:size(globalIDs,2)
    for k = 1:size(globalIDs,2)
        cellsInOtherSessions = sum(globalIDs(:,setdiff(1:size(globalIDs,2),k))~=0,2);
        cellsInOtherSessions = cellsInOtherSessions(cellsInOtherSessions&globalIDs(:,k));
        cellMat(j,k) = 100*sum(cellsInOtherSessions > floor(size(globalIDs,2)-j))/cellsInDay(k);
    end
end

figure;imagesc(flipud(cellMat));axis square;xlabel('days');ylabel('shared in at least');colorbar;set(gca, 'YTick', 1:size(globalIDs,2));set(gca, 'XTick', 1:size(globalIDs,2));set(gca,'YDir','normal')


figure;bar([100*cellLog./cellsInDay;100*cellLog2./cellsInDay;100*cellLog3./cellsInDay]')
ylabel('%Cells of day')
xlabel('Day')
legend('%Cells aligned to all previous days','%Cells aligned to any other day','%Cells in half the other sessions')


%%
commonIn1And2 = find(prod(globalIDs(:,[1,2]),2));
commonCenters = zeros([length(commonIn1And2),2,2]);
for k = 1:length(commonIn1And2)
        i1 = images1(:,:,globalIDs(commonIn1And2(k),1));
        i2 = images2(:,:,globalIDs(commonIn1And2(k),2));
        [cx1,cy1] = find(max((max(i1)))==i1);
        [cx2,cy2] = find(max((max(i2)))==i2);
        commonCenters(k,1,:) = [cx1,cy1];
        commonCenters(k,2,:) = [cx2,cy2];
end

allImages = {images1,images2};
myColorMap = lines;

backgroundImg = zeros([size(allImages{1},1),size(allImages{1},2),3]);
for k = 1:length(allImages)
    imageT = max(thresholdImages(allImages{k}),[],3)';  
    imageT = (imageT-min(imageT(:)))/(max(imageT(:))-min(imageT(:)));
    
    for i = 1:3
        backgroundImg(:,:,i) = max(cat(3,backgroundImg(:,:,i),imageT*myColorMap(k,i)),[],3);
    end
end

figure;
imagesc(zeros(size(backgroundImg)))
hold on;
h = imagesc(backgroundImg);
% plot(commonCenters(:,:,1),commonCenters(:,:,2),'.','markersize',10);  
% for k = 1:length(commonIn1And2)
%     plot(squeeze(commonCenters(k,:,1)),squeeze(commonCenters(k,:,2)),'w:');  
% end
alphaMask = 0.3*(0~=(sum(backgroundImg,3)));
set(h,'AlphaData',alphaMask)