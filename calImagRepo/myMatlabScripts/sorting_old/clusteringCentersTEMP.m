meanValVect = zeros([1,size(signalImagesCrop,3)]);
    for j = 1:size(signalImagesCrop,3)
        tempEvent = signalImagesCrop(:,:,j);
        meanValVect(j) = nanmean(tempEvent(cellMask));
    end

figure;
maxX = zeros([1,size(signalImagesCropWithCell,3)]);
maxY = maxX;
for j = 1:size(signalImagesCropWithCell,3)
[maxX(j),maxY(j)] = find(signalImagesCropWithCell(:,:,j) == max(max(signalImagesCropWithCell(:,:,j))));
end
idx = fkmeans([maxX',maxY'],3,mahoptions);
myColorMap = lines;
for k = 1:length(unique(idx))
scatter(maxX(idx==k),maxY(idx==k),'markerfacecolor',[1,1,1],'markeredgecolor',myColorMap(k,:));hold on
text(maxX(idx==k)+0.1,maxY(idx==k)+1*rand(size(maxX(idx==k))),num2cell(find(idx==k)+1),'color',myColorMap(k,:))
end

