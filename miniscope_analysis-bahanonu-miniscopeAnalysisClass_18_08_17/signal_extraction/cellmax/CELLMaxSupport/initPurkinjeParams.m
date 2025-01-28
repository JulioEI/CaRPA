function cInit = initPurkinjeParams(imgs, cellWidth)

nCells=ceil(size(imgs,2)/(cellWidth-1));
cInit=zeros(nCells, numel(imgs(:,:,1)));

gaussForConv=calcCellImgs([10 10 cellWidth cellWidth 0], [20 20]);

for cInd=1:nCells

    thisImage=0.01*ones(size(imgs(:,:,1)));
    oneInds=cInd*(cellWidth-1)+(1:cellWidth);
    oneInds(oneInds<1)=[];
    oneInds(oneInds>size(imgs,2))=[];
    thisImage(:,oneInds)=1;
    thisImage=conv2(thisImage, gaussForConv, 'same');
    thisImage=thisImage/sum(thisImage(:));

    cInit(cInd,:)=thisImage(:);

end

cInit(sum(isnan(cInit),2)>0,:)=[];
    
