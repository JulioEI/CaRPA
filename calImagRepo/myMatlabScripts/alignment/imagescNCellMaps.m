function imagescNCellMaps(allImages,myLegend,displayOverlap,myColorMap)
    if nargin < 4
        myColorMap = lines;%rand([size(lines,1),3]);
        if nargin < 3
            displayOverlap = 1;
            if nargin < 2
                myLegend = [];
            end
        end
    end
    
    backgroundImg = zeros([size(allImages{1},1),size(allImages{1},2),3]);
    
    allImages = cellfun(@(x) (x-min(x(:)))./(max(x(:))-min(x(:))),allImages,'UniformOutput',0);
    for k = 1:length(allImages)
        imageT = allImages{k};%max(thresholdImages(allImages{k}),[],3)';  
        for i = 1:3
            backgroundImg(:,:,i) = max(cat(3,backgroundImg(:,:,i),imageT*myColorMap(k,i)),[],3);
        end
    end
    
    %display overlaps
    if displayOverlap
        for i = 1:3
            backgroundImg(:,:,i) = max(cat(3,backgroundImg(:,:,i),(1.2*prod(cat(3,allImages{:}),3))),[],3); 
        end
    end
    
    imagesc(zeros(size(backgroundImg)))
    hold on;
    h = imagesc(backgroundImg);
    % plot(commonCenters(:,:,1),commonCenters(:,:,2),'.','markersize',10);  
    % for k = 1:length(commonIn1And2)
    %     plot(squeeze(commonCenters(k,:,1)),squeeze(commonCenters(k,:,2)),'w:');  
    % end
    alphaMask = 0.3*(0~=(sum(backgroundImg,3)));
    set(h,'AlphaData',alphaMask)

    if ~isempty(myLegend)
        p = [];
        for k = 1:length(allImages) 
            p(k) = plot([nan],[nan],'o','color',myColorMap(k,:));
        end
        legend(p,myLegend)
    end
    
end

