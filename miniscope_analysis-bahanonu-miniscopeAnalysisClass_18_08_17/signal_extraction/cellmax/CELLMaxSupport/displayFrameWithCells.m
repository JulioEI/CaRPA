function hAx=displayFrameWithCells(h, dsImgs,currentFrame,cellROIs,currentCellInd,xLims,yLims,clims)

figure(h)
hold off
imagesc(dsImgs(:,:,currentFrame))
set(gca, 'Fontsize', 14)
title('l fwd | k back | o fwd 20 | i back 20 | p play 20 | d delete | q quit')
xlabel(['Frame ' num2str(currentFrame)])
xlim(xLims)
ylim(yLims)
hAx=gca;
set(gca, 'CLim', clims)
hold on
for cInd=1:currentCellInd
    thisBorder=cellROIs{cInd}.border;
    if cInd==currentCellInd
        plot(thisBorder(:,1), thisBorder(:,2), 'w')
    else
        plot(thisBorder(:,1), thisBorder(:,2), 'k')
    end
end

