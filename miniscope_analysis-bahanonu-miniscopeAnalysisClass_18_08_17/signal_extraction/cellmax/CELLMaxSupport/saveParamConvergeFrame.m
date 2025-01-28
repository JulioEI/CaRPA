function saveParamConvergeFrame(theseCentroids, trueCentroids, doppCentroids, writerObj, xlims, ylims)

h=figure(27);
if ~isempty(trueCentroids)
    figure(h); hold off; plot(trueCentroids(:,1), trueCentroids(:,2), 'r.', 'Markersize', 14); hold on;
end
plot(theseCentroids(:,1), theseCentroids(:,2), 'k.', 'Markersize', 14); drawnow
if ~isempty(doppCentroids)
    figure(h); plot(doppCentroids(:,1), doppCentroids(:,2), 'bo', 'Markersize', 10); drawnow
end
xlim(xlims); ylim(ylims); drawnow
thisFrame=getframe(h);
writeVideo(writerObj,thisFrame);