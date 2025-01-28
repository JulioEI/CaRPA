function plotTracesOverMovie(pos,vel,mov)
if isstr(mov)
    mov = VideoReader(mov);
end

x = pos(:,1)';
y = pos(:,2)';
z = zeros(size(x));
col = vel';
colSat = col;
colSat(col > quantile(col,.95)) = quantile(col,.95);
colSat(col < quantile(col,.05)) = quantile(col,.05);
figure;surface([x;x],[y;y],[z;z],[colSat;colSat],...
'facecol','no',...
'edgecol','interp',...
'linew',2);
colorbar
axis square;

figure;colormap gray;axis square;colorbar
cMap = hot(100);
cValue = round(99*(colSat - min(colSat))./(max(colSat)-min(colSat)))+1;
h = imagesc(mov.read(1));
hold on;
p = plot(pos(1,1),pos(1,2),'bx','markerSize',20);
for k = 1:mov.numberOfFrames
  h.CData = mov.read(k);caxis([0,1]);
  plot(pos(k,1),pos(k,2),'.','color',cMap(cValue(k),:));
  hold on;
  p.XData = pos(k,1);
  p.YData = pos(k,2);
  drawnow;%pause(0.1);
end
end

