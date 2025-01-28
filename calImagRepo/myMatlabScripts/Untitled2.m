vAvi = VideoWriter('Mouse-3005-21080510_112915-icx-dfof.avi');
vMp4 = VideoWriter('Mouse-3005-21080510_112915-icx-dfof.mp4');
vAvi.FrameRate = 20;
vAvi.Quality = 100;
vMp4.FrameRate = 20;
vMp4.Quality = 100;

calcium = 'D:\Storage\Processing\Mouse3005\Mouse-3005-21080510-icx\Mouse-3005-21080510_112915-icx-dfof.h5';
hinf = hdf5info(calcium);
dataset = hdf5read(hinf.GroupHierarchy.Datasets);

open(vAvi);


upLim = quantile(dataset(:),0.995);
dwLim = quantile(dataset(:),0.005);
for k = 1:size(dataset,3)
   imagesc(dataset(:,:,k),[dwLim,upLim]);
   set(gca,'Visible','off')
   frame = getframe(gcf);
   writeVideo(vAvi,frame);
   writeVideo(vMp4,frame);
end

close(vAvi);
close(vMp4);