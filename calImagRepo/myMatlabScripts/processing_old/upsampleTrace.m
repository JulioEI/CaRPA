mousePath = 'E:\Processing\Mouse2028 - 500itEM\20150303-icx\';
movie20hz = loadMovieList([mousePath,'2015_03_03_p000_mouse2028_NULL000_turboreg_crop_dfof_1.h5']);

% hinf = hdf5info([mousePath,'2015_03_03_p000_mouse2028_NULL000_turboreg_crop_dfof_1.h5']);
% movie20hz = hdf5read(hinf.GroupHierarchy.Datasets);
robustMovie = movie20hz;
robustMovie(isnan(movie20hz))=min(movie20hz(:));
[upScaledPhi, upFilteredTraces] = recalcPhiAndDetectEvents(robustMovie,emAnalysisOutput.cellImages,emAnalysisOutput.dsCellTraces,emAnalysisOutput.CELLMaxoptions,'runEventDetection',0);
sum(isnan(upFilteredTraces(:)))/numel(upFilteredTraces)*100
k = 100;
figure;plot(upFilteredTraces(k,:));hold on;plot(linspace(1,length(upFilteredTraces(k,:)),length(emAnalysisOutput.cellTraces(k,:))),emAnalysisOutput.cellTraces(k,:));legend('20k','5k')
figure;plot(upFilteredTraces(k,:));hold on;plot(linspace(1,length(upFilteredTraces(k,:)),length(emAnalysisOutput.cellTraces(k,:))),emAnalysisOutput.dsCellTraces(k,:));legend('20k','5k')
