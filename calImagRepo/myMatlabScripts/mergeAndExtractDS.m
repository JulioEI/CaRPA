folders = {'D:\Storage\Processing\Mouse2015'};
joinFileTreshold = 60;
for f = 1:length(folders)
    concatFilesFromTime(folders{f},{'dfof_?\d*.h5$'},joinFileTreshold)
    concatFilesFromTime(folders{f},{'downsampled_?\d*.h5$'},joinFileTreshold)
end