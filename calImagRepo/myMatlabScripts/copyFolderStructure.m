source = 'E:\Processing\Mouse2021';
target = 'C:\Users\csc\Desktop\caImagg\AllData\Mouse2021';
tracesEventsRegexp = {'TracesAndEvents?\d*.mat$'};

folders = dir(source);
folders = folders([folders.isdir]);
folders = folders(3:end);
for k = 1:length(folders)
    
    
    files = dir([source,filesep,folders(k).name,filesep,'*.mat']);
    files = files(find(fileSolution.checkRegexp({files.name},tracesEventsRegexp)));
    
    if ~isempty(files)
        mkdir([target,filesep,folders(k).name])
        for file = {files.name}
            copyfile([source,filesep,folders(k).name,filesep,file{1}],[target,filesep,folders(k).name,filesep,file{1}])
        end
    else
        warning(['no files for ',folders(k).name])
    end
    
end