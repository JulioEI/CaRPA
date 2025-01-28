function renameVidsInOrder( folder )
    
    aviFiles = dir([folder,filesep,'*.avi']);
    logFiles = dir([folder,filesep,'*.xml']);
    extension = '.xml';
    if isempty(logFiles)
        logFiles = dir([folder,filesep,'*.txt']);
        extension = '.txt';
    end
    
    if length(logFiles) == length(aviFiles)
        for k = 1:length(logFiles)
            newName = replace(logFiles(k).name,['log',extension],'behavior.avi');
            source = [aviFiles(k).folder,filesep,aviFiles(k).name];
            destination = [aviFiles(k).folder,filesep,newName];
            movefile(source,destination)
        end
    else
        error('Different num of logfiles and avifiles')
    end

end

