sftp = mysftp;
destinationPath = 'C:\Users\csc\Desktop\caImagg\AllData';
tracesEventsRegexp = {'decisions_?\d*.mat$','Sorted_?\d*.mat$'};%{'TracesAndEvents?\d*.mat$'};%{'Sorted_?\d*.mat$','decisions_?\d*.mat$'};
archivePath = '/archive/pjercog/Processed_Data_CalcImaging';
mouseList = {'2028','2029','2023','2010','2011','2012','2026','2019','2021','2025','2024','2022'};
%%
targetFolders = cell([length(mouseList),1]);
for mouseK = 1:length(mouseList)
    sftp.cd(archivePath)
    
    animalFolders = sftp.commandToVar('ls');
    correctAnimalFolder = find(fileSolution.checkRegexp(animalFolders',{mouseList{mouseK}}));
    if isempty(correctAnimalFolder)
       error('Could not find the logs, animal missing')
    elseif length(correctAnimalFolder) > 1
%         s = listdlg('PromptString','Warning,duplicate folders, select the appropiate one',...
%                     'SelectionMode','single',...
%                     'ListString',animalFolders(correctAnimalFolder)');
%         correctAnimalFolder = correctAnimalFolder(s);
        correctAnimalFolder = correctAnimalFolder(2);
    end
    targetFolders{mouseK} = animalFolders{correctAnimalFolder};
end
%
for mouseK = 1:length(targetFolders)
    disp([num2str(mouseK),'/',num2str(length(targetFolders))])
    mousePath = [destinationPath,filesep,'Mouse',mouseList{mouseK}];
    mkdir(mousePath)
    sftp.cd(archivePath)
    sftp.cd(targetFolders{mouseK});
    dayFolders = sftp.commandToVar('ls')';
    currentPath = sftp.commandToVar('pwd');
    for folder = dayFolders
       if isempty(strfind(folder{1},'.')) %Is a folder, not a file
           sftp.cd(currentPath{1})
           sftp.cd(folder{1})
           filesInFolder = sftp.commandToVar('ls')';
           tracesEventsFiles = find(fileSolution.checkRegexp(filesInFolder,tracesEventsRegexp));
           if isempty(tracesEventsFiles)
               warning([folder{1},' has no tracesEvents files'])
           else
              serverFiles = cellfun(@(x) ['"',x,'"'],filesInFolder(tracesEventsFiles),'UniformOutput',0);
              folderPath = [mousePath,filesep,folder{1}];
              mkdir(folderPath)
              sftp.download(serverFiles,folderPath);
           end
       end
    end           
end
