function folderStruct = getFolderStruct( rootFolder,fileRegexp)
if nargin < 2
    fileRegexp = {'downsampled_?\d*.h5$'};
end

dateParser = 'yyyymmdd';
timeParser = 'HHMMSS';

allItems = dir(rootFolder);
allItems = allItems(3:end); %remove WINDOWS metadata folders
allFolders = allItems([allItems.isdir]);
folderStruct = [];
i = 1;
for dayFolder = {allFolders.name}
    regexpOut = regexp(dayFolder{1}, '(\d+)','tokens');
    day = regexpOut{end}{1};%Assumess date is the last sequence of numbers in the folder name
    fullPath = fullfile(rootFolder, dayFolder{1});    

	folderStruct(i).path = fullPath;
	folderStruct(i).folderName = dayFolder{1};
    folderStruct(i).date = datevec(day,dateParser);
%     folderStruct(i).mouse = obj.mouse;

    allInsideitems = dir(fullPath);
    allInsideitems = allInsideitems(3:end); %remove WINDOWS metadata folders
    allFiles = allInsideitems(~[allInsideitems.isdir]);
    concatIdx = find(fileSolution.checkRegexp({allFiles.name},fileRegexp));
    if isempty(concatIdx)
        warning(['cannot find concat files for ' dayFolder{1}])
    else
        concatFiles = {allFiles(concatIdx).name};
        k = 1;
        for concatFile = concatFiles
            noExtensionName = split(concatFile{1},'.');
            regexpOut = regexp(noExtensionName{1}, '(\d+)','tokens');
            session = regexpOut{end}{1};
            folderStruct(i).sessions(k).date = datevec([day,session],[dateParser,timeParser]);
            folderStruct(i).sessions(k).fileName = concatFile{1};
            k = k + 1;
        end
        %Sort sessions in descending time order
        [~,sortedIdx] = sort(cellfun(@datenum,{folderStruct(i).sessions.date}));
        folderStruct(i).sessions = folderStruct(i).sessions(sortedIdx);
        i = i+1; %If no concat files were found, ignore the folder
    end
end
end

