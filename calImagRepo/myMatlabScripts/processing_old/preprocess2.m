pathToMiji = 'C:\Users\csc\Desktop\Fiji.app\scripts';

folders = {'E:\Processing\newTurboreg\Mouse2028'};%{'E:\Processing\Mouse2029','E:\Processing\Mouse2028 - 125itEM','E:\Processing\Mouse2028 - 200itEM','E:\Processing\Mouse2028 - 500itEM','E:\Processing\Mouse2028'};
daysExcluded = {};%{'201502281-icx','201502282-icx','201503011-icx','201503012-icx'};%{'201502281-icx','201502282-icx','201503011-icx','201503012-icx'};
daysIncluded = {'20150305-icx'};

concatFiles = {};
concatFolders = {};
for k = 1:length(folders)
    allFiles = dir(folders{k});
    timePerDays = [];
    predNPerDays = [];
    for dayFolder = {allFiles.name}
        if ~isempty(strfind(dayFolder{1},'-icx')) &&  1==sum(strcmp(dayFolder{1},daysIncluded)) &&  1~=sum(strcmp(dayFolder{1},daysExcluded))%Avoids only specified days
            insideFiles = dir([folders{k},'\',dayFolder{1}]);
            for file = {insideFiles.name}
                if regexp(file{1},'concat')
                    concatFiles = [concatFiles,[folders{k},'\',dayFolder{1},'\',file{1}]];
                    concatFolders = [concatFolders,[folders{k},'\',dayFolder{1}]];
                end
            end
        end
    end
end

%Setup all turboreg areas.
figure;
pLog = zeros(length(concatFiles),4);
i = 0;
for concat = concatFiles
    i = i + 1;
    hinfo = hdf5info(concat{1});
    try
        hReadInfo = hinfo.GroupHierarchy.Datasets(1);
    catch
        hReadInfo = hinfo.GroupHierarchy.Groups.Datasets(1);
    end
    xDim = hReadInfo.Dims(1);
    yDim = hReadInfo.Dims(2);
    % select the first frame from the dataset
    thisFrame = readHDF5Subset(concat{1},[0 0 100],[xDim yDim 1],'datasetName',hReadInfo.Name);  
    
   
    imagesc(thisFrame); axis image; colormap gray; title('Click to drag-n-draw region. Double-click region to continue.');xlabel([num2str(i),' of ',num2str(length(concatFiles))])
    p = round(wait(imrect));
    pLog(i,:) = p;
end

%Start automation
%%
for k = 1:length(concatFolders)
    for i = 10
        disp([num2str(k),' of ',num2str(length(concatFolders))])
        clear obj
        obj = miniscopeAnalysis;
        obj.modelAddNewFolders(concatFolders{k});
        warning('off')
        fileInfo = hdf5info(concatFiles{k});
        datasetName = fileInfo.GroupHierarchy.Datasets.Name;
        obj.modelPreprocessMovie(datasetName,pathToMiji,pLog(k,:))
        disp(i)
        warning('on')
    end
end