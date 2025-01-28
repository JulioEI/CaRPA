%This uses the code of batch_processing_custom to process multiple
%concat_recordings in batch. Search #CUSTOM_CODE! to see changes.

%This code assuemes that each concat file contains only one dataset.


%Grab concat files:

%Right now 1 folder contains 1 day and 1 concat.

days = {'03'};%{'11','09','07','05','03','011','012','281','282'};%'31','29','27','25','23','21','19','17','15','13'};%;
myDir = 'E:\Processing\Mouse2028 - 500itEM\';
pathToMiji = 'C:\Users\csc\Desktop\Fiji.app\scripts';

allFiles = dir(myDir);

concatFiles = {};
concatFolders = {};

for day = days
    
    regList = regexp({allFiles(:).name},['^201\d\d\d',day{1},'-icx$']);

    
    indexAll = find(~cellfun(@isempty,regList));
    
    if length(indexAll) == 1
        insideFiles = dir([myDir,allFiles(indexAll).name]);
        regList = regexp({insideFiles(:).name},'concat');
        indexInside = find(~cellfun(@isempty,regList));  
        
        if length(indexInside) == 1
        
            concatFiles = [concatFiles,[myDir,allFiles(indexAll).name,'\',insideFiles(indexInside).name]];
            concatFolders = [concatFolders,[myDir,allFiles(indexAll).name]];
        
        elseif length(indexInside) > 1
            error('Regexp returned two files')
            
        else
            warning(['Day ',day{1},' ignored, no files found'])
        end
            
    elseif length(indexAll) > 1
        error('Regexp returned two files')
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
        obj = behaviorAnalysis;
        obj.modelAddNewFolders(concatFolders{k});
        warning('off')
        fileInfo = hdf5info(concatFiles{k});
        datasetName = fileInfo.GroupHierarchy.Datasets.Name;
        obj.modelPreprocessMovie(datasetName,pathToMiji,pLog(k,:))
        disp(i)
        warning('on')
    end
end