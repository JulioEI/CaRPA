function cellStruct = getCellStruct(mouseFolders,cellStruct)

fileRegexp = {'TracesAndEvents?\d*.mat$'};

%Dont recompute animals that have been computed
if nargin < 2
    cellStruct = [];
else
    newMouseFolders = mouseFolders;
    mouses = cell([1,length(mouseFolders)]);
    for mouseK = 1:length(mouseFolders)
        folderNameParts = strsplit(mouseFolders{mouseK},filesep);
        regexpOut = regexp(folderNameParts{end}, 'ouse\D*(\d+)','tokens');
        mouse = regexpOut{1}{1};
        mouses{mouseK} = mouse;
    end
    for mouseK = 1:length(cellStruct)
        for mouseI = 1:length(mouses)
            if strcmp(mouses{mouseI},cellStruct(mouseK).mouse)
                newMouseFolders{mouseI} = [];
                break;
            end
        end
    end
    mouseFolders = newMouseFolders(~cellfun(@isempty,newMouseFolders));
end

%Compute the new amimals
for mouseK = 1:length(mouseFolders)
    disp([num2str(mouseK),'/',num2str(length(mouseFolders))])
    folderNameParts = strsplit(mouseFolders{mouseK},filesep);
    regexpOut = regexp(folderNameParts{end}, 'ouse\D*(\d+)','tokens');
    mouse = regexpOut{1}{1};
    newStruct(mouseK).mouse = mouse;      
    
    allItems = dir(mouseFolders{mouseK});
    allItems = allItems(3:end); %remove WINDOWS metadata folders
    allFolders = allItems([allItems.isdir]);
    dayNames = {allFolders.name};
    for dayK = 1:length(dayNames)
        newStruct(mouseK).days(dayK).name = dayNames{dayK};

        fullPath = fullfile(mouseFolders{mouseK}, dayNames{dayK});
        allInsideitems = dir(fullPath);
        allInsideitems = allInsideitems(3:end); %remove WINDOWS metadata folders
        allFiles = allInsideitems(~[allInsideitems.isdir]);
        fileIdx = find(fileSolution.checkRegexp({allFiles.name},fileRegexp));
        if isempty(fileIdx)
            warning(['cannot find files for ' dayNames{dayK}])
        else
            traceFiles = {allFiles(fileIdx).name};
            newStruct(mouseK).days(dayK).cells = zeros([1,length(traceFiles)]);
            newStruct(mouseK).days(dayK).frames = zeros([1,length(traceFiles)]);
            
            for fileK = 1:length(traceFiles)
                out = load([mouseFolders{mouseK},filesep,dayNames{dayK},filesep,traceFiles{fileK}]);
                framesX = size(out.tracesEvents.position,1);
                framesR = size(out.tracesEvents.rawTraces,1);
                if framesX~=framesR
                    warning(['DIFFERENT X AND R IN ',mouseFolders{mouseK},filesep,dayNames{dayK},filesep,traceFiles{fileK}]);
                    disp('x_________')
                    disp(framesX)
                    disp('r_________')
                    disp(framesR)
                    askFixDiff([mouseFolders{mouseK},filesep,dayNames{dayK},filesep,traceFiles{fileK}],out)
                end
                numCells = size(out.tracesEvents.rawTraces,2);
                newStruct(mouseK).days(dayK).cells(fileK) = numCells;
                newStruct(mouseK).days(dayK).frames(fileK) = framesX;
            end
        end
    end
end
cellStruct = [cellStruct,newStruct];
end

function askFixDiff(fullFile,out)
    answer = questdlg('Would you like to fix (cutting) the difference?', ...
        'IMPORTANT', ...
        'Yes please','No thank you','No thank you');
    % Handle response
    switch answer
        case 'Yes please'
            disp([' fixing...'])
            framesX = size(out.tracesEvents.position,1);
            framesR = size(out.tracesEvents.rawTraces,1);
            if framesX < framesR
                framesS = framesX;
                framesL = framesR;
            else
                framesS = framesR;
                framesL = framesX;                
            end
            for field = fields(out.tracesEvents)'
                if size(out.tracesEvents.(field{1}),1) == framesL
                    out.tracesEvents.(field{1}) = out.tracesEvents.(field{1})(1:framesS,:);
                end
            end
            tracesEvents = out.tracesEvents;
            save(fullFile,'tracesEvents')
        case 'No thank you'
            disp([' okay, skypping...'])
    end
end
