function fixTracesTEMP(mouseFolders)

fileRegexp = {'TracesAndEvents?\d*.mat$'};

for mouseK = 1:length(mouseFolders)
    disp([num2str(mouseK),'/',num2str(length(mouseFolders))])
    
    allItems = dir(mouseFolders{mouseK});
    allItems = allItems(3:end); %remove WINDOWS metadata folders
    allFolders = allItems([allItems.isdir]);
    dayNames = {allFolders.name};
    for dayK = 1:length(dayNames)
        disp(dayNames{dayK})
        fullPath = fullfile(mouseFolders{mouseK}, dayNames{dayK});
        allInsideitems = dir(fullPath);
        allInsideitems = allInsideitems(3:end); %remove WINDOWS metadata folders
        allFiles = allInsideitems(~[allInsideitems.isdir]);
        fileIdx = find(fileSolution.checkRegexp({allFiles.name},fileRegexp));
        if isempty(fileIdx)
            warning(['cannot find files for ' dayNames{dayK}])
        else
            traceFiles = {allFiles(fileIdx).name};            
            for fileK = 1:length(traceFiles)
                fullFile = [mouseFolders{mouseK},filesep,dayNames{dayK},filesep,traceFiles{fileK}];
                out = load(fullFile);
                framesX = size(out.tracesEvents.position,1);
                framesR = size(out.tracesEvents.rawTraces,1);
                if framesX~=framesR
                    error(['DIFFERENT X AND R IN ',mouseFolders{mouseK},filesep,dayNames{dayK},filesep,traceFiles{fileK}]);
                end
                x = out.tracesEvents.position;
                outliers = x>(mean(x)+3*std(x));
                if sum(outliers(:,1)) ~= 0 %only for x
                    warning('OUTLIERS DETECTED, MANUAL ACTION REQUIERED')
                    askFixOutliers(fullFile,out)
                    out = load(fullFile);
                end
                if size(out.tracesEvents.cellAnatomicLocat,1) ~= size(out.tracesEvents.rawTraces,2)
                    warning('ERROR IN CELL LOCAT DETECTED, FIXING...')
                    decFile = strrep(fullFile,'TracesAndEvents','emAnalysisSorted');
                    dec = load(decFile);
                    out.tracesEvents.cellAnatomicLocat = out.tracesEvents.cellAnatomicLocat(logical(dec.validCellMax),:);
                    clear tracesEvents
                    tracesEvents = out.tracesEvents;
                    save(fullFile,'tracesEvents')
                end
            end
        end
    end
end
end
function askFixOutliers(fullFile,out)
%     answer = questdlg('Would you like to fix (cutting) the outliers?', ...
%         'IMPORTANT', ...
%         'Yes please','No thank you','No thank you');
    answer = 'Yes please';
    % Handle response
    switch answer
        case 'Yes please'
            disp([' fixing...'])
            x = out.tracesEvents.position;
            answer = [];
            isOk = [];
            while isempty(isOk)
                while isempty(answer)
                    figure;plot(x);
                    prompt = {'Enter lower cut:','Enter higher cut:'};
                    title = 'Input';
                    definput = {'1',num2str(length(x))};
                    opts.WindowStyle = 'normal';
                    answer = inputdlg(prompt,title,[1 40; 1 40],definput,opts);
                    close all;
                end
                lB = str2num(answer{1});
                uB = str2num(answer{2});
                figure;plot(x(lB:uB,:));
                prompt = {'Is this ok?'};
                title = 'Input';
                definput = {'Press OK to confirm, press RETURN to try again.'};
                opts.WindowStyle = 'normal';
                isOk = inputdlg(prompt,title,[1 40],definput,opts);
                answer = [];
                close all;
            end
            
            framesX = size(out.tracesEvents.position,1);
            for field = fields(out.tracesEvents)'
                if size(out.tracesEvents.(field{1}),1) == framesX
                    out.tracesEvents.(field{1}) = out.tracesEvents.(field{1})(lB:uB,:);
                end
            end
            tracesEvents = out.tracesEvents;
            save(fullFile,'tracesEvents')
        case 'No thank you'
            disp([' okay, skypping...'])
    end
end
