classdef carpa < handle
% insert definition here
    properties
        rootFolder;
        mouse;
        spatialDS = 4;
        joinFileTreshold = 120; %Will join the concat files if they are closer than treshold (minutes). To not join put 0.
        analysisType = 'cellmax'; %'pcaica';
        archiveRawCaPath = '/archive/pjercog/Raw_Data_CalcImaging';
        showProcessed = 0;
    end
    
    properties(SetAccess = private)
        sftp;
        folderStruct;
        askForFiles = 1;
    end
    
    properties (Constant)
        %REGEXP STUFF FOR FINDING FILES       
        analysisFileRegexp = {'^Mouse.*emAnalysis_?\d*.mat$','^Mouse.*pcaicaAnalysis_?\d*.mat$'};
        downsampleRegexp = {'^Mouse.*downsampled_?\d*.h5$'};
        dfofRegexp = {'^Mouse.*dfof_?\d*.h5$'};
%         concatRegexp = {'concat_recording.*\d+_\d+.hdf5$','concat_recording.*\d+_\d+.h5$'};%{'concat_recording.*.h5$','concat_recording.*.hdf5$'};
        concatRegexp = {'recording.*\d+_\d+.hdf5$','recording.*\d+_\d+.h5$'};%{'concat_recording.*.h5$','concat_recording.*.hdf5$'};

        decisionsRegexp = {'^Mouse.*Sorted_?\d*.mat$','^Mouse.*decisions_?\d*.mat$'}

        behaviorRegexp = {'^Mouse.*.avi$'}
        logRegexp = {'^Mouse.*.xml$','^Mouse.*.txt$'}
%         {'^Mouse.*.txt$'}
        tracesEventsRegexp = {'^Mouse.*TracesAndEvents?\d*.mat$'}
        
        %Time parsers
        dateParser = 'yyyymmdd'; %Reads and prints time with this format
        timeParser = 'HHMMSS';
    end   
    
    methods
        
        function obj = carpa
           answer = questdlg('What would you like to do', ...
           'Welcome to CaRPA!', ...
           'Work on existing animal','Create file structure for new animal','Download files from server','Work on existing animal');
            switch answer
                case 'Work on existing animal'
                    obj.rootFolder = uigetdir('','Select the animal root folder');
                                      
                    regexpOut = regexp(obj.rootFolder, 'ouse\D*(\d+)','tokens');
                    try mouseCode = regexpOut{1}{1};catch;mouseCode = '0000';end
                    mouseCode = inputdlg('Enter the mouse code:','Mouse code',1,{mouseCode});
                    obj.mouse = mouseCode{1};
                    
                case 'Create file structure for new animal'
                    obj.newStructure;
                    disp('... structure sucessfully created');
                case 'Download files from server'
                    %INTITATE SERVER OBJECT
                    try
                        obj.sftp = mysftp('neurocomp.fcrb.es','4022');
                    catch ME
                        disp(ME.message)
                        disp('server not initiated')
                    end
                    obj.rootFolder = uigetdir('','Select the folder to put the downloaded files');
                    
                    regexpOut = regexp(obj.rootFolder, 'ouse\D*(\d+)','tokens');
                    try mouseCode = regexpOut{1}{1};catch;mouseCode = '0000';end
                    mouseCode = inputdlg('Enter the mouse code:','Mouse code',1,{mouseCode});
                    obj.mouse = mouseCode{1};
            
                    obj.downloadH5fromServer;
                    
                    obj.setExperimentNames;
                    
                otherwise
                    return;
            end
            
            obj.buildFolderStructure;
            
        end
        
        function newStructure(obj)
            obj.rootFolder = uigetdir('','Select the root folder of the recordings: '); 
            %Try to get the mouse name
            nameSep = strsplit(obj.rootFolder,filesep);nameSep = nameSep(~cellfun(@isempty,nameSep));
            folderName = nameSep{end};
            regexpOut = regexp(folderName, 'ouse\D*(\d+)','tokens');
            try mouseCode = regexpOut{1}{1};catch;mouseCode = '0000';end
            
            mouseCode = inputdlg('Enter the mouse code:','Mouse code',1,{mouseCode});
            obj.mouse = mouseCode{1};
            allFiles = dir(obj.rootFolder);allFiles=allFiles(~ismember({allFiles.name},{'.','..'}));
            
            if any([allFiles.isdir])
                error('To create a file structure only calcium files should be in the folder.')
            end
            
           filesPerExperiment = struct;
           k = 1;
           answer = questdlg('Does the mouse perform more than one experiment?', ...
           'Experiment data', ...
           'Yes','No','No');
            switch answer
                case 'No'
                    experimentName = inputdlg('Enter experiment name:','Experiment name',1);if isempty(experimentName);return;end
                    filesPerExperiment(k).name = experimentName{1};
                    filesPerExperiment(k).files = {allFiles.name};
                case 'Yes'
                    experimentName = inputdlg('Enter a first experiment name:','Experiment name',1);if isempty(experimentName);return;end
                    correctFiles = {allFiles.name};
                    while true
                       filesPerExperiment(k).name = experimentName{1};
                       filter = [];
                       for file = correctFiles
                           [~,rawName] = fileparts(file{1});
                           filter = [filter,rawName,'.*;'];
                       end
                       files = uigetfile(filter,'Select all files for this experiment',obj.rootFolder,'MultiSelect', 'on');if ~iscell(files);files = {files};end
                       filesPerExperiment(k).files = files;
                       correctFiles = setdiff(correctFiles,filesPerExperiment(k).files);
                       if isempty(correctFiles)
                           break;
                       end
                       experimentName = inputdlg('Enter another experiment name:','Experiment name',1);
                       k = k + 1;
                    end
                otherwise
                    return;
            end
                        
            for k = 1:length(filesPerExperiment)
                for file = filesPerExperiment(k).files
                    noExtensionName = split(file{1},'.');
                    regexpOut = regexp(noExtensionName{1}, '(\d+)','tokens');
                    day = regexpOut{1}{1};
                    dayFolderName = fullfile(obj.rootFolder,['Mouse-',obj.mouse,'-',day,'-',filesPerExperiment(k).name]);
                    warning('off', 'MATLAB:MKDIR:DirectoryExists');
                    mkdir(dayFolderName)
                    movefile(fullfile(obj.rootFolder,file{1}),fullfile(dayFolderName,file{1}))
                end
            end
        end
                
        function menu(obj)
            set(0, 'DefaultUICOntrolFontSize', 14)
            scnsize = get(0,'ScreenSize');dlgSize = [scnsize(3)*0.2 scnsize(4)*0.7];
            methods = {'PREPROCESS','EXTRACT','FILTER','POSTPROCESS','QUIT'};
            %%%% I removed the last element because display is for some
            %%%% reason changed that I can't figure out why !!!
            s = listdlg('PromptString','Select a method:','SelectionMode','multiple','ListString',methods,'ListSize',dlgSize);
            
            if any(s==length(methods))
                return;
            end
            
            if length(s) == 1
                choice = questdlg('Choose specific days/sessions?', ...
                    'IMPORTANT', ...
                    'Yes','No, process everything','Cancel','No, process everything');
            else
                choice = 'No, process everything';
            end
            
            switch choice
                case 'Yes'
                    obj.askForFiles = 1;
                case 'No, process everything'
                    obj.askForFiles = 0;
                case 'Cancel'
                    return;
            end
            
            if any(s==1)
                obj.preprocess;
            end
            if any(s==2)
                obj.extract;
            end            
            if any(s==3)
                obj.filter;
            end
            if any(s==4)
                obj.postprocess;
            end  
            
        end
        
        function selectedFiles = selectFiles(obj,type,mask)
            showProcessed = obj.showProcessed;
%             showProcessed = 1;
            if nargin < 3
                showProcessed = 1;
            end
            
            dayHasAFileOfType = ~cellfun(@isempty,{obj.folderStruct.(type)});
            files = {obj.folderStruct(dayHasAFileOfType).(type)};
            allFiles = cat(2,files{:});       
            
            if isempty(allFiles)
                selectedFiles = {};
                warning(['No ',type,' files exist'])
                return
            end
            
            if ~showProcessed
                files = carpaUtilities.maskSelectedFiles({obj.folderStruct(dayHasAFileOfType).(type)},{obj.folderStruct(dayHasAFileOfType).(mask)});
                allFiles = cat(2,files{:});                               
                if isempty(allFiles)
                    warning('All files processed, set showProcessed to 1 to reprocess files (Line 190:193)')
                    selectedFiles = {};
                    return
                end
            end
            
            allDays = unique({allFiles.folder});
            
            set(0, 'DefaultUICOntrolFontSize', 14)
            scnsize = get(0,'ScreenSize');dlgSize = [scnsize(3)*0.5 scnsize(4)*0.7];
            
            selectType = 'days';
            while true
                switch selectType
                    case 'days'
                        selection = listdlg('PromptString','Select a method:','SelectionMode','multiple','ListString',[{'SELECT BY INDIVIDUAL SESSIONS INSTEAD',''},allDays],'ListSize',dlgSize);
                        if selection == 1
                            selectType = 'sessions';
                        elseif selection == 2
                            error('very funny')
                        else
                            selectedFiles = files(selection-2);
                            break;
                        end
                    
                    case 'sessions'
                        selection = listdlg('PromptString','Select a method:','SelectionMode','multiple','ListString',[{'SELECT BY WHOLE DAYS INSTEAD',''},{allFiles.dispName}],'ListSize',dlgSize);
                        if selection == 1
                            selectType = 'days';
                        elseif selection == 2
                            error('very funny')
                        else
                            %group by days
                            selectedFiles = {};
                            sFiles = allFiles(selection-2);
                            folderList = {};
                            for k = 1:length(sFiles)
                                if ~any(strcmp(folderList,sFiles(k).folder))
                                    folderList = [folderList, sFiles(k).folder];
                                    selectedFiles = [selectedFiles,{sFiles(k)}];
                                else
                                    selectedFiles{end} = [selectedFiles{end},sFiles(k)];
                                end
                            end
                            break;
                        end                      
                end
            end 
            
        end
        
        function selectedFiles = getAllFilesWithType(obj,type,mask)
            showProcessed = obj.showProcessed;
            if nargin < 3
                showProcessed = 1;
            end
            
            dayHasAFileOfType = ~cellfun(@isempty,{obj.folderStruct.(type)});
            selectedFiles = {obj.folderStruct(dayHasAFileOfType).(type)};
           
            if isempty(selectedFiles)
                warning(['No ',type,' files exist'])
                return
            end
            
            if ~showProcessed
                selectedFiles = carpaUtilities.maskSelectedFiles({obj.folderStruct(dayHasAFileOfType).(type)},{obj.folderStruct(dayHasAFileOfType).(mask)});
                if isempty(selectedFiles)
                    warning('All files processed, set showProcessed to 1 to reprocess files')
                end
            end
            
        end
        
        function preprocess(obj)
            obj.buildFolderStructure;
            if obj.askForFiles
                selectedFiles = selectFiles(obj,'concat','downsample');
            else
                selectedFiles = getAllFilesWithType(obj,'concat','downsample');
            end
            
            if isempty(selectedFiles);warning('No files selected');return;end
            
            disp(repmat('-',[1,100]))
            disp(['Starting the preprocess of your files'])
            disp(repmat('-',[1,100]))
            
            for day = 1:length(selectedFiles)
                disp(repmat('#',[1,100]))
                disp(['day ', num2str(day), ' of ', num2str(length(selectedFiles))])
                obj.preprocessConcatDay(selectedFiles{day})
            end
            obj.buildFolderStructure;
        end
        
        function extract(obj)
            obj.buildFolderStructure;
            if obj.askForFiles
                selectedFiles = selectFiles(obj,'downsample','analysis');
            else
                selectedFiles = getAllFilesWithType(obj,'downsample','analysis');
            end
            
            if isempty(selectedFiles);warning('No files selected');return;end
            
            disp(repmat('-',[1,100]))
            disp(['Starting the extraction of your files'])
            disp(repmat('-',[1,100]))
            
            for day = 1:length(selectedFiles)
                disp(repmat('#',[1,100]))
                disp(['day ', num2str(day), ' of ', num2str(length(selectedFiles))])
                for f = 1:length(selectedFiles{day})
                    disp(repmat('.',[1,50]))
                    disp(['file ', num2str(f), ' of ', num2str(length(selectedFiles{day}))])
                    file = selectedFiles{day}(f);
                    fileForExtraction = [obj.rootFolder,filesep,file.folder,filesep,file.fileName];
                    switch obj.analysisType
                        case 'cellmax'
                           carpaUtilities.cellMaxExtraction(fileForExtraction,file.experiment);
                        case 'pcaica'
                           carpaUtilities.pcaicaExtraction(fileForExtraction,file.experiment);
                        otherwise
                           error('options are either cellmax or pcaica')
                    end
                end
            end
            obj.buildFolderStructure;                        
        end

        
        function filter(obj)
            obj.buildFolderStructure;
            if obj.askForFiles
                selectedFiles = selectFiles(obj,'analysis','decisions');
            else
                selectedFiles = getAllFilesWithType(obj,'analysis','decisions');
            end
            if isempty(selectedFiles);warning('No files selected');return;end
            
            tempMovieFileNames = cellfun(@(x) {x.fileName}, getAllFilesWithType(obj,'downsample'),'UniformOutput',0);
            allMovieFileNames = cat(2,tempMovieFileNames{:});
            
            cInfObjs = [];
            disp(repmat('-',[1,100]))
            disp('Initializing the manual sort for all the selected files, this may take a while...')
            disp(repmat('-',[1,100]))
            for day = 1:length(selectedFiles)
                disp(repmat('#',[1,100]))
                disp(['day ', num2str(day), ' of ', num2str(length(selectedFiles))])
                for f = 1:length(selectedFiles{day})
                    disp(repmat('.',[1,50]))
                    disp(['file ', num2str(f), ' of ', num2str(length(selectedFiles{day}))])
                    file = selectedFiles{day}(f);
                    
                    analysisFile = [obj.rootFolder,filesep,file.folder,filesep,file.fileName];
                    movieName = allMovieFileNames(carpaUtilities.compareNamesTillExperiment(file.fileName,allMovieFileNames,file.experiment));
                    if isempty(movieName)
                        warning(['No movie found for file ', file.fileName, ' skipping...'])
                        continue;
                    end
                    if length(movieName) > 1
                        warning(['Two compatible files found for file ', file.fileName, 'skipping...'])
                        continue;                        
                    end
                    movieFile = [obj.rootFolder,filesep,file.folder,filesep,movieName{1}];
                        
                    if ~isempty(find(carpaUtilities.checkRegexp({file.fileName},obj.analysisFileRegexp(1)))) && strcmp(obj.analysisType,'cellmax') %em
                        analysisEM = load(analysisFile);
                        cInf = cellInfo(movieFile,permute(analysisEM.emAnalysisOutput.cellImages,[3,1,2]),analysisEM.emAnalysisOutput.scaledProbability,[]);
                        typeAnalysis = 'em';
                        clear analysisEM
                    elseif ~isempty(find(carpaUtilities.checkRegexp({file.fileName},obj.analysisFileRegexp(2)))) && strcmp(obj.analysisType,'pcaica')%ica
                        analysisICA = load(analysisFile);
                        cInf = cellInfo(movieFile,permute(analysisICA.pcaicaAnalysisOutput.IcaFilters,[3,1,2]),analysisICA.pcaicaAnalysisOutput.IcaTraces,[]);                        
                        typeAnalysis = 'ica';
                        clear analysisICA
                    else
                        continue;
                    end
                    
                    cInf.predictClassifier('probTresh',.95);
                    cInfObjs = [cInfObjs,cInf];
                    
                end
            end
            
            disp(repmat('-',[1,100]))
            disp(['Starting with the manual sort for the ',num2str(length(cInfObjs)) ,' files'])
            disp(repmat('-',[1,100]))
            for k = 1:length(cInfObjs)
                disp(['MANUAL SORT FOR FILE ',num2str(k),' of ',num2str(length(cInfObjs))])
                try
                    cInfObjs(k).manualSort('tresholds',{'getOverlap','>=.4','getglobalSNR','>1'});
                catch
                    warning('Something went wrong, skiping file...')
                    continue;
                end
                cInfObjs(k).validCells(cInfObjs(k).validCells == 3) = 0;
                cInfObjs(k).saveDecisions(typeAnalysis,'skipConfirmation',1);
            end
            
            obj.buildFolderStructure;
        end

        function postprocess(obj)
            obj.buildFolderStructure;
            if obj.askForFiles
                selectedFiles = selectFiles(obj,'decisions','tracesEvents');
            else
                selectedFiles = getAllFilesWithType(obj,'decisions','tracesEvents');
            end
            
            if isempty(selectedFiles);warning('No files selected');return;end

            tempDFOFFiles = getAllFilesWithType(obj,'dfof');
            if isempty(tempDFOFFiles);error('NO DFOF FILES FOUND');end
            tempDFOFFileNames = cellfun(@(x) {x.fileName}, tempDFOFFiles,'UniformOutput',0);
            allDFOFFileNames = cat(2,tempDFOFFileNames{:});
            if isempty(allDFOFFileNames);error('There are no dfof files in folder');end

            tempBEHAVIORFiles = getAllFilesWithType(obj,'behavior');
            if isempty(tempBEHAVIORFiles);error('NO BEHAVIOR FILES FOUND');end
            tempBEHAVIORFileNames = cellfun(@(x) {x.fileName}, tempBEHAVIORFiles,'UniformOutput',0);
            allBEHAVIORFileNames = cat(2,tempBEHAVIORFileNames{:});
            if isempty(allBEHAVIORFileNames);error('There are no behavior files in folder');end

           tempLOGFiles = getAllFilesWithType(obj, 'logs');
            if isempty(tempLOGFiles)
                % If no log files found, create a new log file in the same folder as 'obj'
                extractedPart = '';
                pattern = '(?<=^.*?-.*?)(\d{8}_\d{6})(?=-)';  % Original pattern for basic case
                
                for day = 1:length(selectedFiles)
                    % Apply regexp to extract the desired part
                    disp(class(selectedFiles{day}.fileName));
                    
                    % First, try to match the pattern for the first date-time string
                    matches = regexp(selectedFiles{day}.fileName, pattern, 'match');
                    
                    % If no match, try a more complex pattern for files with multiple time stamps
                    if isempty(matches)
                        pattern = '(?<=^.*?-.*?)(\d{8}_\d{6}(&\d{6})*)(?=-)';  % Modified pattern to account for multiple time components
                        matches = regexp(selectedFiles{day}.fileName, pattern, 'match');
                    end
                    
                    % Check if the match was found
                    if ~isempty(matches)
                        extractedPart = matches{1};  % Extract the part with date and time(s)
                        disp(extractedPart);  % Output: 20220524_100253 or multiple time stamps
                        
                        % Now handle cases where there might be multiple time components
                        timeParts = strsplit(extractedPart, '&');  % Split by '&' to get each time component
                        
                        % For simplicity, we'll take the first time component for date-time parsing
                        firstTimePart = timeParts{1};
                        
                        % Extract date and time from the matched string (we'll only use the first time component here)
                        dateStr = firstTimePart(1:8);  % '20220524'
                        timeStr = firstTimePart(10:end);  % '100253'
                        
                        % Parse the date and time into separate components
                        year = str2double(dateStr(1:4));  % '2022'
                        month = str2double(dateStr(5:6));  % '05'
                        dayStr = str2double(dateStr(7:8));  % '24'
                        
                        hour = str2double(timeStr(1:2));  % '10'
                        minute = str2double(timeStr(3:4));  % '02'
                        second = str2double(timeStr(5:6));  % '53'
                        
                        % Create the final date-time array
                        dateTimeArray = [year, month, dayStr, hour, minute, second];
            
                        % Create the log entry in the appropriate structure
                        obj.folderStruct(day).logs(end+1).date = dateTimeArray;
                        obj.folderStruct(day).logs(end).folder = obj.folderStruct.folderName;
                        obj.folderStruct(day).logs(end).experiment = obj.folderStruct.experiment;
                        obj.folderStruct(day).logs(end).fileName = sprintf('Mouse-%s-%s-%s-log.xml', obj.folderStruct.mouse, extractedPart, obj.folderStruct.experiment);
                        % obj.folderStruct(day).logs(end).dispName = ['log_','day:',dateStr,'_sess:',char(join(sessions,'-'))];
                        obj.folderStruct(day).logs(end).dispName = ['log_','day:',dateStr,'_sess:',timeStr];
                    else
                        disp('No date and time string found');
                    end
            
                    disp('No log files found. A new log file has been created.');
                end
            end

            pause(1);  % Pause for 1 second to ensure file creation
            tempLOGFiles = getAllFilesWithType(obj, 'logs');
            tempLOGFileNames = cellfun(@(x) {x.fileName}, tempLOGFiles,'UniformOutput',0);
            allLOGFileNames = cat(2,tempLOGFileNames{:});
            if isempty(allLOGFileNames);error('There are no log files in folder');end

            tempEXTRACTFiles = getAllFilesWithType(obj,'analysis');
            if isempty(tempEXTRACTFiles);error('NO EXTRACT FILES FOUND');end
            tempEXTRACTFileNames = cellfun(@(x) {x.fileName}, tempEXTRACTFiles,'UniformOutput',0);
            allEXTRACTFileNames = cat(2,tempEXTRACTFileNames{:});
            if isempty(allEXTRACTFileNames);error('There are no analysis files in folder');end
            
            ALL = struct;
            
            disp(repmat('-',[1,100]))
            disp('Preparing all the files...')
            disp(repmat('-',[1,100]))
            for day = 1:length(selectedFiles)
                disp(repmat('#',[1,100]))
                disp(['day ', num2str(day), ' of ', num2str(length(selectedFiles))])
                for f = 1:length(selectedFiles{day})
                    disp(repmat('.',[1,50]))
                    disp(['file ', num2str(f), ' of ', num2str(length(selectedFiles{day}))])
                    file = selectedFiles{day}(f);
                    
                    decisionsFile = [obj.rootFolder,filesep,file.folder,filesep,file.fileName];
                    
                    dfofName = allDFOFFileNames(carpaUtilities.compareNamesTillExperiment(file.fileName,allDFOFFileNames,file.experiment));
                    if isempty(dfofName)
                        warning(['No dfof movie found for file ', file.fileName, 'skipping...'])
                        continue;
                    end
                    if length(dfofName) > 1
                        warning(['Two compatible files found for file ', file.fileName, 'skipping...'])
                        continue;                        
                    end
                    dfofFile = [obj.rootFolder,filesep,file.folder,filesep,dfofName{1}];
                    
                    extractName = allEXTRACTFileNames(carpaUtilities.compareNamesTillExperiment(file.fileName,allEXTRACTFileNames,file.experiment));
                    if isempty(extractName)
                        warning(['No extraction file found for file ', file.fileName, 'skipping...'])
                        continue;
                    end
                    if length(extractName) == 2
                        emFile = find(carpaUtilities.checkRegexp(extractName,obj.analysisFileRegexp(1)));
                        icaFile = find(carpaUtilities.checkRegexp(extractName,obj.analysisFileRegexp(2)));
                        if strcmp(obj.analysisType,'cellmax') %em
                            extractName = extractName(emFile);
                        elseif strcmp(obj.analysisType,'pcaica')%ica
                            extractName = extractName(icaFile);
                        end
                    end
                    if length(extractName) > 1
                        warning(['Two compatible files found for file ', file.fileName, 'skipping...'])
                        continue;                        
                    end
                    extractFile = [obj.rootFolder,filesep,file.folder,filesep,extractName{1}];
                    
                    %LOGS
                    logNames = {};
                    for k = 1:size(file.date,1)
                        lName = ['Mouse-',obj.mouse,'-',datestr(file.date(k,:),obj.dateParser),'_',datestr(file.date(k,:),obj.timeParser),'-',file.experiment,'-log'];
                        logNames = [logNames,allLOGFileNames(carpaUtilities.compareNamesTillExperiment(lName,allLOGFileNames,file.experiment))];
                    end
                    if isempty(logNames)
                        disp(repmat('-',[1,100]))
                        disp('The sesions in this file are...')
                        disp(file.date)
                        disp(repmat('-',[1,100]))
                        h = msgbox(num2str(file.date),'SESSIONS IN FILE:','Help');
                        newLogNames = uigetfile({'*.txt;*.xml';},'Select all log sessions of this file: ',[obj.rootFolder,filesep,file.folder],'MultiSelect', 'on');
                        delete(h)
                        if ~iscell(newLogNames)
                            newLogNames = {newLogNames};
                        end
                        
                        if isempty(newLogNames)
                            warning(['No logs selected for file ', file.fileName, 'skipping...'])
                            continue;
                        end
                        
                        %Rename log files
                        if length(newLogNames) == size(file.date,1)
                            for k = 1:length(newLogNames)
                                ext = strsplit(newLogNames{k},'.');
                                saveName = ['Mouse-',obj.mouse,'-',datestr(file.date(k,:),obj.dateParser),'_',datestr(file.date(k,:),obj.timeParser),'-',file.experiment,'-log.',ext{end}];
                                movefile([obj.rootFolder,filesep,file.folder,filesep,newLogNames{k}],[obj.rootFolder,filesep,file.folder,filesep,saveName])
                            end
                            obj.buildFolderStructure;
                            tempLOGFileNames = cellfun(@(x) {x.fileName}, getAllFilesWithType(obj,'logs'),'UniformOutput',0);
                            allLOGFileNames = cat(2,tempLOGFileNames{:});
                            logNames = {};
                            for k = 1:size(file.date,1)
                                lName = ['Mouse-',obj.mouse,'-',datestr(file.date(k,:),obj.dateParser),'_',datestr(file.date(k,:),obj.timeParser),'-',file.experiment,'-log'];
                                logNames = [logNames,allLOGFileNames(carpaUtilities.compareNamesTillExperiment(lName,allLOGFileNames,file.experiment))];
                            end
                        else
                            logNames = newLogNames;
                        end
                    end
                    logFiles = cellfun(@(x) [obj.rootFolder,filesep,file.folder,filesep,x], logNames, 'UniformOutput',0);
                    
                    %% BEHAVIOR
                    behaviorNames = {};
                    for k = 1:size(file.date,1)
                        bName = ['Mouse-',obj.mouse,'-',datestr(file.date(k,:),obj.dateParser),'_',datestr(file.date(k,:),obj.timeParser),'-',file.experiment,'-behavior.avi'];
                        behaviorNames = [behaviorNames,allBEHAVIORFileNames(carpaUtilities.compareNamesTillExperiment(bName,allBEHAVIORFileNames,file.experiment))];
                    end
                    if isempty(behaviorNames)
                        disp(repmat('-',[1,100]))
                        disp('The sesions in this file are...')
                        disp(file.date)
                        disp(repmat('-',[1,100]))
                        h = msgbox(num2str(file.date),'SESSIONS IN FILE:','Help');
                        newBehaviorNames = uigetfile('*.avi','Select all behavior sessions of this file: ',[obj.rootFolder,filesep,file.folder],'MultiSelect', 'on');
                        if ~iscell(newBehaviorNames)
                            newBehaviorNames = {newBehaviorNames};
                        end
                        delete(h)
                        
                        if isempty(newBehaviorNames)
                            warning(['No behavior selected for file ', file.fileName, 'skipping...'])
                            continue;
                        end
                        
                        if length(newBehaviorNames) ~= length(logFiles)%size(file.date,1)
                            
                            box_text = {};
                            frames_logs = zeros([1,length(logFiles)]);
                            for log_k = 1:length(logFiles)
                                [~, frames_logs(log_k)] = readFramesFromLog(logFiles{log_k});
                                box_text = [box_text;[strrep(strrep(num2str(file.date(log_k,:)),'    ','.'),'..','.'),' : ',num2str(frames_logs(log_k))]];
                            end
                            
                            box_text = [box_text;'---------------'];
                            frames_behavior = zeros([1,length(newBehaviorNames)]);
                            for beh_k = 1:length(newBehaviorNames)
                                mov_temp = VideoReader([obj.rootFolder,filesep,file.folder,filesep,newBehaviorNames{beh_k}]);
                                frames_behavior(beh_k) = mov_temp.NumberOfFrames;
                                box_text = [box_text;[newBehaviorNames{beh_k},' : ',num2str(frames_behavior(beh_k))]];
                            end   
                            
                            h2 = msgbox(box_text,'FRAMES IN LOGS AND BEHAVIOR:','Help');
                            
                            %
                            
                            if length(newBehaviorNames) > size(file.date,1) %Merge behaviors
                                answer = questdlg('Merge or delete behavior files?','WARNING: More behavior than calcium frames','Merge','Delete','Skip','Merge');
                                behavior_name_list = cellfun(@(x) fullfile(obj.rootFolder,file.folder,x),newBehaviorNames, 'UniformOutput',0);
                                switch answer
                                    case 'Merge'
                                        new_behavFiles = carpaUtilities.merge_list(behavior_name_list,fullfile(obj.rootFolder,file.folder),obj.behaviorRegexp);
                                        behaviorFiles = new_behavFiles;
                                        
                                    case 'Delete'
                                         new_behavFiles = carpaUtilities.delete_list(behavior_name_list,fullfile(obj.rootFolder,file.folder),obj.behaviorRegexp);
                                         behaviorFiles = new_behavFiles;
                                        
                                    case 'Skip'
                                        continue
                                end
                                
                            elseif length(newBehaviorNames) < size(file.date,1) %Merge calcium
                                answer = questdlg('Merge or delete calcium files?','WARNING: More calcium than behavior frames','Merge','Delete','Skip','Merge');
                                switch answer
                                    case 'Merge'
                                        disp('Merging calcium')
                                        new_logFiles = carpaUtilities.merge_list(logFiles,fullfile(obj.rootFolder,file.folder),obj.logRegexp);
                                        logFiles = new_logFiles; 
                                        
                                    case 'Delete'
                                        deleted_logFiles = carpaUtilities.delete_list(logFiles,fullfile(obj.rootFolder,file.folder),obj.logRegexp);
                                        [~,idx] = setdiff(logFiles,deleted_logFiles);
                                        
                                        new_logFiles = logFiles;
                                        for k = 1:length(idx)
                                            new_logFiles{idx(k)} = {new_logFiles{idx(k)}};
                                        end
                                        logFiles = new_logFiles; 
                                        
                                    case 'Skip'
                                        continue
                                end 
                                behaviorFiles = cellfun(@(x) [obj.rootFolder,filesep,file.folder,filesep,x], newBehaviorNames, 'UniformOutput',0);
                            end
                            delete(h2)
                        else
                            %Rename behavior files

                            for k = 1:length(newBehaviorNames)
                                saveName = ['Mouse-',obj.mouse,'-',datestr(file.date(k,:),obj.dateParser),'_',datestr(file.date(k,:),obj.timeParser),'-',file.experiment,'-behavior.avi'];
                                movefile([obj.rootFolder,filesep,file.folder,filesep,newBehaviorNames{k}],[obj.rootFolder,filesep,file.folder,filesep,saveName])
                            end
                            obj.buildFolderStructure;
                            tempBEHAVIORFileNames = cellfun(@(x) {x.fileName}, getAllFilesWithType(obj,'behavior'),'UniformOutput',0);
                            allBEHAVIORFileNames = cat(2,tempBEHAVIORFileNames{:});
                            behaviorNames = {};
                            for k = 1:size(file.date,1)
                                bName = ['Mouse-',obj.mouse,'-',datestr(file.date(k,:),obj.dateParser),'_',datestr(file.date(k,:),obj.timeParser),'-',file.experiment,'-behavior.avi'];
                                behaviorNames = [behaviorNames,allBEHAVIORFileNames(carpaUtilities.compareNamesTillExperiment(bName,allBEHAVIORFileNames,file.experiment))];
                            end
                            behaviorFiles = cellfun(@(x) [obj.rootFolder,filesep,file.folder,filesep,x], behaviorNames, 'UniformOutput',0);
                        end
                    else
                        behaviorFiles = cellfun(@(x) [obj.rootFolder,filesep,file.folder,filesep,x], behaviorNames, 'UniformOutput',0);
                    end
                   
                    ALL(end+~isempty(fields(ALL))).extractFile = extractFile;
                    ALL(end).decisionsFile = decisionsFile;
                    ALL(end).dfofFile = dfofFile;
                    ALL(end).logFiles = logFiles;
                    ALL(end).behaviorFiles = behaviorFiles;
                
                end
            end
                          
            background_default = '-0.8';
            answer = inputdlg('Please set the background threshold', 'Background threshold', [1,35], {background_default});
            if isempty(answer)
                answer = {background_default};
            end
            backgroundTreshold = str2double(answer{1});
            
            disp(repmat('-',[1,100]))
            disp('Starting the postprocessing process...')
            disp(repmat('-',[1,100]))
            for k = 1:length(ALL)
                disp(repmat('#',[1,100]))
                disp(repmat('.',[1,50]))
                disp(['Postprocessing file ', num2str(k), ' of ', num2str(length(ALL))])
                %try
                carpaUtilities.getTracesEventsTraj(ALL(k).extractFile,ALL(k).decisionsFile,ALL(k).dfofFile,ALL(k).logFiles,ALL(k).behaviorFiles,selectedFiles{k}.experiment,backgroundTreshold)
                %catch ME
                %     disp(ME.message)
                %     warning('ERROR GENERATING TRACES EVENTS')
                %end
            end
            obj.buildFolderStructure;
        end
        
        function preprocessConcatDay(obj,files)
            joinedSessionsIdx = ones([1,length(files)]);%Sessions with the same number will be joined
            disp('Concatenating files now...');

            for k = 1:(length(files)-1)
                if etime(files(k+1).date,files(k).date)/60 > obj.joinFileTreshold
                    joinedSessionsIdx(k+1:end) = joinedSessionsIdx(k+1:end)+1;
                end
            end
            for k = unique(joinedSessionsIdx) %Process the joined sessions
                %disp(repmat('-',[1,50]))
                fileNums = joinedSessionsIdx == k;
                fileNames = {files(fileNums).fileName};
                disp(fileNames)
                fileFull = cellfun(@(x) [obj.rootFolder,filesep,files(1).folder,filesep,x],fileNames,'UniformOutput',false);
                parsedSessionDate = cellfun(@(x) datestr(x,obj.timeParser),{files(fileNums).date},'UniformOutput',false);

                saveName = ['Mouse-',obj.mouse,'-',datestr(files(1).date,obj.dateParser),'_',char(join(parsedSessionDate,'&')),'-',files(1).experiment];
                [movieDS,movieDFOF,~,~] = carpaUtilities.preprocessMovie(fileFull);
                carpaUtilities.createHDF5([obj.rootFolder,filesep,files(1).folder,filesep,saveName,'-dfof.h5'],movieDFOF)
                carpaUtilities.createHDF5([obj.rootFolder,filesep,files(1).folder,filesep,saveName,'-dfof-downsampled.h5'],movieDS)
            end
        end 
                        
        function downloadH5fromServer(obj)
           disp(['Downloading mouse',obj.mouse,' h5files from the server'])
           obj.sftp.cd(obj.archiveRawCaPath)
           animalFolders = obj.sftp.commandToVar('ls');
           scnsize = get(0,'ScreenSize');dlgSize = [scnsize(3)*0.2 scnsize(4)*0.7];
            correctAnimalFolder = listdlg('PromptString','Select animal to download:',...
        'SelectionMode','single',...
        'ListString',animalFolders,'ListSize',dlgSize);
               
           currentPath = obj.sftp.commandToVar('pwd');
           obj.sftp.cd(currentPath{1})
           obj.sftp.cd(animalFolders{correctAnimalFolder})
           dayFolders = obj.sftp.commandToVar('ls')';
           currentPath = obj.sftp.commandToVar('pwd');
           
           s = listdlg('PromptString','Select day(s) to download:',...
                'SelectionMode','multiple',...
                'ListString',dayFolders,'ListSize',dlgSize);
            
           for folder = dayFolders(s)
               if isempty(strfind(folder{1},'.')) %Is a folder, not a file
                   obj.sftp.cd(currentPath{1})
                   obj.sftp.cd(folder{1})
                   filesInFolder = obj.sftp.commandToVar('ls')';
                   fileExtensions = cellfun(@(x) x{end}, cellfun(@(x) strsplit(x,'.'),filesInFolder,'UniformOutput',0),'UniformOutput',0);
                   rawAndH5Files = find(~cellfun(@isempty,strfind(fileExtensions,'raw')) | ~cellfun(@isempty,strfind(fileExtensions,'h5')) | ~cellfun(@isempty,strfind(fileExtensions,'hdf5')));
                   %Get the different raw sessions
                   allSessions = struct('session',{},'files',{});
                   for file = filesInFolder(rawAndH5Files)
                       fileSession = regexp(file{1}, 'recording.*\d+_(\d+)','tokens');
                       if ~isempty(fileSession)
                           sessionIdx = find(strcmp(fileSession{1}{1},{allSessions(:).session}));
                           if isempty(sessionIdx)
                               allSessions(length(allSessions) + 1).session = fileSession{1}{1};
                               allSessions(length(allSessions)).files = file{1};
                           else
                               allSessions(sessionIdx).files = [allSessions(sessionIdx).files, file];
                           end
                       end
                   end
                   for sessionIdx = 1:length(allSessions)
                       destinationPath = [obj.rootFolder,filesep,folder{1}];
                       
                       mkdir(destinationPath)
                       
                       if ~iscell(allSessions(sessionIdx).files)
                           allSessions(sessionIdx).files = {allSessions(sessionIdx).files};
                       end
                       sessionFileExtensions = cellfun(@(x) x{end}, cellfun(@(x) strsplit(x,'.'),allSessions(sessionIdx).files,'UniformOutput',0),'UniformOutput',0);
                       sessionHdf5files = find(~cellfun(@isempty,strfind(sessionFileExtensions,'h5')) | ~cellfun(@isempty,strfind(sessionFileExtensions,'hdf5')));
                       allH5FilesPerSession = allSessions(sessionIdx).files(sessionHdf5files);

                       if length(allH5FilesPerSession) > 1
                            namesWithoutConcat = cellfun(@(x) erase(x,'concat_'),allH5FilesPerSession,'UniformOutput',0);
                            if isequal(namesWithoutConcat{:})
                                sessionHdf5files = find(~cellfun(@isempty,strfind(allH5FilesPerSession,'concat_')));
                            end
                       end
                       
                       if length(sessionHdf5files) == 1
                           if isempty(dir([destinationPath,filesep,'concat_',allSessions(sessionIdx).files{sessionHdf5files}]))%file already processed
                               obj.sftp.download(allSessions(sessionIdx).files{sessionHdf5files},destinationPath);
                               fileName = allSessions(sessionIdx).files{sessionHdf5files};
                               destinationName = [destinationPath,filesep,fileName];

                               %Check if the file is valid
                               try
                                   hinf = hdf5info(destinationName);
                               catch
                                   hinf = [];
                               end

                               if ~isempty(hinf)
                                   %Check if is decompressed, decompress if not
                                   if isempty(strfind(fileName,'concat'))

                                       datasetIdx = find(cellfun(@(x) length(x), {hinf.GroupHierarchy.Datasets(:).Dims}) == 3);

                                       sizeLimit = 10*1e9;
                                       fileInfo = dir(destinationName);

                                       if fileInfo.bytes >  sizeLimit %Read part by part
                                           origDims = hinf.GroupHierarchy.Datasets(datasetIdx).Dims;
                                           sliceVect = ceil([0,origDims(3)/(fileInfo.bytes/sizeLimit)*(1:ceil((fileInfo.bytes/sizeLimit)))]);
                                           sliceVect(end) = origDims(3);

                                           datasetDS = [];
                                           for slice = 1:length(sliceVect)-1
                                               dataset = h5read(destinationName,hinf.GroupHierarchy.Datasets(datasetIdx).Name,[1 1 sliceVect(slice)+1],[origDims(1) origDims(2) sliceVect(slice+1)-sliceVect(slice)]);
                                               datasetDS = cat(3,datasetDS,carpaUtilities.downsampleSpace(double(dataset),obj.spatialDS));
                                           end

                                       else
                                            dataset = hdf5read(hinf.GroupHierarchy.Datasets(datasetIdx));
                                            datasetDS = carpaUtilities.downsampleSpace(double(dataset),obj.spatialDS);
                                       end
                                       clear dataset
                                       try
                                          carpaUtilities.createHDF5([destinationPath,filesep,'concat_',fileName],datasetDS)
                                          delete(destinationName)
                                       catch ME
                                          disp(ME.message)
                                       end
                                   end
                                   continue
                               else
                                   delete(destinationName)
                                   allSessions(sessionIdx).files(sessionHdf5files) = [];
                               end
                           else
                               warning(['File ',[destinationPath,filesep,'concat_',allSessions(sessionIdx).files{sessionHdf5files}],' already processed'])
                               continue
                           end
                       elseif length(sessionHdf5files) > 1
                           error(['Two hdf5 for one sesison in ', allSessions(sessionIdx).files{sessionHdf5files}])
                       end
                       
                       namePart = strsplit(allSessions(sessionIdx).files{1},'-');
                       if ~isempty(dir([destinationPath,filesep,'concat_',namePart{1},'.*']))
                           continue
                       end
                       downloadFiles = cell(size(allSessions(sessionIdx).files));
                       i = 1;
                       for file = allSessions(sessionIdx).files
                           downloadFiles{i} = ['''',file{1},''''];
                           i = i + 1;
                       end
                       obj.sftp.download(downloadFiles,destinationPath);

                       destinationFiles = allSessions(sessionIdx).files;
                       filesInOrder = [length(destinationFiles),1:(length(destinationFiles)-1)];
                       destinationFiles = destinationFiles(filesInOrder);
                       [~,firstRawName] = fileparts(destinationFiles{1});
                       concatName = [destinationPath,filesep,'concat_',firstRawName,'.hdf5'];
                       datasetDS = []; 
                       for k = 1:length(destinationFiles)
                           rawFile = [destinationPath,filesep,destinationFiles{k}];
                           [~,rawName] = fileparts(rawFile);
                           hdf5FileName = [destinationPath,filesep,rawName,'.hdf5'];
                           system(['InscopixDecompress.bat ',rawFile,' ',hdf5FileName])

                           hdf5Files = dir([destinationPath,filesep,rawName,'*.hdf5']);
                           fullDataset = [];
                           for hdf5File = hdf5Files'
                               hinf = hdf5info([hdf5File.folder,filesep,hdf5File.name]);
                               dataset = hdf5read(hinf.GroupHierarchy.Datasets(find(cellfun(@(x) length(x), {hinf.GroupHierarchy.Datasets(:).Dims}) == 3)));
                               fullDataset = cat(3,fullDataset,carpaUtilities.downsampleSpace(double(dataset),obj.spatialDS));
                               clear dataset
                           end
                           datasetDS = cat(3,datasetDS,fullDataset);
                           clear fullDataset
                           cellfun(@(x) delete([destinationPath,filesep,x]),{hdf5Files.name})
                       end
                       carpaUtilities.createHDF5(concatName,datasetDS)
                       cellfun(@(x) delete([destinationPath,filesep,x]),allSessions(sessionIdx).files)
                       clear datasetDS
                   end     
               end
           end
        end
        
        function buildFolderStructure(obj)
            allItems = dir(obj.rootFolder);
            allItems = allItems(3:end); %remove WINDOWS metadata folders
            if isempty(allItems)
                obj.downloadH5fromServer;
                allItems = dir(obj.rootFolder);
                allItems = allItems(3:end); %remove WINDOWS metadata folders
            end
            allFolders = allItems([allItems.isdir]);
            obj.folderStruct = [];
            i = 0;
            for dayFolder = {allFolders.name}
                i = i + 1;
                allDigits = regexp(dayFolder{1}, '(\d+)','tokens');
                day = allDigits{end}{1};%Assumess date is the last sequence of numbers in the folder name
                fullPath = fullfile(obj.rootFolder, dayFolder{1});    
                                
                
                obj.folderStruct(i).path = fullPath;
                obj.folderStruct(i).folderName = dayFolder{1};
                obj.folderStruct(i).date = datevec(day,obj.dateParser);
                obj.folderStruct(i).mouse = obj.mouse;
                
                experiment = char(join(cellfun(@(x) x, regexp(dayFolder{1},'-([a-z | A-Z | _]+)','tokens')),'-'));
                obj.folderStruct(i).experiment = experiment;
                
                allInsideitems = dir(fullPath);
                allInsideitems = allInsideitems(3:end); %remove WINDOWS metadata folders
                allFiles = allInsideitems(~[allInsideitems.isdir]);
                
                obj.folderStruct(i).concat = [];
                obj.folderStruct(i).dfof = [];
                obj.folderStruct(i).downsample = [];
                obj.folderStruct(i).analysis = [];
                obj.folderStruct(i).decisions = [];
                obj.folderStruct(i).behavior = [];
                obj.folderStruct(i).logs = [];
                obj.folderStruct(i).tracesEvents = [];
                obj.folderStruct(i).other = [];
                
                for k = 1:length(allFiles)
                    
                    try
                        noExtensionName = split(allFiles(k).name,'.');
                        allDigits = regexp(noExtensionName{1}, '(\d+)','tokens');
                        sessions = cellfun(@(x) x, allDigits(3:end)); 
                        sMat = [];
                        for s = 1:length(sessions)
                            sMat = [sMat;datevec([day,sessions{s}],[obj.dateParser,obj.timeParser])];
                        end
                    end
                    
                    if carpaUtilities.checkRegexp({allFiles(k).name},obj.concatRegexp)
                        
                        if isempty(sessions)%usually concat has no mouse code at the begining
                            sessions = allDigits{end};
                            sMat = datevec([day,sessions{1}],[obj.dateParser,obj.timeParser]);
                        end
                        
                        obj.folderStruct(i).concat(end+1).date = sMat;
                        obj.folderStruct(i).concat(end).folder = dayFolder{1};
                        obj.folderStruct(i).concat(end).experiment = experiment;
                        obj.folderStruct(i).concat(end).fileName = allFiles(k).name; 
                        obj.folderStruct(i).concat(end).dispName =  ['concat_','day:',day,'_sess:',char(join(sessions,'-'))];
                    
                    elseif carpaUtilities.checkRegexp({allFiles(k).name},obj.dfofRegexp)

                        obj.folderStruct(i).dfof(end+1).date = sMat;
                        obj.folderStruct(i).dfof(end).folder = dayFolder{1};
                        obj.folderStruct(i).dfof(end).experiment = experiment;
                        obj.folderStruct(i).dfof(end).fileName = allFiles(k).name;
                        obj.folderStruct(i).dfof(end).dispName = ['dfof','day:',day,'_sess:',char(join(sessions,'-'))];
     
                    elseif carpaUtilities.checkRegexp({allFiles(k).name},obj.downsampleRegexp)

                        obj.folderStruct(i).downsample(end+1).date = sMat;
                        obj.folderStruct(i).downsample(end).folder = dayFolder{1};
                        obj.folderStruct(i).downsample(end).experiment = experiment;
                        obj.folderStruct(i).downsample(end).fileName = allFiles(k).name;
                        obj.folderStruct(i).downsample(end).dispName = ['downsampled_','day:',day,'_sess:',char(join(sessions,'-'))];
                                                
                    elseif carpaUtilities.checkRegexp({allFiles(k).name},obj.analysisFileRegexp)
                        
                        obj.folderStruct(i).analysis(end+1).date = sMat;
                        obj.folderStruct(i).analysis(end).folder = dayFolder{1};
                        obj.folderStruct(i).analysis(end).experiment = experiment;
                        obj.folderStruct(i).analysis(end).fileName = allFiles(k).name;
                        obj.folderStruct(i).analysis(end).dispName = ['analysis_','day:',day,'_sess:',char(join(sessions,'-'))];
               
                        
                    elseif carpaUtilities.checkRegexp({allFiles(k).name},obj.decisionsRegexp)
                        
                        obj.folderStruct(i).decisions(end+1).date = sMat;
                        obj.folderStruct(i).decisions(end).folder = dayFolder{1};
                        obj.folderStruct(i).decisions(end).experiment = experiment;
                        obj.folderStruct(i).decisions(end).fileName = allFiles(k).name;
                        obj.folderStruct(i).decisions(end).dispName = ['decisions_','day:',day,'_sess:',char(join(sessions,'-'))];
               
                        
                    elseif carpaUtilities.checkRegexp({allFiles(k).name},obj.behaviorRegexp)
                        
                        obj.folderStruct(i).behavior(end+1).date = sMat;
                        obj.folderStruct(i).behavior(end).folder = dayFolder{1};
                        obj.folderStruct(i).behavior(end).experiment = experiment;
                        obj.folderStruct(i).behavior(end).fileName = allFiles(k).name;
                        obj.folderStruct(i).behavior(end).dispName = ['behavior_','day:',day,'_sess:',char(join(sessions,'-'))];
               
                        
                    elseif carpaUtilities.checkRegexp({allFiles(k).name},obj.logRegexp)
                        
                        obj.folderStruct(i).logs(end+1).date = sMat;
                        obj.folderStruct(i).logs(end).folder = dayFolder{1};
                        obj.folderStruct(i).logs(end).experiment = experiment;
                        obj.folderStruct(i).logs(end).fileName = allFiles(k).name;
                        obj.folderStruct(i).logs(end).dispName = ['log_','day:',day,'_sess:',char(join(sessions,'-'))];
               
                    
                    elseif carpaUtilities.checkRegexp({allFiles(k).name},obj.tracesEventsRegexp)
                        
                        obj.folderStruct(i).tracesEvents(end+1).date = sMat;
                        obj.folderStruct(i).tracesEvents(end).folder = dayFolder{1};
                        obj.folderStruct(i).tracesEvents(end).experiment = experiment;
                        obj.folderStruct(i).tracesEvents(end).fileName = allFiles(k).name;
                        obj.folderStruct(i).tracesEvents(end).dispName = ['tracesEvents_','day:',day,'_sess:',char(join(sessions,'-'))];
                        
                    else 
                        obj.folderStruct(i).other(end+1).fileName = allFiles(k).name;
                        obj.folderStruct(i).other(end).experiment = experiment;
                        obj.folderStruct(i).other(end).folder = dayFolder{1};
                    end
                end
                
                %Sort concat in descending time order
                if ~isempty(obj.folderStruct(i).concat)
                    [~,sortedIdx] = sort(cellfun(@datenum,{obj.folderStruct(i).concat.date}));
                    obj.folderStruct(i).concat = obj.folderStruct(i).concat(sortedIdx);
                end
                
            end
        end
        
        function setExperimentNames(obj)
            
            if isempty(obj.rootFolder)
                obj.rootFolder = uigetdir('','Select the folder with all the sessions');
            end
            
            if isempty(obj.mouse)
                regexpOut = regexp(obj.rootFolder, 'ouse\D*(\d+)','tokens');
                try mouseCode = regexpOut{1}{1};catch;mouseCode = '0000';end
                mouseCode = inputdlg('Enter the mouse code:','Mouse code',1,{mouseCode});
                obj.mouse = mouseCode{1};
            end
            
            allFolders = dir(obj.rootFolder);
            
            allFolders = allFolders(3:end);
            allFolders = allFolders([allFolders.isdir]);
            
           foldersPerExperiment = struct;
           k = 1;
           scnsize = get(0,'ScreenSize');dlgSize = [scnsize(3)*0.2 scnsize(4)*0.7];
           answer = questdlg('Does the mouse perform more than one experiment?', ...
           'Experiment data', ...
           'Yes','No','No');
            switch answer
                case 'No'
                    experimentName = inputdlg('Enter experiment name:','Experiment name',1);if isempty(experimentName);return;end
                    foldersPerExperiment(k).name = experimentName{1};
                    foldersPerExperiment(k).folders = {allFolders.name};
                case 'Yes'
                    experimentName = inputdlg('Enter a first experiment name:','Experiment name',1);if isempty(experimentName);return;end
                    correctFolders = {allFolders.name};
                    while true
                       foldersPerExperiment(k).name = experimentName{1};
                        
                       foldersIdx = listdlg('PromptString','Select all folders for this experiment','SelectionMode','multiple','ListString',correctFolders,'ListSize',dlgSize);
                       folders = correctFolders(foldersIdx);
                       foldersPerExperiment(k).folders = folders;
                       correctFolders = setdiff(correctFolders,foldersPerExperiment(k).folders);
                       if isempty(correctFolders)
                           break;
                       end
                       experimentName = inputdlg('Enter another experiment name:','Experiment name',1);
                       k = k + 1;
                    end
                otherwise
                    return;                    
            end
            
            for k = 1:length(foldersPerExperiment)
                for folder = foldersPerExperiment(k).folders
                    nameParts = char(join(cellfun(@(x) x, regexp(folder{1},'-([a-z | A-Z]+)','tokens')),'-'));
                    noExpName = erase(folder{1},nameParts);
                    newName = [noExpName,foldersPerExperiment(k).name];
                    movefile([obj.rootFolder,filesep,folder{1}],[obj.rootFolder,filesep,newName])
                end
            end
            
        end
         
    end
    
end