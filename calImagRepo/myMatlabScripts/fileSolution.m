classdef fileSolution < handle
% insert definition here
    properties
        %FOLDER INFORMATION
        rootFolder = '';
        experiment = '';
        folderName = '';
        mouse = '0000';
        folderStruct = [];
        
        %OTHER OBJECT INFO
        sftp = '';
        %Time parsers
        dateParser = 'yyyymmdd'; %Reads and prints time with this format
        timeParser = 'HHMMSS';
        
        %OPTIONS
        decideIfMissing = 1;
        extractIfMissing = 1;
        processIfMissing = 1;
        postProcessIfMissing = 1;
        joinFileTreshold = 60; %Will join the concat files if they are closer than treshold (minutes). To not join put 0.
        analysisType = 'cellmax';%'cellmax'
        spatialDS = 4;
        
        %REGEXP STUFF FOR FINDING FILES       
        analysisFileRegexp = {'emAnalysis_?\d*.mat$','pcaicaAnalysis_?\d*.mat$'};
        downsampleRegexp = {'downsampled_?\d*.h5$'};
        dfofRegexp = {'dfof_?\d*.h5$'};
        concatRegexp = {'concat_recording.*\d+_\d+.hdf5$','concat_recording.*\d+_\d+.h5$'};%{'concat_recording.*.h5$','concat_recording.*.hdf5$'};
        decisionsRegexp = {'Sorted_?\d*.mat$','decisions_?\d*.mat$'}
        behaviorRegexp = {'.avi$'}
        logRegexp = {'.xml$','.txt$'}
        tracesEventsRegexp = {'TracesAndEvents?\d*.mat$'}
        
        %SERVER PATHS
        archiveBehaviorPath = '/archive/pjercog/Raw_Data_CalcImaging/Behavior';
        archiveRawCaPath = '/archive/pjercog/Raw_Data_CalcImaging';
    end
    
    methods
        function obj = fileSolution
            %SAVE EXPERIMENT, FOLDER AND MOUSE
%             obj.experiment = 'linear_track';%'Global_remapping';%
%             obj.rootFolder = 'D:\Storage\Processing\Mouse-2023 - test';%'E:\Processing\Mouse2019\';%
            
            obj.rootFolder = uigetdir('','Select the folder to process: ');        
            prompt = {'Experiment name:'};
            dlg_title = 'Enter the experiment name: ';
            num_lines = 1;
            defaultans = {'null'};
            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            
            if isempty(answer)
                obj.experiment = 'icx';
            else
                obj.experiment = answer{1};
            end
            
            nameSep = strsplit(obj.rootFolder,filesep);nameSep = nameSep(~cellfun(@isempty,nameSep));
            obj.folderName = nameSep{end};
            regexpOut = regexp(obj.folderName, 'ouse\D*(\d+)','tokens');
            try
                obj.mouse = regexpOut{1}{1};
            end
            disp(['Mouse ',obj.mouse,' (run changeNames method to modify)'])
            
            %INTITATE SERVER OBJECT
            try
                obj.sftp = mysftp;
            catch ME
                disp(ME.message)
                disp('server not initiated')
            end
            
            %LOAD REG FILES INSIDE THE FOLDER
            obj.getLocalConcatFiles;
            
        end
        
        function changeNames(obj)
            prompt = {'Enter experiment name','Enter mouse number'};
            dlg_title = 'Enter the new names: ';
            num_lines = 1;
            defaultans = {obj.experiment,obj.mouse};
            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            obj.experiment = answer{1};
            obj.mouse = answer{2};
        end
        
        function copyOnlyExtensionFiles(obj,extensions)
            if nargin < 2
                extensions = {'avi','mat'};
            end
            if ~iscell(extensions)
                extensions = {extensions};
            end
            folders = obj.chooseDays;
            copyRootPath = uigetdir('','Select the root folder to copy the data');
            copyFolderName = [obj.folderName,'_data'];
            for folder = folders
                for extension = extensions
                    allItems = dir([folder{1}.path,filesep,'*.',extension{1}]);
                    destinationPath = [copyRootPath,filesep,copyFolderName,filesep,folder{1}.folderName];
                    mkdir(destinationPath)
                    for k = 1:length(allItems)
                        copyfile([allItems(k).folder,filesep,allItems(k).name],[destinationPath,filesep,allItems(k).name])
                    end
                end
            end
        end
        
        function downloadH5fromServer(obj)
           disp(['Downloading mouse',obj.mouse,' h5files from the server'])
           obj.sftp.cd(obj.archiveRawCaPath)
           animalFolders = obj.sftp.commandToVar('ls');
           
           correctAnimalFolder = find(obj.checkRegexp(animalFolders',{obj.mouse}));
           if isempty(correctAnimalFolder)
               error('Could not find the logs, animal missing')
           elseif length(correctAnimalFolder) > 1
                s = listdlg('PromptString','Warning,duplicate folders, select the appropiate one',...
                            'SelectionMode','single',...
                            'ListString',animalFolders(correctAnimalFolder)');
                correctAnimalFolder = correctAnimalFolder(s);
           end
           currentPath = obj.sftp.commandToVar('pwd');
           obj.sftp.cd(currentPath{1})
           obj.sftp.cd(animalFolders{correctAnimalFolder})
           dayFolders = obj.sftp.commandToVar('ls')';
           currentPath = obj.sftp.commandToVar('pwd');
           
           s = listdlg('PromptString','Select day(s) to download:',...
                'SelectionMode','multiple',...
                'ListString',dayFolders);
            
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
                                               datasetDS = cat(3,datasetDS,fileSolution.downsampleSpace(double(dataset),obj.spatialDS));
                                           end

                                       else
                                            dataset = hdf5read(hinf.GroupHierarchy.Datasets(datasetIdx));
                                            datasetDS = fileSolution.downsampleSpace(double(dataset),obj.spatialDS);
                                       end
                                       clear dataset
                                       try
                                          obj.createHDF5([destinationPath,filesep,'concat_',fileName],datasetDS)
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
                               fullDataset = cat(3,fullDataset,fileSolution.downsampleSpace(double(dataset),obj.spatialDS));
                               clear dataset
                           end
                           datasetDS = cat(3,datasetDS,fullDataset);
                           clear fullDataset
                           cellfun(@(x) delete([destinationPath,filesep,x]),{hdf5Files.name})
                       end
                       obj.createHDF5(concatName,datasetDS)
                       cellfun(@(x) delete([destinationPath,filesep,x]),allSessions(sessionIdx).files)
                       clear datasetDS
%                        cellfun(@(y) cellfun(@(x) delete([destinationPath,filesep,x]),y),{allSessions(sessionIdx).files})
                   end
                                       
%                    hdf5files = find(~cellfun(@isempty,strfind(fileExtensions,'h5')) | ~cellfun(@isempty,strfind(fileExtensions,'hdf5')));
%                    logFiles = find(~cellfun(@isempty,strfind(fileExtensions,'txt')) | ~cellfun(@isempty,strfind(fileExtensions,'xml')));
%                    if isempty(hdf5files) %There arent h5files
%                        warning(['WARNING NO H5 FILES FOUND FOR ',folder{1}])
%                        destinationPath = [obj.rootFolder,filesep,folder{1},'_NO_FILES'];
%                        mkdir(destinationPath)
%                        statusErrors = statusErrors + 1;
%                        continue;
%                    elseif length(logFiles) > length(hdf5files)
%                        warning(['WARNING THERE MAY BE UNCOMPRESSED RAW FILES FOR ',folder{1}])
%                        destinationPath = [obj.rootFolder,filesep,folder{1},'_MISSING_FILES'];
%                        mkdir(destinationPath)
%                        statusErrors = statusErrors + 1;
%                    else
%                        destinationPath = [obj.rootFolder,filesep,folder{1}];
%                        mkdir(destinationPath)
%                    end
%                    
%                    downloadFiles = cell(size(hdf5files));
%                    i = 1;
%                    for file = {filesInFolder{hdf5files}}
%                        downloadFiles{i} = ['''',file{1},''''];
%                        i = i + 1;
%                    end
%                    obj.sftp.download(downloadFiles,destinationPath);
               end
           end
        end
        
        function getLocalConcatFiles(obj)
            allItems = dir(obj.rootFolder);
            allItems = allItems(3:end); %remove WINDOWS metadata folders
            if isempty(allItems)
                obj.downloadH5fromServer;
                allItems = dir(obj.rootFolder);
                allItems = allItems(3:end); %remove WINDOWS metadata folders
            end
            allFolders = allItems([allItems.isdir]);
            obj.folderStruct = [];
            i = 1;
            for dayFolder = {allFolders.name}
                regexpOut = regexp(dayFolder{1}, '(\d+)','tokens');
                day = regexpOut{end}{1};%Assumess date is the last sequence of numbers in the folder name
                fullPath = fullfile(obj.rootFolder, dayFolder{1});    
                
                %Adjust the folder's name
                correctFolderName = [obj.rootFolder,filesep,'Mouse-',obj.mouse,'-',day,'-',obj.experiment];
                if strcmp(correctFolderName,fullPath) == 0
                    movefile(fullPath,correctFolderName)
                    fullPath = correctFolderName;
                end
                
                obj.folderStruct(i).path = fullPath;
                obj.folderStruct(i).folderName = dayFolder{1};
                obj.folderStruct(i).date = datevec(day,obj.dateParser);
                obj.folderStruct(i).mouse = obj.mouse;
                
                allInsideitems = dir(fullPath);
                allInsideitems = allInsideitems(3:end); %remove WINDOWS metadata folders
                allFiles = allInsideitems(~[allInsideitems.isdir]);
                concatIdx = find(obj.checkRegexp({allFiles.name},obj.concatRegexp));
                if isempty(concatIdx)
                    warning(['cannot find concat files for ' dayFolder{1}])
                else
                    concatFiles = {allFiles(concatIdx).name};
                    k = 1;
                    for concatFile = concatFiles
                        noExtensionName = split(concatFile{1},'.');
                        regexpOut = regexp(noExtensionName{1}, '(\d+)','tokens');
                        session = regexpOut{end}{1};
                        obj.folderStruct(i).sessions(k).date = datevec([day,session],[obj.dateParser,obj.timeParser]);
                        obj.folderStruct(i).sessions(k).fileName = concatFile{1};
                        k = k + 1;
                    end
                    %Sort sessions in descending time order
                    [~,sortedIdx] = sort(cellfun(@datenum,{obj.folderStruct(i).sessions.date}));
                    obj.folderStruct(i).sessions = obj.folderStruct(i).sessions(sortedIdx);
                    i = i+1; %If no concat files were found, ignore the folder
                end
            end
        end
        
        function checkProcessedFiles(obj)
            
            folders = obj.chooseDays;
            
            for folder = folders   
                try
                    disp(repmat('#',[3,50]))
                    disp(['Checking: ', folder{1}.folderName])

                    allFiles = fileSolution.getFilesInFolder(folder{1}.path);

                    %Check if concat < dfof downsample files
                    disp(repmat('%',[2,50]))
                    disp(['Checking for dfof file(s) in ',folder{1}.folderName])
                    dfofFiles = {allFiles(find(obj.checkRegexp({allFiles.name},obj.dfofRegexp))).name};
                    numSessionsInDfof = sum(cellfun(@length,strfind(dfofFiles,'&'))+1); %Sessions are separated by '&'
                    if numSessionsInDfof < length(folder{1}.sessions) %this will process all sessions again (TODO?: process only not found sessions)
                        disp(['Missing dfof file(s) in ',folder{1}.folderName])
                        if obj.processIfMissing
                            disp(['Processing ...',folder{1}.folderName])
                            obj.preprocessConcat(folder{1})
                            allFiles = fileSolution.getFilesInFolder(folder{1}.path);
                        end
                    end

                    %Check if dfof > downsample
                    %TODO

                    %Check if dfof > analysis file
                    disp(repmat('%',[2,50]))
                    disp(['Checking for analysis file(s) in ',folder{1}.folderName])
                    if sum(obj.checkRegexp({allFiles.name},obj.analysisFileRegexp)) < sum(obj.checkRegexp({allFiles.name},obj.downsampleRegexp))
                        disp(['Missing analysis file(s) in ',folder{1}.folderName])
                        if obj.extractIfMissing
                            disp(['Extracting ...',folder{1}.folderName])
                            for k = find(obj.checkRegexp({allFiles.name},obj.downsampleRegexp))%this will extract all sessions again (TODO?: process only not found sessions)
                                fileForExtraction = fullfile(folder{1}.path,allFiles(k).name);
                                switch obj.analysisType
                                    case 'cellmax'
                                       fileSolution.cellMaxExtraction(fileForExtraction,obj.experiment);
                                    case 'pcaica'
                                       fileSolution.pcaicaExtraction(fileForExtraction,obj.experiment);
                                    otherwise
                                       error('extraction method no valid')
                                end
                            end
                            allFiles = fileSolution.getFilesInFolder(folder{1}.path);
                        end
                    end

                    %Check if there are decisions files
                    disp(repmat('%',[2,50]))
                    disp(['Checking for decision file in ',folder{1}.folderName])
                    if sum(obj.checkRegexp({allFiles.name},obj.analysisFileRegexp)) > sum(obj.checkRegexp({allFiles.name},obj.decisionsRegexp))
                        disp(['Missing decision file in ',folder{1}.folderName])
                        if obj.decideIfMissing %Experimental
                            disp(['Deciding ...',folder{1}.folderName])
                            for k = find(obj.checkRegexp({allFiles.name},obj.analysisFileRegexp))%this will process all sessions again (TODO?: process only not found sessions)
                                load([folder{1}.path,filesep,allFiles(k).name])
                                if ~isempty(find(obj.checkRegexp({allFiles(k).name},obj.analysisFileRegexp(1)))) %em
                                   emName = allFiles(k).name;
                                   nameParts = strsplit(emName,'-');
                                   template_name = char(strjoin(nameParts(1:end-1),'-'));
                                   fileName = allFiles(find(obj.checkRegexp({allFiles(:).name},{[template_name,'-dfof-downsampled.h5']}))).name;
                                   fileSolution.decideCells([folder{1}.path,filesep,fileName],emAnalysisOutput,'em')   
                                   %CHANGE THIS!
                                   %[~,fileName,fileExt] = fileparts(emAnalysisOutput.CELLMaxoptions.movieFilename);
                                   %moviePath = fullfile(folder{1}.path,[fileName,fileExt]);
                                   %fileSolution.decideCells(moviePath,emAnalysisOutput,'em')                                     
                                elseif ~isempty(find(obj.checkRegexp({allFiles(k).name},obj.analysisFileRegexp(2)))) %ica
                                   [~,fileName,fileExt] = fileparts(pcaicaAnalysisOutput.movieFilename);
                                   moviePath = fullfile(folder{1}.path,[fileName,fileExt]);
                                   fileSolution.decideCells(moviePath,pcaicaAnalysisOutput,'ica')
                                end
                                
                            end
                            allFiles = fileSolution.getFilesInFolder(folder{1}.path);
                        end
                    end

                    %Check if there is a extraction log, download if not
                    disp(repmat('%',[2,50]))
                    disp(['Checking log files in ',folder{1}.folderName])
                    if sum(obj.checkRegexp({allFiles.name},obj.logRegexp)) < sum(obj.checkRegexp({allFiles.name},obj.concatRegexp))
                        disp(['Missing log file in ',folder{1}.folderName])
                        try 
                            obj.fetchLogs(folder{1})
                            allFiles = fileSolution.getFilesInFolder(folder{1}.path);
                        catch ME
                            disp(ME.message)
                        end

                        %Rename logs
                        regFiles = {allFiles(find(obj.checkRegexp({allFiles.name},obj.logRegexp))).name};

                        for sessions = folder{1}.sessions
                            sessionTime = datestr(sessions.date,obj.timeParser);
                            sessionDate = datestr(sessions.date,obj.dateParser);
                            for regFile = regFiles
                                ext = split(regFile{1},'.');
                                if ~isempty(regexp(ext{1},sessionTime, 'once'))
                                    newName = fullfile(folder{1}.path,['Mouse-',folder{1}.mouse,'-',sessionDate,'_',sessionTime,'-',obj.experiment,'-log.',ext{end}]);
                                    movefile(fullfile(folder{1}.path,regFile{1}),newName)
                                end
                            end
                        end
                        allFiles = fileSolution.getFilesInFolder(folder{1}.path);                    

                    end

                    %Check if there is a behavior file, download if not
                    disp(repmat('%',[2,50]))
                    disp(['Checking for behavior files in ',folder{1}.folderName])
                    if 0%sum(obj.checkRegexp({allFiles.name},obj.behaviorRegexp)) < sum(obj.checkRegexp({allFiles.name},obj.concatRegexp))
                        disp(['Missing behavior file in ',folder{1}.folderName])
                        try 
                           obj.fetchBehavior(folder{1})
                           allFiles = fileSolution.getFilesInFolder(folder{1}.path);
                        catch ME
                            disp(ME.message)
                        end  

                        %ASSIGN EACH VIDEO OR VIDEOS TO THE CONCAT FILES\
                        for i = find(obj.checkRegexp({allFiles.name},obj.concatRegexp))
                                calcInf = hdf5info([folder{1}.path,filesep,allFiles(i).name]);
                                framesCalcN = calcInf.GroupHierarchy.Datasets.Dims(3);
                                calcTime = regexp(calcInf.Filename,'(\d+_\d+).*.h5$','tokens');%regexp(calcInf.Filename,'(\d+_\d+).h5$','tokens');
                                calcLog = allFiles(find(obj.checkRegexp({allFiles.name},{[calcTime{1}{1},'.*.txt$'],[calcTime{1}{1},'.*.xml$']}))).name;
                                missingFramesN = length(readFramesFromLog([folder{1}.path,filesep,calcLog]));
                                calcLen = framesCalcN + missingFramesN;
                                flag = [];
                                for j = find(obj.checkRegexp({allFiles.name},obj.behaviorRegexp))
                                    videoName = [folder{1}.path,filesep,allFiles(j).name];
                                    mov = VideoReader(videoName);
                                    if mov.NumberOfFrames == calcLen
                                        flag = [flag,j];
                                    end
                                end
                                if length(flag) > 1
                                    error(['Behavior videos ',allFiles(flag).name,' have the same length as calcium, ',allFiles(i).name])
                                elseif isempty(flag)
                                    error(['No behavior videos with the length of calcium, ',allFiles(i).name])
                                else
                                    videoName = [folder{1}.path,filesep,allFiles(flag).name];
                                    newName = fullfile(folder{1}.path,['Mouse-',folder{1}.mouse,'-',calcTime{1}{1},'-',obj.experiment,'-behavior.avi']);
%                                     try
                                        copyfile(videoName,newName)
%                                     catch ME
%                                         disp(ME.message)
%                                     end                      
                                end
                        end  
                        cellfun(@(x) delete([folder{1}.path,filesep,x]),{allFiles(find(obj.checkRegexp({allFiles.name},obj.behaviorRegexp))).name})
                        allFiles = fileSolution.getFilesInFolder(folder{1}.path);        
                    end
                        %CHECK THAT THERE IS A BEHAVIOR FILE FOR EVERY CONCAT
%                         if sum(obj.checkRegexp({allFiles.name},obj.behaviorRegexp)) ~= sum(obj.checkRegexp({allFiles.name},obj.concatRegexp))
%                             warning(['Different number of behavior movies and concat files for folder ',folder{1}.folderName,' revise'])
%                             for i = find(obj.checkRegexp({allFiles.name},obj.dfofRegexp))
%                                 calc = fileSolution.readHdf5([folder{1}.path,filesep,allFiles(j).name]);
%                                 for j = find(obj.checkRegexp({allFiles.name},obj.behaviorRegexp))
%                                     mov = VideoReader([folder{1}.path,filesep,allFiles(j).name]);
%                                 end
%                             end

%                         end

                        %ASSIGN BEHAVIOR TO SESSIONS

%                         %Get names
%                         behavFiles = {allFiles(find(obj.checkRegexp({allFiles.name},obj.behaviorRegexp))).name};
% 
%                         %Extract times
%                         behavFilesSimple = cellfun(@fileSolution.simplifyName,behavFiles,'UniformOutput',0);
%                         deleteDateStr = datestr(folder{1}.date,'yymmdd');
%                         behavFilesSimpleNoDate = cellfun(@(x) erase(x,deleteDateStr),behavFilesSimple,'UniformOutput',0);
%                         timesBehavCell = regexp(behavFilesSimpleNoDate,'(\d+)','tokens');
%                         timesBehav = cellfun(@(x) datevec(x{1}{1},obj.timeParser),timesBehavCell,'UniformOutput',0);   
% 
%                         %Sort behavior times
%                         [~,sortedIdx] = sort(cellfun(@datenum,timesBehav));
%                         behavFiles = behavFiles(sortedIdx);
% 
%                         %Assign times
%                         k = 1;
%                         for behavFile = behavFiles
%                             sessionTime = datestr(folder{1}.sessions(k).date,obj.timeParser);
%                             sessionDate = datestr(folder{1}.date,obj.dateParser);
%                             newName = fullfile(folder{1}.path,['Mouse-',folder{1}.mouse,'-',sessionDate,'_',sessionTime,'-',obj.experiment,'-behavior.avi']);
%                             movefile(fullfile(folder{1}.path,behavFile{1}),newName)
%                             k = k + 1;
%                         end
%                         allFiles = fileSolution.getFilesInFolder(folder{1}.path);
%                     end

                    %Check for traceEvent files
                    disp(repmat('%',[2,50]))
                    disp(['Checking for traceEvent file in ',folder{1}.folderName])
                    if sum(obj.checkRegexp({allFiles.name},obj.analysisFileRegexp)) > sum(obj.checkRegexp({allFiles.name},obj.tracesEventsRegexp))
                        disp(['Missing traceEvent file in ',folder{1}.folderName])
                        if obj.postProcessIfMissing
                            disp(['Post-processing ...',folder{1}.folderName])
                            for k = find(obj.checkRegexp({allFiles.name},obj.analysisFileRegexp))%this will process all sessions again (TODO?: process only not found sessions)
                                 %Load extraction file
                                 extractionFileName = fullfile(folder{1}.path,allFiles(k).name);
                                 [~,fileName,~] = fileparts(extractionFileName);
                                 extractionFile = load(extractionFileName);
                                 dummyField = fields(extractionFile);
                                 extractionFile = extractionFile.(dummyField{1});

                                 %Load decisons
                                 decisionsFileRegexp = cellfun(@(x) join({fileName,x},''),obj.decisionsRegexp,'UniformOutput',0);

                                 decisions = load(fullfile(folder{1}.path,allFiles(find(obj.checkRegexp({allFiles.name},decisionsFileRegexp))).name));
                                 dummyField = fields(decisions);
                                 decisions = decisions.(dummyField{1});

                                 %Get calcium name
                                 emName = allFiles(k).name;
                                 fileName = allFiles(find(obj.checkRegexp({allFiles(:).name},{[emName(10:end-25),'.+','-dfof-downsampled.h5']}))).name;
                                 fileName = fileName(1:end-3);
                                 fileExt = '.h5';
                                 %CHANGE THIS!
%                                  try
%                                     [~,fileName,fileExt] = fileparts(extractionFile.CELLMaxoptions.movieFilename); %em
%                                  catch
%                                     [~,fileName,fileExt] = fileparts(extractionFile.movieFilename); %ica
%                                  end
                                
                                 fileNameDfof = erase(fileName,'-downsampled'); %TODO?: Maybe search over files to find the downsampled with the obj.downsampleRegexp?
                                 calcium = fullfile(folder{1}.path,[fileNameDfof,fileExt]);

                                 %Get times in filename
                                 numbersInFileName = regexp(fileName,'(\d*)','tokens');
                                 timesInFilename = numbersInFileName(3:end);%First is mouse, second is session

                                 %Get log Files
                                 logFiles = {};
                                 for time = timesInFilename
                                    logRegexp = cellfun(@(x) join({time{1}{1},x},'.*'),obj.logRegexp,'UniformOutput',0);
%                                     logRegexp = logRegexp(1); %If both the txt and xml are present take only one of these
                                    logFiles = [logFiles,fullfile(folder{1}.path,allFiles(find(obj.checkRegexp({allFiles.name},logRegexp))).name)];
                                    %logFiles = logFiles(1);
                                 end

                                 %Get behav Files
                                 behavFiles = {};
                                 for time = timesInFilename
                                    behavRegexp = cellfun(@(x) join({time{1}{1},x},'.*'),obj.behaviorRegexp,'UniformOutput',0);
                                    behavFiles = [behavFiles,fullfile(folder{1}.path,allFiles(find(obj.checkRegexp({allFiles.name},behavRegexp))).name)];
                                 end

                                 %Get traces and events
                                 try
                                    fileSolution.getTracesEventsTraj(extractionFile,decisions,calcium,logFiles,behavFiles,obj.experiment)
                                 catch ME
                                     disp(ME.message)
                                     warning('ERROR GENERATING TRACES EVENTS')
                                     pause;
                                     disp('eh')
                                 end
                            end
                            allFiles = fileSolution.getFilesInFolder(folder{1}.path);
                        end
                    end

                catch ME
                    disp(ME.message)
                    fileID = fopen([obj.rootFolder,filesep,'errors.txt'],'a+');
                    fprintf(fileID,'%s\n',folder{1}.path);
                    fprintf(fileID,'%s\n',ME.message);
                    fprintf(fileID,'\n');
                    fclose(fileID);
                end 
            end
        end
        
        
        function fetchLogs(obj,folder)
           disp(['Downloading ',folder.folderName,' logs from the server'])
           obj.sftp.cd(obj.archiveRawCaPath)
            
           animalFolders = obj.sftp.commandToVar('ls');
           correctAnimalFolder = find(obj.checkRegexp(animalFolders',{obj.mouse}));
           if isempty(correctAnimalFolder)
               error('Could not find the logs, animal missing or duplicated')
           elseif length(correctAnimalFolder) > 1
               disp('Warning,duplicate folders:')
               disp({animalFolders{correctAnimalFolder}})
               disp('Searching for files in both folders:')
           end
           
           currentPath = obj.sftp.commandToVar('pwd');
           for animalFolder = {animalFolders{correctAnimalFolder}}
               obj.sftp.cd(currentPath{1})
               obj.sftp.cd(animalFolder{1})
           
               sessionFolders = obj.sftp.commandToVar('ls')';
               date = datestr(folder.date,obj.dateParser);
               dateVariations = {date,date(3:end)};
               correctSessionFolder = find(obj.checkRegexp(sessionFolders,dateVariations));
               if length(correctSessionFolder) < 1
                   errorMsg = 'Could not find the logs, day missing';
                   continue;
               else
                   obj.sftp.cd(sessionFolders{correctSessionFolder})
               end

               sessionFiles = obj.sftp.commandToVar('ls')';
               dayPlusSessionRegexp = cellfun(@(x) join({date,'.*',x},''),cellfun(@(x) datestr(x,obj.timeParser),{folder.sessions(:).date},'UniformOutput',false),'UniformOutput',false);

               sessionLogRegexp = {};
               for dayPlusSession = dayPlusSessionRegexp
                   sessionLogRegexp = [sessionLogRegexp,cellfun(@(x) join({dayPlusSession{1}{1},'.*',x},''),obj.logRegexp,'UniformOutput',false)];
               end
               correctsessionFiles = find(obj.checkRegexp(sessionFiles,sessionLogRegexp));
               if length(correctsessionFiles) < 1
                   errorMsg = 'Could not find the logs, file missing';
                   continue;
               else
                   obj.sftp.download({sessionFiles{correctsessionFiles}},folder.path)
                   disp('sucess')
                   return;
               end            
           end
           error(errorMsg)
        end
        
        function fetchBehavior(obj,folder)
           disp(['Downloading ',folder.folderName,' behavior from the server'])
           obj.sftp.cd(obj.archiveBehaviorPath)

           correctExperimentFolder = [];
           for experimentFolder = obj.sftp.commandToVar('ls')'
               if strcmpi(obj.simplifyName(obj.experiment),obj.simplifyName(experimentFolder{1}))
                   correctExperimentFolder = experimentFolder{1};
                   break;
               end
           end

           if isempty(correctExperimentFolder)
               error('Could not find behavior, experiment missing')
           else
               obj.sftp.cd(correctExperimentFolder)
           end

           animalFolders = obj.sftp.commandToVar('ls');
           correctAnimalFolder = find(obj.checkRegexp(animalFolders',{obj.mouse}));
           if length(correctAnimalFolder) ~= 1
               error('Could not find behavior, animal missing or duplicated')
           else
                obj.sftp.cd(animalFolders{correctAnimalFolder})
           end

           obj.sftp.cd('Behavioral-videos') %Might not be the same for each animal
           behavioralVideos = obj.sftp.commandToVar('ls')';
           behavioralVideosSimplified = cellfun(@(x) obj.simplifyName(x),behavioralVideos,'UniformOutput',0);

           date = datestr(folder.date,obj.dateParser);
           dateVariations = {date,date(3:end)};
           correctBehaviorFile = find(obj.checkRegexp(behavioralVideosSimplified,dateVariations));
           if length(correctBehaviorFile) < 1
               error('Could not find behavior, day missing')
           else
               downloadVideos = cell(size(correctBehaviorFile));
               i = 1;
               for video = {behavioralVideos{correctBehaviorFile}}
                   downloadVideos{i} = ['''',video{1},''''];
                   i = i + 1;
               end
               obj.sftp.download(downloadVideos,folder.path);
               disp('sucess')
           end                                 
        end
        
        function preprocessConcat(obj,folder)
            joinedSessionsIdx = ones([1,length(folder.sessions)]);%Sessions with the same number will be joined
            for k = 1:(length(folder.sessions)-1)
                if etime(folder.sessions(k+1).date,folder.sessions(k).date)/60 > obj.joinFileTreshold
                    joinedSessionsIdx(k+1:end) = joinedSessionsIdx(k+1:end)+1;
                end
            end
            for k = unique(joinedSessionsIdx) %Process the joined sessions
                disp(repmat('-',[1,50]))
                fileNums = joinedSessionsIdx == k;
                fileNames = {folder.sessions(fileNums).fileName};
                disp(fileNames)
                fileFull = cellfun(@(x) [folder.path,filesep,x],fileNames,'UniformOutput',false);
                parsedSessionDate = cellfun(@(x) datestr(x,obj.timeParser),{folder.sessions(fileNums).date},'UniformOutput',false);

                saveName = ['Mouse-',folder.mouse,'-',datestr(folder.date,obj.dateParser),'_',char(join(parsedSessionDate,'&')),'-',obj.experiment];
                [movieDS,movieDFOF,~,~] = fileSolution.preprocessMovie(fileFull);
                obj.createHDF5([folder.path,filesep,saveName,'-dfof.h5'],movieDFOF)
                obj.createHDF5([folder.path,filesep,saveName,'-dfof-downsampled.h5'],movieDS)
            end
        end
        
        function days = chooseDays(obj,type)
            if nargin < 2
                type = 'multiple';
            end
            days = {};
            dayList = cellfun(@(x) datestr(x,obj.dateParser),{obj.folderStruct.date},'UniformOutput',0);
            s = listdlg('PromptString','Select day(s):',...
                'SelectionMode',type,...
                'ListString',dayList);
            for si = s
                if length(obj.folderStruct(si).sessions) >= 1
                    days = [days,obj.folderStruct(si)];
                else
                    warning('The selected day has no avalaible sessions.')
                end
            end
            if strcmp(type,'single')
                days = days{1};
            end
        end
        
        function preprocessAllConcat(obj)
            for folder = obj.folderStruct
                disp(repmat('#',[1,50]))
                disp(['Processing: ', folder.folderName])
                obj.preprocessConcat(folder);
            end
        end
    end
    
    methods(Static)
              
        function [movieDS,movieDFOF,movieRegFilt,movieReg] = preprocessMovie(movie,varargin)
            movie = fileSolution.readHdf5(movie);
            
            freqLow = fileSolution.parseInput({varargin,'freqLow',1});
            freqHigh = fileSolution.parseInput({varargin,'freqHigh',4});
            filterReg = fileSolution.parseInput({varargin,'filterReg','divideByLowpass'});
%             filterDFOF = fileSolution.parseInput({varargin,'filterDFOF','kalmanTime'});
            dsTimeFactor = fileSolution.parseInput({varargin,'dsTimeFactor',4});
            
            movieReg = double(fileSolution.registerMovie(movie));
            cropRegMovie = double(fileSolution.cropRegisteredMovie(movieReg));
            
            movieRegFilt = double(fileSolution.filterMovie(cropRegMovie,filterReg,freqLow,freqHigh));
            movieDFOF = fileSolution.dfofMovie(movieRegFilt);
%             movieDFOFFilt = fileSolution.filterMovie(movieDFOF,filterDFOF);
            movieDS = fileSolution.downsampleTime(movieDFOF,dsTimeFactor);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%PROCESSING FUNCTIONS%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function cropRegMovie = cropRegisteredMovie(movie)
           %For each pixel, look at wheather in time it has been 0 at least once
            cropMask = ~all(movie,3);
            
            %Ensure we are not removing regions unconected to the borders
            rProps = regionprops(cropMask,'PixelList');
            for k = 1:length(rProps)
                pixelList = rProps(k).PixelList;
                if isempty(intersect(pixelList,[1,size(cropMask,1),size(cropMask,2)]))
                    cropMask(pixelList(:,1),pixelList(:,2)) = 0;
                end
            end
            
            %If there is still something to crop, crop it
            if any(cropMask(:))
                %Find largest rectangle
                bounds = FindLargestRectangles(~cropMask);
                cropRegMovie = movie(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),:);
            else
                cropRegMovie = movie;
            end
        end
        
        function movieREG = registerMovie(movie)
            disp('Registring...')
            movieREG = register1p(movie); %Eliminate first 10 samples?
        end
            
        function movieF = filterMovie(movie, filter, freqLow, freqHigh)
            disp('Filtering...')
            switch filter
                case 'divideByLowpass'
                    movieF = normalizeMovie(single(movie),'normalizationType','lowpassFFTDivisive','freqLow',freqLow,'freqHigh',freqHigh,'waitbarOn',0,'bandpassMask','gaussian');
                case 'bandpass'
                    movieF = normalizeMovie(single(movie),'normalizationType','fft','freqLow',freqLow,'freqHigh',freqHigh,'bandpassType','bandpass','showImages',0,'bandpassMask','gaussian');
                case 'myHighpass'
                    hLarge = fspecial('average', 40);
                    hSmall = fspecial('average', 2); 
                    movieF = zeros(size(movie));
                    parfor t = 1:size(movie,3)
                        movieF(:,:,t) = filter2(hSmall,movie(:,:,t)) - filter2(hLarge, movie(:,:,t));
                    end
                case 'kalmanTime'
                    movieF = Kalman_Stack_Filter(movie,.8);
                case 'wiener'
                    movieF = zeros(size(movie));
%                     H = fspecial('log');
                    parfor k = 1:size(movie,3)
                        movieF(:,:,k) = wiener2(movie(:,:,k),[5 5]);%imfilter(movie(:,:,k),H,'replicate');
                    end
                otherwise
                    warning(['Filter ',filter, 'not implemented'])
            end
        end

        function movieDF = dfofMovie(movie)
            disp('DFOFing...')
            movieDF = movie./mean(movie,3) - 1;
        end
        
        function movieDS = downsampleTime(movie,factor)
            disp('Downsampling in time...')
            downX = size(movie,1);
            downY = size(movie,2);
            downZ = floor(size(movie,3)/factor);
            for frame=1:downY
               downsampledFrame = imresize(squeeze(movie(:,frame,:)),[downX downZ],'bilinear');		  
               movie(1:downX,frame,1:downZ) = downsampledFrame;
            end
            movie(:,:,(downZ+1):end) = 0;
            thisMovieTmp = movie(:,:,1:downZ);
            movieDS = thisMovieTmp;        
        end
        
        function movieDS = downsampleSpace(movie,factor)
            disp('Downsampling in space...')
            downX = floor(size(movie,1)/factor);
            downY = floor(size(movie,2)/factor);
            downZ = size(movie,3);
            for frame=1:downZ
               downsampledFrame = imresize(squeeze(movie(:,:,frame)),[downX downY],'bilinear');
               movie(1:downX,1:downY,frame) = downsampledFrame;
            end
            movieDS = movie(1:downX,1:downY,:);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%UTILITY FUNCTIONS%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        function allFiles = getFilesInFolder(path)
            allItems = dir(path);
            allItems = allItems(3:end); %remove WINDOWS metadata folders
            allFiles = allItems(~[allItems.isdir]);
        end
        
        function movie = readHdf5(movie,varargin)
            if isstr(movie)
                hinf = hdf5info(movie);
                movie = hdf5read(hinf.GroupHierarchy.Datasets);
            end
            if iscell(movie)
                fullMovie = [];
                if isstr(movie{1}) 
                    %Check out if dimensions are compatible
                    movieSizes = zeros([length(movie),3]);
                    for k = 1:length(movie)
                        hinf = hdf5info(movie{k});
                        movieSizes(k,:) = hinf.GroupHierarchy.Datasets.Dims;
                    end
                    if length(unique(movieSizes(:,1))) > 1 || length(unique(movieSizes(:,2))) > 1
                        warning('MOVIE DIMENSIONS ARE NOT COMPATIBLE')
                        %Quick and dirty fix
                        dimX = min(movieSizes(:,1));
                        dimY = min(movieSizes(:,2));
                        offsetX = [floor((movieSizes(:,1)-dimX)/2),floor((movieSizes(:,1)-dimX)/2)+mod((movieSizes(:,1)-dimX),2)];
                        offsetY = [floor((movieSizes(:,2)-dimY)/2),floor((movieSizes(:,2)-dimY)/2)+mod((movieSizes(:,2)-dimY),2)];
                        for k = 1:length(movie)
                            hinf = hdf5info(movie{k});
                            sliceMovie = hdf5read(hinf.GroupHierarchy.Datasets);
                            sliceMovie = sliceMovie((1+offsetX(k,1)):(end-offsetX(k,2)),(1+offsetY(k,1)):(end-offsetY(k,2)),:);                           %sliceMovie = 
                            fullMovie = cat(3,fullMovie,sliceMovie);
                        end
                        
                        movie = fullMovie;
                        return;
                    end
                    
                    for moviePart = movie
                        hinf = hdf5info(moviePart{1});
                        sliceMovie = hdf5read(hinf.GroupHierarchy.Datasets);
                        fullMovie = cat(3,fullMovie,sliceMovie);
                    end
                else
                    for moviePart = movie
                        fullMovie = cat(3,fullMovie,moviePart{1});
                    end
                end
                movie = fullMovie;
            end
        end
       
        function simplifiedName = simplifyName(name)
            words = regexp(name,'([^\W_]*)','tokens');
            simplifiedName = '';
            for word = words
                simplifiedName = [simplifiedName,word{1}{1}];
            end
        end
        
        function match = checkRegexp(sourceText,regexpList)
            match = zeros(size(sourceText));
            k = 1;
            for text = sourceText
                regexpMatch = 0;
                for regexpSeq = regexpList
                    if regexp(text{1},regexpSeq{1})
                        regexpMatch = regexpMatch + 1;
                    end
                end
                if regexpMatch > 0
                    match(k) = 1;
                end
                k = k + 1;
            end
        end
        
        function playMovie(movie,varargin)
            
            deltaT = fileSolution.parseInput({varargin,'deltaT',1/20});
            movie = fileSolution.readHdf5(movie);
            
            upLim = quantile(movie(:),0.995);
            dwLim = quantile(movie(:),0.005);
            
            figure;
            for k = 1:size(movie,3)
                imagesc(movie(:,:,k),[dwLim,upLim]);
                pause(deltaT)
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function getTracesEventsTraj(extractionFile,decisions,calcium,missingFrames,behavior,experiment)
            tracesEvents = getTracesAndEvents(extractionFile, decisions, calcium, missingFrames,behavior);
            if ~isempty(tracesEvents)
                [pathstr,name,~] = fileparts(calcium);
                saveNamePart = regexp(name,'(.*\d+)','tokens');
                saveName = [pathstr,filesep,char(saveNamePart{1}),'-',experiment,'-TracesAndEvents.mat'];
                save(saveName,'tracesEvents')
            end
        end
        
        function pcaicaExtraction(file,experiment,varargin)
            
            [pathstr,name,~] = fileparts(file);
            
            nPC = fileSolution.parseInput({varargin,'nPC',750});
            nIC = fileSolution.parseInput({varargin,'nIC',500});
            inputDatasetName = fileSolution.parseInput({varargin,'inputDatasetName','/1'});
            pcaicaOutputUnits = fileSolution.parseInput({varargin,'pcaicaOutputUnits','fl'});
            
            [PcaOutputSpatial, PcaOutputTemporal, PcaOutputSingularValues, PcaInfo] = run_pca({file}, nPC, 'movie_dataset_name',inputDatasetName);
            
            calcInf = hdf5info(file);
            movieDims.x = calcInf.GroupHierarchy.Datasets.Dims(1);
            movieDims.y = calcInf.GroupHierarchy.Datasets.Dims(2);
            
            [IcaFilters, IcaTraces, IcaInfo] = run_ica(PcaOutputSpatial, PcaOutputTemporal, PcaOutputSingularValues, movieDims.x, movieDims.y, nIC, 'output_units',pcaicaOutputUnits);
            IcaTraces = permute(IcaTraces,[2 1]);
            
            pcaicaAnalysisOutput.IcaInfo = IcaInfo;
            pcaicaAnalysisOutput.PcaInfo = PcaInfo;
            pcaicaAnalysisOutput.IcaFilters = IcaFilters;
            pcaicaAnalysisOutput.IcaTraces = IcaTraces;
            pcaicaAnalysisOutput.nPCs = nPC;
            pcaicaAnalysisOutput.nICs = nIC;
            pcaicaAnalysisOutput.movieFilename = file;
            
            saveNamePart = regexp(name,'(.*\d+)','tokens');
            saveName = [pathstr,filesep,char(saveNamePart{1}),'-',experiment,'-pcaicaAnalysis.mat'];
            save(saveName,'pcaicaAnalysisOutput','-v7.3');
            
        end
        
        function cellMaxExtraction(file,experiment,varargin)
            
            [pathstr,name,~] = fileparts(file);
            
            emOptions.useParallel = fileSolution.parseInput({varargin,'useParallel',1});
            emOptions.movieDatasetName = fileSolution.parseInput({varargin,'movieDatasetName','/1'});
            
            CELLMaxoptions.initMethod = 'grid';
            CELLMaxoptions.gridSpacing = 9;
            CELLMaxoptions.gridWidth = 4;
            CELLMaxoptions.inputSizeManual = 0;
            CELLMaxoptions.subsampleMethod = 'resampleRemaining';
            CELLMaxoptions.percentFramesPerIteration = 0.5000;
            CELLMaxoptions.percentRemainingSubsample = 0.7500;
            CELLMaxoptions.maxSqSize = 101;
            CELLMaxoptions.threshForElim = 0.0050;
            CELLMaxoptions.localICimgs = [];
            CELLMaxoptions.localICtraces = [];
            CELLMaxoptions.minIters = 200;
            CELLMaxoptions.maxIters = 460;
            CELLMaxoptions.numSigmasThresh = 0.5000;
            CELLMaxoptions.nParallelWorkers = inf;
            CELLMaxoptions.generateNovelSeed = 1;
            CELLMaxoptions.numFramesRandom = 2000;
            CELLMaxoptions.readMovieChunks = 1;
            CELLMaxoptions.movieFilename = file;
            emOptions.CELLMaxoptions = fileSolution.parseInput({varargin,'CELLMaxoptions',CELLMaxoptions});
            
            startTime = tic;
            [emAnalysisOutput, ~] = CELLMax_Wrapper(file,'options',emOptions);
            
            emOptions.CELLMaxoptions.sqSizeX = [];
            emOptions.CELLMaxoptions.sqSizeY = [];
            emAnalysisOutput.dsCellTraces = emAnalysisOutput.cellTraces;
            emOptions.CELLMaxoptions.numSignalsDetected = size(emAnalysisOutput.dsCellTraces,1);
            emOptions.versionCellmax = emAnalysisOutput.versionCellmax;
            emOptions.time.startTime = startTime;
            emOptions.time.endTime = toc(startTime);
            
            saveNamePart = regexp(name,'(.*\d+)','tokens');
            saveName = [pathstr,filesep,char(saveNamePart{1}),'-',experiment,'-emAnalysis.mat'];
            save(saveName,'emAnalysisOutput','-v7.3','emOptions');
            rmdir(fullfile(pathstr,'tmpImages'),'s')
        end
        
        function decideCells(moviePath,analysisOutput,type)
            cInf = cellInfo(moviePath,permute(analysisOutput.cellImages,[3,1,2]),analysisOutput.scaledProbability,[]);
            cInf.setTresholds('showProgress',0,'tresholds',{'getOverlap','>0.5','getglobalSNR','>1','getCellShapePropieties.Eccentricity','<=0.9','getCellShapePropieties.Area','>100','getsScore','<0.02'});
            cInf.saveDecisions(type,'skipConfirmation',1)
        end
        
        function createHDF5(filePath,movie,varargin)
            datasetName = fileSolution.parseInput({varargin,'datasetName','/1'});
            try
                createHdf5File(filePath,datasetName,movie);
            catch
                disp('filename may be too large or directory may not exist, fix manually')
                pause;
            end
        end
    end   
    
    methods(Static, Access = private)

        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%UTILITY FUNCTIONS%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        function [varargout] = parseInput(varargin)
            %{{inputs},optional(arg1 default1, arg2 default2, ...)}
            tmpVar = varargin{1};
            if iscell(tmpVar{1})
                argumentList = tmpVar(2:2:end);
                defaultVal = tmpVar(3:2:end);
                input = tmpVar{1};
                inputArg = input(1:2:end);
                inputVal = input(2:2:end);
                for k = 1:length(argumentList)
                    idx = strcmp(argumentList{k},inputArg);
                    if sum(idx) == 0
                        varargout{k} = defaultVal{k};
                    else
                        varargout{k} = inputVal{idx};
                    end
                end                
            else
                for k = 2:2:length(tmpVar)
                    varargout{k} = tmpVar{k};
                end
            end
        end
    end
        
end
