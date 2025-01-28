classdef singleAnimalAnalysis < dataAnalysis
    %This is a parent class giving basic utility to navegate the folder
    %structure of single animal, along particular animal parameters. To be
    %inherited from more specific classes.
    
    properties
        
        %Folder structures
        rootFolder = '';
        dayStruct = '';
        mouse = '';
        folderName = '';
        
    end
    
    methods
        
        function obj = singleAnimalAnalysis(folder)
            
            if nargin < 1 || isempty(folder)
                obj.rootFolder = uigetdir('','Select the animal to analyze: ');
            else
                obj.rootFolder = folder;
            end
            
            nameSep = strsplit(char(obj.rootFolder),filesep);nameSep = nameSep(~cellfun(@isempty,nameSep));
            obj.folderName = nameSep{end};
            regexpOut = regexp(obj.folderName, 'ouse\D*(\d+)','tokens');
            try
                obj.mouse = regexpOut{1}{1};
            end
            obj.getAnalyzableSessions;
        end
        
        function sessions = chooseSessions(obj,sessionVar)
            %Enter sessions with cell e.g.{1,{1,2},3}means the
            %first and second sessions of day 1 and all the sessions of day
            %3 or use the gui
            sessions = [];
            if nargin < 2 || isempty(sessionVar)
                %If no sessions are specified bring up the gui
                %dayList = cellfun(@(x) datestr(x,obj.dateParser),{obj.dayStruct.date},'UniformOutput',0);
                dayList = {obj.dayStruct.folderName};
                s = listdlg('PromptString','Select day(s):',...
                    'SelectionMode','multiple',...
                    'ListString',dayList,'listSize',[460,500]);
                for si = s
                    if length(obj.dayStruct(si).sessions) == 1
                        sessions = [sessions,obj.dayStruct(si).sessions];
                    elseif length(obj.dayStruct(si).sessions) > 1
                        timeList = {};
                        for fileN = {obj.dayStruct(si).sessions(:).tracesEventsFileName}
                            regexpOut = regexp(fileN{1}, '(\d+&\d+)','tokens');
                            if isempty(regexpOut)
                                regexpOut = regexp(fileN{1}, '(\d+)','tokens');     
                            end
                            timeList = [timeList,regexpOut{end}{1}];
                        end

                        s2 = listdlg('PromptString',['Select session(s) for day: ',dayList(si)],...
                            'SelectionMode','multiple',...
                            'ListString',timeList);
                        for s2i = s2
                            sessions = [sessions,obj.dayStruct(si).sessions(s2i)];
                        end
                    else
                        warning('The selected day has no avalaible sessions. Run fileSolution to fix')
                    end
                end
            elseif isstruct(sessionVar)
                sessions = sessionVar;
            elseif isstr(sessionVar)
                if strcmp(sessionVar,'all')
                    for day = obj.dayStruct
                        for session = day.sessions
                            sessions = [sessions,session];
                        end
                    end
                else
                    error('Session variable not understood')
                end
            else
                if isscalar(sessionVar)
                    sessionVar = {sessionVar};
                end
                for k = 1:length(sessionVar)
                    if iscell(sessionVar{k})
                        sessIdx = [sessionVar{k}{:}]; %Take selected sessions
                        for i = sessIdx
                            sessions = [sessions,obj.dayStruct(sessionVar{k-1}).sessions(i)];
                        end
                    else
                        if k ~= length(sessionVar)
                            if iscell(sessionVar{k+1})
                                continue;
                            end
                        end
                        sessIdx = 1:length(obj.dayStruct(sessionVar{k}).sessions);
                        for i = sessIdx
                            sessions = [sessions,obj.dayStruct(sessionVar{k}).sessions(i)];
                        end
                    end
                end
            end
        end
        
        function [outField,nSessions,sessions,field] = getField(obj,varargin)
            sessions = obj.chooseSessions(dataAnalysis.parseInput({varargin,'sessions',[]}));
            message = dataAnalysis.parseInput({varargin,'message','Select output:'});
            nSessions = length(sessions);
            
            decideFieldForEachFile = dataAnalysis.parseInput({varargin,'decideForEachFile',0});
            
            [outField] = cell([1,length(sessions)]);
             for k = 1:length(sessions)
                %Get position and responses from output file           
                output = obj.loadOutput(sessions(k));
                
                %Get field
                if decideFieldForEachFile || ~exist('field','var')
                    field = dataAnalysis.parseInput({varargin,'field',@() obj.selectOutput(output,message)});
                    if iscell(field);field = field{k};end
                end
                outField{k} = eval(['output.',field]);
             end
        end
        
        function [x,r,nSessions,sessions,rField] = getXR(obj,varargin)
            %OPTIONAL INPUTS: session, rField, decideRForEachFile
            
            sessions = obj.chooseSessions(dataAnalysis.parseInput({varargin,'sessions',[]}));
            nSessions = length(sessions);
            
            decideRForEachFile = dataAnalysis.parseInput({varargin,'decideRForEachFile',0});
            
            [x,r] = deal(cell([1,length(sessions)]));
                        
            for k = 1:length(sessions)
                %Get position and responses from output file           
                output = obj.loadOutput(sessions(k));
                
                %Get response field
                if decideRForEachFile || ~exist('rField','var')
                    rField = dataAnalysis.parseInput({varargin,'rField',@() obj.selectOutput(output,'Select cell response')});
                    if iscell(rField);rField = rField{k};end
                end

                %Get traces and events and convert spikes into one-hot embedding   
                x{k} = output.position;
                r{k} = dataAnalysis.ind2OneHot(eval(['output.',rField]),length(x));

                %Check that length calcium == length behav
                if length(x{k}) ~= size(r{k},1)
                    warning(['Calcium traces have a length of ',num2str(size(r{k},1)),'. However, behavioral trajectories have a length of ',num2str(length(x{k}))])
                end
            end
        end
        
        function output = loadOutput(obj,sessions)
            if nargin < 2
                sessions = obj.chooseSessions;
            end
            output = cell([1,length(sessions)]);
            for k = 1:length(sessions)
%                 disp(['...loading ',sessions(k).tracesEventsFileName,' ...'])
                fullPath = fullfile(sessions(k).rootFolder,sessions(k).folderName,sessions(k).tracesEventsFileName);
                output{k} = load(fullPath);dummyField = fields(output{k});output{k} = output{k}.(dummyField{1});%load regardless of struct name 
                %If positions missing, get them from behav
                if ~isfield(output{k},'position')
                    singleAnimalAnalysis.appendTrajectoriesToOutput(sessions(k))
                    output{k} = load(fullPath);dummyField = fields(output{k});output{k} = output{k}.(dummyField{1});%load regardless of struct name 
                end
            end
            if length(sessions) == 1;output = output{1};end
        end

    end
    
    methods (Static)
        
        function scores = getClassifierConfidence(sessions,PCdecisions)           
           %Runs the decoder over the selected sessions to get the scores of the desired cells.
           % If PCdecisions is empty the decisions in the decision file are used. 
           
           scores = {};
           for k = 1:length(sessions)
                decisions = load(fullfile(sessions(k).rootFolder,sessions(k).folderName,sessions(k).decisionsFileName));
                emFile = load(fullfile(sessions(k).rootFolder,sessions(k).folderName,sessions(k).analysisFileFileName));
                [~,fileName,fileExt] = fileparts(emFile.emAnalysisOutput.movieFilename);
                moviePath = fullfile(sessions(k).rootFolder,sessions(k).folderName,[fileName,fileExt]);
                cInf = cellInfo(moviePath,permute(emFile.emAnalysisOutput.cellImages,[3,1,2]),emFile.emAnalysisOutput.scaledProbability,[]);
                sessionScores = cInf.predictClassifier;
                if nargin < 2 || isempty(PCdecisions)
                    scores{k} = sessionScores(logical(decisions.validCellMax));                    
                else
                    scores{k} = sessionScores(PCdecisions{k}(logical(decisions.validCellMax)));
                end
           end
       end
       
    end
    
    methods(Access = protected)
        
        function getAnalyzableSessions(obj)
            allItems = dir(obj.rootFolder);
            allItems = allItems(3:end); %remove WINDOWS metadata folders
            allFolders = allItems([allItems.isdir]);
            i = 1;
            for dayFolder = {allFolders.name}
                regexpOut = regexp(dayFolder{1}, '(\d+)','tokens');
                day = regexpOut{end}{1};%Assumess date is the last sequence of numbers in the folder name
                fullPath = fullfile(obj.rootFolder, dayFolder{1});
                
                
                allInsideitems = dir(fullPath);
                allInsideitems = allInsideitems(3:end); %remove WINDOWS metadata folders
                allFiles = allInsideitems(~[allInsideitems.isdir]);
                tracesEventsIdx = find(obj.checkRegexp({allFiles.name},obj.tracesEventsRegexp));
                if isempty(tracesEventsIdx)
                    warning(['cannot find traces/events files for ' dayFolder{1}])
                else
                    
                    obj.dayStruct(i).path = fullPath;
                    obj.dayStruct(i).folderName = dayFolder{1};
                    obj.dayStruct(i).date = datevec(day,obj.dateParser);
                    obj.dayStruct(i).mouse = obj.mouse;
                    
                    tracesEventsFiles = {allFiles(tracesEventsIdx).name};
                    k = 1;
                    for tracesEventsFile = tracesEventsFiles
                        noExtensionName = split(tracesEventsFile{1},'.');
                        regexpOut = regexp(noExtensionName{1}, '(\d+&\d+)','tokens');
                        if isempty(regexpOut)
                            regexpOut = regexp(noExtensionName{1}, '(\d+)','tokens');
                            session = regexpOut{3}{1};
                        else
                            session = regexpOut{end}{1};
                        end
                        sessionTimes = regexp(session, '(\d+)','tokens');
                        sessionTimes = cat(2,sessionTimes{:});
                        
                        formatTimes = cellfun(@(x) datevec([day,x],[obj.dateParser,obj.timeParser]),sessionTimes,'UniformOutput',0);
                        obj.dayStruct(i).sessions(k).date = cat(1,formatTimes{:});
                        obj.dayStruct(i).sessions(k).folderName = dayFolder{1};
                        obj.dayStruct(i).sessions(k).rootFolder = obj.rootFolder;
                        obj.dayStruct(i).sessions(k).tracesEventsFileName = tracesEventsFile{1};
                        
                        %If they exist, Find also cellmap and decision files
                        try
                            analysisFileRegexp = cellfun(@(x) [session,'.*',x],obj.analysisFileRegexp,'UniformOutput',0);
                            analysisFileIdx = find(obj.checkRegexp({allFiles.name},analysisFileRegexp));
                            obj.dayStruct(i).sessions(k).analysisFileFileName = allFiles(analysisFileIdx).name;

                            decisionsFileRegexp = cellfun(@(x) [session,'.*',x],obj.decisionsRegexp,'UniformOutput',0);
                            decisionsFileIdx = find(obj.checkRegexp({allFiles.name},decisionsFileRegexp));
                            obj.dayStruct(i).sessions(k).decisionsFileName = allFiles(decisionsFileIdx).name;
                        catch
                            
                        end
                        %
                        
                        k = k + 1;
                    end
                    i = i+1; %If no concat files were found, ignore the folder
                end
            end 
        end
        
    end
    
    methods(Static, Access = private)
        
        function appendTrajectoriesToOutput(session)
            %get behavior files
            fullPath = fullfile(session.rootFolder,session.folderName,session.tracesEventsFileName);
            allItems = dir(fullfile(session.rootFolder,session.folderName));
            allItems = allItems(3:end); %remove WINDOWS metadata folders
            allFiles = allItems(~[allItems.isdir]);
            behaviorFiles = cell([1,size(session.date,1)]);
            for k = 1:size(session.date,1)
                thisSess = session.date(k,:);
                day = datestr(thisSess,dataAnalysis.dateParser);
                time = datestr(thisSess,dataAnalysis.timeParser);
                behavFiles = {allFiles(find(session.checkRegexp({allFiles.name},{[day,'.*',time,'.*behavior']}))).name};
                behaviorFiles{k} = fullfile(session.rootFolder,session.folderName,behavFiles{1});
            end
            appendPosToTracesEventsFile(fullPath,behaviorFiles)
        end
        
    end
    
end

