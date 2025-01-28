
myDir = 'E:\Processing\Mouse2001\'; %Last \ is important!!!
excludedSessions = {};%{'20140114'};%{'20140113','20140114'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileID = fopen('settings.prj','w');

globalLog = '(iUserDict\nUserDict\np0\n(dp1\nS''data''\np2\n(dp3\nS''log''\np4\n(lp5\n(cdatetime\ndatetime\np6\n';
xline = '(S''none''\n';

allSesions = dir(myDir);
i = 7;
for folder = flip({allSesions.name})
    sessionName = regexp(folder,'^Mouse.*-(\d*)-icx$','tokens');
    sessionName = sessionName{1};
    if ~isempty(sessionName) && ~sum(strcmp(sessionName{1},excludedSessions))
        sessionName = sessionName{1};
        disp(['Writting session ',sessionName{1}])
        folderDir = strcat(myDir,folder{1});
        allRecordings = dir(folderDir);
        
        %Getting all the raw files for each of the times.
        rawFiles = containers.Map;
        allFiles = {allRecordings.name};
        for file = allFiles
            singleNumber = regexp(file,'\w+\d+_(\d+)-*\d*.raw','tokens');
            singleNumber = singleNumber{1};
            if ~isempty(singleNumber)
                index = singleNumber{1}{1};
                if ~isKey(rawFiles,singleNumber{1})
                   rawFiles(index)= file;
                else
                   rawFiles(index)= [rawFiles(index);file];
                end
            end
        end
        
        %Searching for the last raw, getting txt data
        for key = flip(rawFiles.keys)
            
            if i ~= 7 %First one does not have a(g
                globalLog = [globalLog,'a(g6\n'];
            end
            
            globalLog = [globalLog,xline,'p',int2str(i),'\n'];
            i = i + 1;
            globalLog = [globalLog,'tp',int2str(i),'\n'];
            i = i + 1;
            globalLog = [globalLog,'Rp',int2str(i),'\n'];
            i = i + 1;      
            
            raws = rawFiles(key{1});
            if size(raws,1) > 1
                logRaw = raws(end-1);
            else
                logRaw = raws(1);
            end
            if strcmp(file,'')
            cropname = regexp(logRaw,'(recording_\d+_\d+)','tokens');
            cropname = cropname(1);
            if ~isempty(cropname{1})
                cropname = cropname{1};
                filePath = strcat(folderDir,'\',cropname{1},'.txt');
                try
                    text = fileread(filePath{1});
                    
                    frames = regexp(text,'FRAMES: (\d*)','tokens');
                    frames = frames(1);
                    time = regexp(text,'TIME: (\d*:\d*)','tokens');
                    time = time(1);
                    dropedFrames =  regexp(text,'DROPPED COUNT: (\d*)','tokens');
                    dropedFrames = dropedFrames(1);
                               
                catch
                    disp(['WARNING, .TXT NOT FOUND FOR FILE ',cropname{1}{1},'.txt'])
                    disp('THIS MAY CAUSE PROBLEMS LATER ON')
                    frames = {{'null'}};
                    time = {{'null'}};
                    dropedFrames = {{'null'}};
                end
                recordingInfo = ['S''Record ',frames{1}{1},' frames in ',time{1}{1},'. Dropped ',dropedFrames{1}{1},' frames.\n'];
                nameLog = ['S''Recorded (',logRaw{1},')''\n'];
                
                globalLog = [globalLog,recordingInfo];
                globalLog = [globalLog,'p',int2str(i),'\n'];
                i = i + 1;
                globalLog = [globalLog,'I00\n'];
                globalLog = [globalLog,'tp',int2str(i),'\n'];
                i = i + 1;
                globalLog = [globalLog,'a(g6\n',xline,'p',int2str(i),'\n'];
                i = i + 1;        
                globalLog = [globalLog,'tp',int2str(i),'\n'];
                i = i + 1;
                globalLog = [globalLog,'Rp',int2str(i),'\n'];
                i = i + 1;
            
            else
                snapshotName = raws{1};
                nameLog = ['S''Snapshot (',snapshotName(1:end-4),')''\n'];
            end
            
            globalLog = [globalLog,nameLog,'p',int2str(i),'\n'];
            i = i + 1;
            globalLog = [globalLog,'I00\n'];
            globalLog = [globalLog,'tp',int2str(i),'\n'];
            i = i + 1;            
        end
    end
end

globalLog = [globalLog,'asS''settings\n','p',int2str(i),'\n'];
i = i + 1;
globalLog = [globalLog,'(dp',int2str(i),'\nsS''counter''\n'];
i = i + 1;
globalLog = [globalLog,'p',int2str(i),'\nI224\nsS''created''\n'];
i = i + 1;
globalLog = [globalLog,'p',int2str(i),'S''',datestr(datetime),'''\n'];
i = i + 1;
globalLog = [globalLog,'p',int2str(i),'\nsS''rois''\n'];
globalLog = [globalLog,'p',int2str(i),'\n'];
i = i + 1;
globalLog = [globalLog,'(lp',int2str(i),'\n'];
i = i + 1;
globalLog = [globalLog,'(dp',int2str(i),'\n'];
i = i + 1;
globalLog = [globalLog,'(S''height''\n','p',int2str(i),'\nI100\n'];
i = i + 1;
globalLog = [globalLog,'(S''width''\n','p',int2str(i),'\nI100\n'];
i = i + 1;
globalLog = [globalLog,'(S''top''\n','p',int2str(i),'\nI0\n'];
i = i + 1;
globalLog = [globalLog,'(S''enable''\n','p',int2str(i),'\nI100\n'];
i = i + 1;
globalLog = [globalLog,'(S''left''\n','p',int2str(i),'\nI0\n'];
i = i + 1;
globalLog = [globalLog,'(sasS''record_format''\n','p',int2str(i),'\n'];
i = i + 1;
globalLog = [globalLog,'(lp',int2str(i),'\n'];
i = i + 1;
globalLog = [globalLog,'(S''File Type''\n','p',int2str(i),'\n'];
i = i + 1;
globalLog = [globalLog,'(S''recording''\n','p',int2str(i),'\n'];
i = i + 1;
globalLog = [globalLog,'tp',int2str(i),'\n'];
i = i + 1;
globalLog = [globalLog,'a(S''Separator''\n','p',int2str(i),'\nS''_''\n'];
i = i + 1;
globalLog = [globalLog,'p',int2str(i),'\n'];
i = i + 1;
globalLog = [globalLog,'tp',int2str(i),'\n'];
i = i + 1;
globalLog = [globalLog,'a(S''Date''\n','p',int2str(i),'\nS''YYYYMMDD''\n'];
i = i + 1;
globalLog = [globalLog,'p',int2str(i),'\n'];
i = i + 1;
globalLog = [globalLog,'tp',int2str(i),'\n'];
i = i + 1;
globalLog = [globalLog,'a(g',int2str(i),'\n'];
i = i + 1;
globalLog = [globalLog,'g',int2str(i),'\n'];
i = i + 1;
globalLog = [globalLog,'tp',int2str(i),'\n'];
i = i + 1;
globalLog = [globalLog,'a(S''Time''\n','p',int2str(i),'\nS''HHMMSS (24 hour)''\n'];
i = i + 1;
globalLog = [globalLog,'p',int2str(i),'\n'];
i = i + 1;
globalLog = [globalLog,'tp',int2str(i),'\n'];
i = i + 1;
globalLog = [globalLog,'asS''name''\n','p',int2str(i),'\nVMouse',myDir(end-4:end-1),'\n'];
i = i + 1;
globalLog = [globalLog,'p',int2str(i),'\nssb.'];


fprintf(fileID,globalLog);
fclose(fileID);

