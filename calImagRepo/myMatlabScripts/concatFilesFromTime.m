function concatFilesFromTime(folder,fileRegexp,joinFileTreshold)

    dateParser = 'yyyymmdd';
    timeParser = 'HHMMSS';
    
    folder = getFolderStruct(folder,fileRegexp);
    for d = 1:length(folder)
        disp(repmat('%',[3,50]))
        disp(['Day ',num2str(d),'/',num2str(length(folder))])
        joinedSessionsIdx = ones([1,length(folder(d).sessions)]);
        for k = 1:(length(folder(d).sessions)-1)
            if etime(folder(d).sessions(k+1).date,folder(d).sessions(k).date)/60 > joinFileTreshold
                joinedSessionsIdx(k+1:end) = joinedSessionsIdx(k+1:end)+1;
            end
        end

        for k = unique(joinedSessionsIdx) %Process the joined sessions
                disp(repmat('-',[1,50]))
                fileNums = joinedSessionsIdx == k;
                fileNames = {folder(d).sessions(fileNums).fileName};
                disp(fileNames)
                nameParts = strsplit(fileNames{1},'_');
                regexpOut = regexp(fileNames{1}, '.*\d+-(.*)$','tokens');
                nameEnd = regexpOut{1}{1};
                
                fileFull = cellfun(@(x) [folder(d).path,filesep,x],fileNames,'UniformOutput',false);
                
                if length(fileFull) > 1
                    parsedSessionDate = cellfun(@(x) datestr(x,timeParser),{folder(d).sessions(fileNums).date},'UniformOutput',false);
                    saveName = [folder(d).path,filesep,nameParts{1},'_',char(join(parsedSessionDate,'&')),'-',nameEnd];
                    try
                        fileSolution.createHDF5(saveName,fileSolution.readHdf5(fileFull))
                    catch
                        disp(['Unable to create file ',saveName,' skipping...']);
                        continue;
                    end
                    cellfun(@delete,fileFull)
                else
                    warning(['Only single file ',fileFull{1},' to join, skipping...']);
                end
        end
    end
end

