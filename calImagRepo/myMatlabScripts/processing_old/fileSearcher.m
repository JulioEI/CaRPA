
myDir = 'X:\Mouse2001\Mouse2001-20140116-icx';
%myDir = 'X:\Mouse2002\Mouse2002-20140116-icx';
%myDir = 'X:\Mouse2002\Mouse2002-20140117-icx';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allFiles = dir(myDir);
files = containers.Map;

%myRegexp = '^recording_\d*_(\d*)\.(\w*)$'; %This one exlcudes the parts of raw files ([...]-001.raw)
%myRegexp = '^recording_\d*_(\d*).*\.(\w*)$'; %This one does not exclude the parts of raw files
myRegexp = '\w*_\d*_(\d*).*\.(\w*)$'; %With snapshots and others

for k = 1:size(allFiles)
	name = regexp(allFiles(k).name,myRegexp,'tokens');
    if ~isempty(name)
        name = name{1};
        session = str2double(name{1});
        type = name{2};
        if ~isKey(files,type)
            files(type)= session;
        elseif ~ismember(session,files(type))
            files(type) = [files(type),session];
        end
    end
end

fprintf('\nIn the folder: %s\n',myDir)
for key = files.keys
    sessions = files(key{1});
    fprintf('\n%d %s:',size(sessions,2),key{1})
    fprintf('\n\t%d',sort(sessions))
    fprintf('\n')
end
