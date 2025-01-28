% pathToMiji = 'C:\Users\csc\Desktop\Fiji.app\scripts';
addpath('C:\Users\csc\Desktop\NoRMCorre-master\');
% addpath(genpath('C:\Users\csc\Desktop\caImagg\miniscope_analysis-bahanonu-miniscopeAnalysisClass_18_08_17CUSTOM\'));

folders = {'E:\Processing\Mouse2029\'};%{'E:\Processing\Mouse2029','E:\Processing\Mouse2028 - 125itEM','E:\Processing\Mouse2028 - 200itEM','E:\Processing\Mouse2028 - 500itEM','E:\Processing\Mouse2028'};
%daysExcluded = {};%{'201502281-icx','201502282-icx','201503011-icx','201503012-icx'};%{'201502281-icx','201502282-icx','201503011-icx','201503012-icx'};
%daysIncluded = {'20150305-icx'};

%%
%Puts the files into dirs
for k = 1:length(folders)
    allFiles = dir(folders{k});
    allFiles = allFiles(3:end); %remove WINDOWS metadata folders
    if (sum([allFiles.isdir]))==0 %if the folder has no folders inside
        for file = {allFiles.name}
            if regexp(file{1},'concat')
                dayStrTmp = regexp(file{1},'concat_recording_(\d+)_\d+.h5','tokens');
                dayNameTmp = [folders{k},'\',dayStrTmp{1}{1},'-icx'];
                mkdir(dayNameTmp);
                sourceTmp = [folders{k},'\',file{1}];
                destinationTmp = [dayNameTmp,'\',file{1}];
                movefile(sourceTmp,destinationTmp);
            end
        end
    end
end

%%

concatFiles = {};
concatFolders = {};
for k = 1:length(folders)
    allFiles = dir(folders{k});
    timePerDays = [];
    predNPerDays = [];
    for dayFolder = {allFiles.name}
        if ~isempty(strfind(dayFolder{1},'-icx'))% &&  1==sum(strcmp(dayFolder{1},daysIncluded)) &&  1~=sum(strcmp(dayFolder{1},daysExcluded))%Avoids only specified days
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

%Reg
% lastName = concatFolders{1};
% repeatCt = 1;
for k = 1:length(concatFolders)
    disp('%%%%%%%%%%%%%%%%%%%%%')
    disp(['Processing ',concatFolders{k},'  (',num2str(k),' of ',num2str(length(concatFolders)),')'])
%     if ~strcmp(lastName,concatFolders{k})
%         lastName = concatFolders{k};
%         repeatCt = 1;
%     end
    Mpr = register1p(concatFiles{k});
    mkdir([concatFolders{k},'\originalConcat'])
    movefile(concatFiles{k},[concatFolders{k},'\originalConcat'])
    newName = concatFiles{k};
    newName = [newName(1:end-3),'_REG',newName(end-2:end)];
%     newName = [newName(1:end-3),'_REG_p',num2str(repeatCt),newName(end-2:end)];
    createHdf5File(newName,'img',Mpr);
    clear Mpr
    clear newName
%     repeatCt = repeatCt + 1;
end

% Start automation
% %%
% days = unique(concatFolders);
% for k = 1:length(days)
%     disp('%%%%%%%%%%%%%%%%%%%%%')
%     disp(['Processing ',days{k},'  (',num2str(k),' of ',num2str(length(days)),')'])
%     tic;
%     clear obj
%     obj = miniscopeAnalysis;
%     obj.modelAddNewFolders(days{k});
%     warning('off')
%     obj.modelPreprocessMovie('img',pathToMiji)
%     disp(i)
%     warning('on')
%     toc;
% end