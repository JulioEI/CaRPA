function checkCaBeh(folderName)

    if nargin < 1
        folderName = uigetdir('E:\Processing\');
    end
    
    fileID = fopen([folderName,filesep,'checkCaBehLog.txt'],'w');  
    processFolder(folderName)
    fclose(fileID);
    
    function processFolder(folder)
        disp(['Looking at... ', folder])
        fprintf(fileID,'%s\n',repmat('%',[1,100]));
        fprintf(fileID,'%s\n',folder);
        fprintf(fileID,'%s\n',repmat('-',[1,100]));
        allfiles = dir(folder);allfiles=allfiles(~ismember({allfiles.name},{'.','..'}));
        if prod(cellfun(@isempty,strfind({allfiles.name},'TracesAndEvents'))) && ~prod(cellfun(@isempty,strfind({allfiles.name},'concat_recording')))
           fprintf(fileID,'%s\n',' WARNING, NO TRACE EVENT FILE WAS FOUND FOR THIS DAY');    
        end
        for k = 1:length(allfiles)
            if allfiles(k).isdir
                processFolder([allfiles(k).folder,filesep,allfiles(k).name])
            else
                [~,name,extension] = fileparts(allfiles(k).name);
                if ~isempty(strfind([name,extension],'dfof.h5')) || ~isempty(strfind(name,'concat_recording'))
                    hinf = hdf5info([allfiles(k).folder,filesep,allfiles(k).name]);
                    framesCa = hinf.GroupHierarchy.Datasets.Dims(3);
                    fprintf(fileID,'%s\n',['Calcium ', allfiles(k).name]);
                    fprintf(fileID,'%s\n',['    ',num2str(framesCa),' frames']);
                    
                elseif ~isempty(strfind(name,'TracesAndEvents'))
                    out = load([allfiles(k).folder,filesep,allfiles(k).name]);
                    framesX = size(out.tracesEvents.position,1);
                    framesR = size(out.tracesEvents.rawTraces,1);
                    numCells = size(out.tracesEvents.rawTraces,2);
                    fprintf(fileID,'%s\n',['TracesEvents ', allfiles(k).name]);
                    if framesX~=framesR
                        fprintf(fileID,'%s\n',['    ',' !!!!!!WARNING X AND R DIFFER IN DIMENSIONS!!!!!!']);    
                    end
                    fprintf(fileID,'%s\n',['    ',num2str(framesX),' frames positon']);
                    fprintf(fileID,'%s\n',['    ',num2str(framesR),' frames responses']);
                    fprintf(fileID,'%s\n',['    ',num2str(numCells),' cells detected']);
                    
                elseif strcmp(extension,'.avi')
                    try
                        mov = VideoReader([allfiles(k).folder,filesep,allfiles(k).name]);
                    catch
%                         Miji;%MIJ.run
%                         IJ=ij.IJ();
%                         macro_path = 'C:\Users\csc\Desktop\Fiji.app\macros\aviChangeCodex.txt';
%                         movie = [allfiles(k).folder,filesep,allfiles(k).name];
%                         IJ.runMacroFile(java.lang.String(macro_path),java.lang.String(movie));

                        fprintf(fileID,'%s\n',['COULD NOT READ MOVIE FOR MOVIE ',allfiles(k).name]);
                        continue;
                    end
                    framesBeh = mov.NumberOfFrames;
                    fprintf(fileID,'%s\n',['Behavior ', allfiles(k).name]);
                    fprintf(fileID,'%s\n',['    ',num2str(framesBeh),' frames']);         
                end
            end
        end 
        fprintf(fileID,'\n');
    end
end

