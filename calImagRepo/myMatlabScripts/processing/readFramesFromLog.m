function [ missingFrames, totalFrames ] = readFramesFromLog(logFile)

    [~,~,ext] = fileparts(logFile);

    switch ext
        case '.xml'
            missingFrames = getNumDataInLabelInXml('dropped',logFile);
            totalFrames = getNumDataInLabelInXml('frames',logFile);
        case {'.txt','.log'}            
            missingFrames = getNumDataInLineInTxt('DROPPED:',logFile);
            totalFrames = getNumDataInLineInTxt('FRAMES:',logFile);
        otherwise
            error(['File extension: ',ext,' not recognized'])
    end

end

function numData = getNumDataInLineInTxt(lineName,logFile)
    frameFile = importdata(logFile,'');
    frameText = frameFile{cellfun(@(x) ~isempty(x),cellfun(@(x) strfind(x,lineName),frameFile,'UniformOutput',false))};
    numData = cellfun(@str2double,regexp(frameText,'-?\d*','match'));
end

function numData = getNumDataInLabelInXml(findLabel,logFile)
    %findLabel is the X in name = 'X' inside atrr
    frameFile = xmlread(logFile);
    allListitems = frameFile.getElementsByTagName('attr'); %Gets attr labels
    for k = 0:allListitems.getLength-1 %java
       thisListitem = allListitems.item(k);
       if strcmp(thisListitem.getAttribute('name'),findLabel)
           frameText = char(thisListitem.getFirstChild.getData);
           numData = cellfun(@str2double,regexp(frameText,'\d*','match')); 
           break;
       end
    end  
end