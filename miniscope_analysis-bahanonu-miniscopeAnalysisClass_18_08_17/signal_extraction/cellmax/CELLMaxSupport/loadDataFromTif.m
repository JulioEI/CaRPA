function [imgs, nFramesByTrial] = loadDataFromTif(pathname, regExp)

if isempty(regExp)
    regExp='recording_*_DFOF.tif';
end

files=dir([pathname filesep regExp]);
imgs=[];

if ~isempty(files)
    nFramesByTrial=zeros(1,length(files));
    for fInd=1:length(files)
        info=imfinfo([pathname filesep files(fInd).name]);
        nFramesByTrial(fInd)=length(info);
        if fInd==1
           imgSize=[info(1).Height info(1).Width];
        elseif info(1).Height~=imgSize(1) || info(1).Width~=imgSize(2)
           error('Tifs are not the same size!')
        end
    end
    nFramesTotal=sum(nFramesByTrial);

    imgs=zeros([imgSize, nFramesTotal], 'single');

    for fInd=1:length(files)
        imgs(:,:,sum(nFramesByTrial(1:(fInd-1)))+(1:nFramesByTrial(fInd)))=loadTifSlow([pathname filesep files(fInd).name]);
    end
else
    disp(['Warning: No files found in ' pathname]);
    nFramesByTrial=[];
end
    
    