function imgs=loadTifSlow(fName)

info=imfinfo(fName);

nFrames=length(info);

if ~isfield(info(1), 'SampleFormat') || strcmp(info(1).SampleFormat, 'Unsigned integer')
    imgs=zeros([info(1).Height, info(1).Width, nFrames],'uint16');
else
    imgs=zeros([info(1).Height, info(1).Width, nFrames],'single');
end

for fr=1:nFrames
    imgs(:,:,fr)=imread(fName, fr);
end
