function out = load_tif_movie(filename,downsample_xy,varargin)
    %loads filename movie, downsamples in space by factor downsample_xy
    % biafra ahanonu
    % parts adapted from
        % jerome lecoq for spikee
        % http://www.mathworks.com/matlabcentral/answers/108021-matlab-only-opens-first-frame-of-multi-page-tiff-stack
    % updating: 2013.10.22
    % inputs
        %
    % outputs
        %
    % changelog
        % 2014.01.07 [10:14:49] - removed low-level reading, seems to have problems with some versions of tifflib, will fix later.
        % 2015.11.12 Added reading of TIFFs from ImageJ bigger than 4GB.
    % TODO
        %

    %========================
    options.exampleOption = 'doSomething';
    options.Numberframe = [];
    % get options
    options = getOptions(options,varargin);
    % display(options)
    % unpack options into current workspace
    % fn=fieldnames(options);
    % for i=1:length(fn)
    %     eval([fn{i} '=options.' fn{i} ';']);
    % end
    %========================

    %First load a single frame of the movie to get generic information
    TifLink = Tiff(filename, 'r'); %Create the Tiff object
    TmpImage = TifLink.read();%Read in one picture to get the image size and data type
    TifLink.close(); clear TifLink

    LocalImage = imresize(TmpImage, 1/downsample_xy); clear TmpImage; %Resize
    SizeImage=size(LocalImage);%xy dimensions
    ClassImage= class(LocalImage); clear LocalImage; %Get the class of the movie

    % Pre-allocate the movie
    if isempty(options.Numberframe)
        out.Numberframe=size(imfinfo(filename),1);% Number of frames
        framesToGrab = 1:out.Numberframe;
    else
        out.Numberframe = options.Numberframe;
        framesToGrab = out.Numberframe;
    end
    out.Movie =zeros(SizeImage(1),SizeImage(2),out.Numberframe,ClassImage);

    % determine whether standard or non-standard TIFF
    if numel(framesToGrab)>1
        standardTIFF()
    else
        nonstandardTIFF()
    end

    function standardTIFF()
        info = imfinfo(filename);
        num_images = numel(info);

        reverseStr = '';
        for frame = framesToGrab
            out.Movie(:,:,frame) = imread(filename, frame);
            reverseStr = cmdWaitbar(frame,numel(framesToGrab),reverseStr,'inputStr','loading non-ImageJ tif','waitbarOn',1,'displayEvery',50);
        end
    end
    function nonstandardTIFF()
        fileInfo = imfinfo(filename);
        if isfield(fileInfo,'ImageDescription')==1
            numFramesStr = regexp(fileInfo.ImageDescription, 'images=(\d*)', 'tokens');
            nFrames = str2double(numFramesStr{1}{1});
        else
            nFrames = size(fileInfo,1);
            fileInfo = fileInfo(1);
        end
        % open file handle, read file using low level I/O functions
        fileID = fopen(filename , 'rb');
        % The StripOffsets field provides the offset to the first strip. Based on
        % the INFO for this file, each image assumed to consist of 1 strip.
        fseek(fileID, fileInfo.StripOffsets, 'bof');
        % It is assumed that the images are stored one after the other.
        frameNo = nFrames;
        out.Movie =zeros(fileInfo.Height,fileInfo.Width,frameNo,ClassImage);

        switch fileInfo.ByteOrder
            case 'little-endian'
                byteorder = 'ieee-le';
            case 'big-endian'
                byteorder = 'ieee-be';
            otherwise
                % body
        end

        fileType = ClassImage;

        reverseStr = '';
        for frameNo = 1:nFrames
            out.Movie(:,:,frameNo) = fread(fileID, [fileInfo.Width fileInfo.Height], fileType, 0, byteorder)';
            reverseStr = cmdWaitbar(frameNo,nFrames,reverseStr,'inputStr','loading ImageJ tif','waitbarOn',1,'displayEvery',50);
        end

        % close handle to file
        fclose(fileID);
    end
    function archiveFunction()
        % % We use low-level access to the tifflib library file to avoid duplicating
        % % Access to the Tif properties while reading long list of directories in Tiffs
        % FileID = tifflib('open',filename,'r');
        % rps = tifflib('getField',FileID,Tiff.TagID.RowsPerStrip);
        % hImage = tifflib('getField',FileID,Tiff.TagID.ImageLength);
        % rps = min(rps,hImage);

        % nrows = 0;
        % for j=framesToGrab
        %     tifflib('setDirectory',FileID,1+j-1);
        %     % Go through each strip of data.
        %     for r = 1:rps:hImage
        %         row_inds = r:min(hImage,r+rps-1);
        %         stripNum = tifflib('computeStrip',FileID,r);
        %         display([num2str(r) ' | ' num2str(row_inds) ' | ' num2str(stripNum)])
        %         if downsample_xy~=1
        %             TmpImage(row_inds,:) = tifflib('readEncodedStrip',FileID,stripNum);
        %         else
        %             nrows = nrows + size(tifflib('readEncodedStrip',FileID,stripNum),1)
        %             out.Movie(row_inds,:,j)= tifflib('readEncodedStrip',FileID,stripNum);
        %         end
        %     end
        %     if downsample_xy~=1
        %         out.Movie(:,:,j)=imresize(TmpImage,1/downsample_xy);
        %     end
        % end
        % tifflib('close',FileID);
    end
end