function output = correctMotion_DFT(hdf5_in_name,hdf5_out_name,f_type,varargin)

    % correctMotion_DFT performs motion registration using the DFT method. Each
    % frame in the input file is registered to the first frame and saved in an
    % hdf5.
    %
    % Inputs:
    %   hdf5_in_name: Name of input hdf5 file (including leading directory.
    %   hdf5_out_name: Name of output hdf5 file (including leading directory).
    %   f_type: Filter type (0-4)
    %           0: no filte
    %           1,2: mean subtract (r = 20), blur (r = 2, r = 20)
    %           3: median filter
    %           4: mean subtract, complement
    % Outputs:
    %   output: Structure containing the outputs of the dft registration
    %           [error,phase,row_shift,col_shift]
    %
    % Jessica Maxey, March 23, 2015
    %
    % changelog
    % 2015.10.09 - added varargin options and put structure in place to allow either path or matrix input of movie to register - biafra

    %========================
    options.dataset = '/Data/Images';
    % get options
    options = getOptions(options,varargin);
    % display(options)
    % unpack options into current workspace
    % fn=fieldnames(options);
    % for i=1:length(fn)
    %   eval([fn{i} '=options.' fn{i} ';']);
    % end
    %========================
    dataset = options.dataset;

    inputMovieClass = class(hdf5_in_name);
    if strcmp(inputMovieClass,'char')
        h5att = h5info(hdf5_in_name,dataset);
        imageStackSize = h5att.Dataspace.Size;
        rows = imageStackSize(1);
        cols = imageStackSize(2);
        frames = imageStackSize(3);
        baseImage = im2single(h5read(hdf5_in_name,dataset,[1,1,1],[rows,cols,1]));
    else

    end

    switch(f_type)
        case(0)
            transform = @(A) transform_0(A);
        case(1)
            hDisk  = fspecial('disk', 20);
            hDisk2 = fspecial('disk', 2);
            transform = @(A) transform_1(A,hDisk,hDisk2);
        case(2)
            hDisk  = fspecial('disk', 20);
            hDisk2 = fspecial('disk', 20);
            transform = @(A) transform_2(A,hDisk,hDisk2);
        case(3)
            med_filt = 5;
            transform = @(A) transform_3(A,med_filt);
        case(4)
            hDisk  = fspecial('disk', 20);
            hDisk2 = fspecial('disk', 2);
            transform = @(A) transform_4(A,hDisk,hDisk2);
    end

    im_ref = transform(baseImage);

    % Specify ROI
    %------------------------------------------------------------
    h_roi = figure;
    imagesc(im_ref); axis image; colormap gray;
    title('Select ROI');
    h_rect = imrect;
    mask_rect = round(getPosition(h_rect));
    startRow = mask_rect(2);
    stopRow = mask_rect(2)+mask_rect(4);
    startCol = mask_rect(1);
    stopCol = mask_rect(1)+mask_rect(3);
    cropBase = im_ref(startRow:stopRow,startCol:stopCol);
    close(h_roi);

    baseFFT = fft2(cropBase);

    chunkSize = [rows,cols,1];
    h5create(hdf5_out_name,dataset,[Inf,Inf,Inf],'ChunkSize',chunkSize,'Datatype','single');
    h5write(hdf5_out_name,dataset,baseImage,[1,1,1],[rows,cols,1]);

    h = waitbar(0,'Motion Correction in Progress...');
    output = zeros(4,frames);
    for f=2:frames
        waitbar(f/frames);
        image = im2single(h5read(hdf5_in_name,dataset,[1,1,f],[rows,cols,1]));
        tImage = transform(image);
        cropImage = tImage(startRow:stopRow,startCol:stopCol);
        imageFFT = fft2(cropImage);
        [output(:,f), result] = dftregistration(baseFFT,imageFFT,100);

        %%% Apply shift parameters to the non-blurred/cropped image
        iFFT = fft2(image);
        [nr,nc]=size(iFFT);
        Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
        Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
        [Nc,Nr] = meshgrid(Nc,Nr);
        iFFT = iFFT.*exp(i*2*pi*(-output(3)*Nr/nr-output(4)*Nc/nc));
        regImageFFT = iFFT*exp(i*output(2));
        regImage = abs(ifft2(regImageFFT));

        h5write(hdf5_out_name,dataset,regImage,[1,1,f],[rows,cols,1]);
    end
    close(h);
end

function A_tr = transform_0(A)
    A_tr = A;
end

function A_tr = transform_1(A, ssm_filter, asm_filter)
    A_tr = A - imfilter(A, ssm_filter, 'replicate');
    A_tr = imfilter(A_tr, asm_filter);
end

function A_tr = transform_2(A, ssm_filter, asm_filter)
    A_tr = A - imfilter(A, ssm_filter, 'replicate');
    A_tr = imfilter(A_tr, asm_filter);
end

function A_tr = transform_3(A,med_filt)
    A_tr = medfilt2(A,[med_filt, med_filt]);
    [h,w] = size(A);
    % A_tr = A_tr(2:h-1,2:w-1);
    A_tr = A_tr(1:h,1:w);
end

function A_tr = transform_4(A, ssm_filter, asm_filter)
    A_tr = A - imfilter(A, ssm_filter, 'replicate');
    A_tr = imcomplement(A_tr);
end