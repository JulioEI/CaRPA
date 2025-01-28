function [thresholdedImages boundaryIndices] = thresholdImages(inputImages,varargin)
    % thresholds input images and makes them binary if requested
    % biafra ahanonu
    % started: 2013.10.xx
    % adapted from SpikeE
    %
    % inputs
        %
    % outputs
        %

    % changelog
        % updated: 2013.11.04 [15:30:05] added try...catch block to get around some errors for specific filters
        % 2014.01.14 refactored so it now can handle multiple images instead of just one
        % 2014.01.16 [16:30:36] fixed error after refactoring where thresholdedImage dims were not a 3D matrix, caused assignment errors.
        % 2014.03.13 slight change to support double and other non-integer images
        % 2017.01.14 [20:06:04] - support switched from [nSignals x y] to [x y nSignals]
    % TODO
        %

    %========================
    options.threshold = 0.5;
    options.waitbarOn = 1;
    options.binary = 0;
    % 1 = open workers, 0 = do not open workers
    options.parallel = 1;
    % 1 = get boundary index, 0 = do nothing
    options.getBoundaryIndex = 0;
    % image filter: none, median,
    options.imageFilter = 'none';
    % normalize images
    options.normalizationType = [];
    % get options
    options = getOptions(options,varargin);
    % display(options)
    % unpack options into current workspace
    % fn=fieldnames(options);
    % for i=1:length(fn)
    %     eval([fn{i} '=options.' fn{i} ';']);
    % end
    %========================

    nImages = size(inputImages);
    inputDims = size(inputImages);
    inputDimsLen = length(inputDims);
    if inputDimsLen==3
        nImages = size(inputImages,3);
    elseif inputDimsLen==2
        nImages = 1;
        tmpImage = inputImages; clear inputImages;
        inputImages(:,:,1) = tmpImage;
        options.waitbarOn = 0;
    else
        return
    end
    % loop over all images and threshold
    reverseStr = '';
    % pre-allocate for speed
    thresholdedImages = zeros(size(inputImages));
    boundaryIndices = cell([nImages 1]);
    manageParallelWorkers('parallel',options.parallel);
    for imageNo=1:nImages

        thisFilt = squeeze(inputImages(:,:,imageNo));
        switch options.imageFilter
            case 'median'
                thisFilt = medfilt2(thisFilt,[3 3]);
            otherwise
                % body
        end
        % threshold
        maxVal=nanmax(thisFilt(:));
        cutoffVal = maxVal*options.threshold;
        % cutoffVal
        % cutoffVal = quantile(thisFilt(:),options.threshold);
        replaceVal = 0;
        thisFilt(thisFilt<cutoffVal)=replaceVal;
        thisFilt(isnan(thisFilt))=replaceVal;

        % make image binary
        if options.binary==1
            thisFilt(thisFilt>=cutoffVal)=1;
        else
            % normalize
            thisFilt=thisFilt/maxVal;
        end

        % Remove any pixels not connected to the image max value if there is a filter with max values at the edge, try...catch to get around errors
        try
            [indx indy] = find(thisFilt==1); %Find the maximum
            B = bwlabeln(thisFilt);
            thisFilt(B~=B(indx,indy)) = 0;
        catch
        end
        % size(thisFilt)
        % size(thresholdedImage)
        thresholdedImages(:,:,imageNo)=thisFilt;

        if options.binary==1&options.getBoundaryIndex==1
            [B,L] = bwboundaries(thisFilt);
            for iNo = 1:length(B)
                boundaryIndices{imageNo} = [boundaryIndices{imageNo} sub2ind(size(thisFilt),B{iNo}(:,1),B{iNo}(:,2))'];
            end
            boundaryIndices{imageNo} = boundaryIndices{imageNo}(:)';
        end
        % within loop
        if (mod(imageNo,20)==0|imageNo==nImages)&options.waitbarOn==1
            reverseStr = cmdWaitbar(imageNo,nImages,reverseStr,'inputStr','thresholding images');
        end
    end
    % ensure backwards compatibility
    if nImages==1&inputDimsLen<3
        thresholdedImages = squeeze(thresholdedImages);
    end
    if nImages>1&~isempty(options.normalizationType)
        % thresholdedImages = permute(normalizeMovie(permute(thresholdedImages,[2 3 1]),'normalizationType','zeroToOne'),[3 1 2]);
        thresholdedImages = normalizeMovie(thresholdedImages,'normalizationType','zeroToOne');
    end
end