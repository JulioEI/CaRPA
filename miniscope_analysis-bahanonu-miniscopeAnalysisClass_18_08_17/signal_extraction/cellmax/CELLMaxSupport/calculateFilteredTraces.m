function cellTraces = calculateFilteredTraces(imgs, cellImages, varargin)
    % example function with outline for necessary components
    % Lacey Kitsch
    % updated by: biafra ahanonu
    % started: xxxx.xx.xx
    % inputs
        %
    % outputs
        %

    % changelog
        % 2017.01.18 - updated to add support for movie chunking - biafra
    % TODO
        %

    %========================
    options.oneCentered=0;
    options.readMovieChunks = 0;
    options.movieFilename = '';
    options.movieDatasetName = '/Data/Images';

    options = getOptions(options, varargin);
    %========================

    if options.readMovieChunks==0
        nFrames=size(imgs,3);
        numPixels=size(imgs,1)*size(imgs,2);
        nPixelsFrame = numel(imgs(:,:,1));
    else
        movieDims = loadMovieList(options.movieFilename,'getMovieDims',1,'inputDatasetName',options.movieDatasetName,'displayInfo',0);
        yDim = movieDims.one;
        xDim = movieDims.two;
        zDim = movieDims.three;
        frameSize = [yDim xDim];
        nElementsDim = yDim*xDim;

        nFrames=zDim;
        numPixels=yDim*xDim;
        nPixelsFrame = numPixels;
    end

    if numel(size(cellImages))==3
        nCells=size(cellImages,3);
        cellImages=reshape(cellImages, [nPixelsFrame nCells]);
    elseif size(cellImages,1)*size(cellImages,2)==numPixels
        nCells=1;
        cellImages=reshape(cellImages, [nPixelsFrame nCells]);
    else
        nCells=size(cellImages,2);
    end
    cellTraces=zeros(nCells,nFrames);
    for cInd=1:nCells
        cellImages(:,cInd)=cellImages(:,cInd)/max(cellImages(:,cInd));
    end

    frInc=200;
    reverseStr = '';
    for fr=1:frInc:nFrames
        if mod(fr,25)==0
            reverseStr = cmdWaitbar(fr,nFrames,reverseStr,'inputStr','calculating traces frame-by-frame','waitbarOn',1,'displayEvery',25);
        end

        endFrame=min(nFrames,fr+frInc-1);
        frameInds=fr:endFrame;

        if options.readMovieChunks==0
            thisImgs = imgs(:,:,frameInds);
        else
            thisImgs = loadMovieList(options.movieFilename,'frameList',frameInds,'inputDatasetName',options.movieDatasetName,'displayInfo',0);
            thisImgs(isnan(thisImgs)) = 0;
        end

        if options.oneCentered
            thisImgs=thisImgs-1;
        end
        %thisImgs(thisImgs<0)=0;
        thisImgs=reshape(thisImgs,[nPixelsFrame length(frameInds)])';
        theseTraces=thisImgs*cellImages;
        cellTraces(:,frameInds)=theseTraces';
    end
end