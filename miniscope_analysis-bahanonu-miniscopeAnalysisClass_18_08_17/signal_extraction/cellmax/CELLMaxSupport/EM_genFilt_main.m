function [cellImages, cellTraces, cellParams, noiseSigma, eventTimes, eventTrigImages, scaledPhi, options] =...
    EM_genFilt_main(DFOF, varargin)

    % Written by Lacey Kitch in 2013-2014

    % Executes EM Cell-finding method on input movie. Default options should
    % work fine for most mini-microscope movies. Main options to consider
    % changing are options.initMethod (initialization method)
    % and options.maxDeltaParams (for convergence speed)


    % Inputs:
    % DFOF: the data. should be a 3D matrix, space x space x time.
    %       - Should be motion corrected and DFOF
    %       - Running on 5Hz data is best for noise elimination and runtime.
    % varargin : use EM_genFilt_main(DFOF, 'options', options) to input an
    %   options structure. See below for definition of options and defaults.

    % Outputs:
    % cellImages: images of the estimated shape of each cell.
    %    array size: nypixels x nxpixels x total # cells.
    % cellTraces: estimated fluorescence values for all cells.
    %    array size: total # cells x # frames.
    % cellParams: centroid of estimated cell shape
    % noiseSigma: the std dev of the noise of the whole movie
    % eventTrigImages: defunct output, is empty
    % eventTimes: cell containing events for each trace, as calculated by
    %   detectEvents

    % varargin: options struct.
    %   options.suppressProgressFig - set to 1 to suppress the figure that
    %       displays progress in the movie
    %   options.initMethod - method of initialization.
    %       'grid' (default) - initializes with an evenly spaced grid of
    %           gaussian blobs.
    %       'random' - initializes with random images
    %       'vertStripes' - initializes with vertical blurry stripes (for
    %           imaging purkinje cells)
    %       'ica' - initializes to images from ICA output. requires you to
    %           input these in options.icImgs and options.icTraces.
    %   options.maxDeltaParams - conversion threshold
    %   options.maxSqSize - max size of the chunk of data that the algorithm will work with at
    %       one time. sqSize=50 means the algorithm will run on no larger than a
    %       50x50 square chunk. Optimal for speed, RAM, and accuracy is about 60.
    %       Using less than 35 will be inaccurate.
    %   options.recalculateFinalTraces - Set to 1 (default) to calculate traces
    %       at end.
    %	options.icImgs - images of ICA output to use in EM initialization. (optional)
    %	options.icTraces - traces of ICA output to use in EM initialization. (optional)
    %   options.doEventDetect - set to 0 to skip event detection at the end
    %       (default is 1)
    %   options.optionsED - options structure for event detection (see
    %       detectEvents)
    %   options.recalculateFinalTraces - set to 0 to skip recalculating the
    %       final traces after conflict resolution at the end, and to skip trace
    %       calculation at the end of each iteration - ONLY set this to 0 if you
    %       will be recalculating the traces with a different movie (ie full time
    %       resolution) later (default is 1)



    %%%%%%%%%% Set Options %%%%%%%%%%%%
    options = getDefaultCELLMaxOptions();

    %%%% replace default options with input options
    options=getOptions(options, varargin);



    %%%%%%%%%% EM %%%%%%%%%%%%

    % check movie class
    if ~isa(DFOF, 'double') && ~isa(DFOF, 'single')
        warning('Movie is not single or double class, unknown behavior could result')
    end
    if isa(DFOF, 'double')
        DFOF=single(DFOF);
    end

    % if set to yes, display max projection to select expected size of cells
    if options.inputSizeManual
        h=figure;
        maxProj=max(DFOF,[],3);
        dispylims=max(round((size(DFOF,1)/2)-options.maxSqSize),1):min(round((size(DFOF,1)/2)+options.maxSqSize), size(DFOF,1));
        dispxlims=max(round((size(DFOF,2)/2)-options.maxSqSize),1):min(round((size(DFOF,2)/2+options.maxSqSize)), size(DFOF,2));
        imagesc(maxProj(dispylims, dispxlims)); hax=gca;
        title('Select a region covering one cell.')
        cell1=imellipse(hax);
        img1=createMask(cell1);
        options.gridWidth=max(ceil(sqrt(sum(img1(:)==1)/pi)),3);
        title('Select the closest neighboring cell.')
        cell2=imellipse(hax);
        pos1=getPosition(cell1);
        pos2=getPosition(cell2);
        options.gridSpacing=max(ceil(norm(pos1(1:2)-pos2(1:2)))+1,5);
        options.sqOverlap=2*options.gridSpacing;
        options.borderRemoveBuffer=ceil(options.sqOverlap/2);
        close(h); drawnow
    end
    % options.maxSqSize=12*options.gridSpacing;
    options.maxSqSize(options.maxSqSize<60)=60; options.maxSqSize(options.maxSqSize>90)=90;

    % if present in options, restrict spatial extent of movie
    if ~isempty(options.xLims) && ~isempty(options.yLims)
        if length(options.xLims)==2
            options.xLims=options.xLims(1):options.xLims(2);
        end
        if length(options.yLims)==2
            options.yLims=options.yLims(1):options.yLims(2);
        end
        DFOF=DFOF(options.yLims,options.xLims,:);
    end


    % calculate the std dev of the noise (used for discretizing the movie)
    [noiseSigma, movieMean] = fitNoiseSigma(DFOF,options);
    testFrame=DFOF(:,:,1);
    if abs(movieMean)<abs(movieMean-1) || min(testFrame(:))<0
        options.oneCentered=0;
        disp('Detected that movie is centered around 0. Finding cells...')
    else
        options.oneCentered=1;
        disp('Detected that movie is centered around 1. Finding cells...')
    end
    options.noiseSigma=noiseSigma;


    % determine size for field of view chunking
    [sqSizeX, sqSizeY, nSqHorz, nSqVert, sqIncX, sqIncY] = findOptimalSquareSize(options.maxSqSize,options.sqOverlap,DFOF);
    if ~options.suppressProgressFig
        progressSquare=zeros(nSqVert,nSqHorz);
    else
        progressSquare=[];
    end


    % store all results together
    maxCellsPerSq=round(1.3*sqSizeX*sqSizeY/(options.gridSpacing^2));
    cellParams=nan(maxCellsPerSq*nSqHorz*nSqVert,2);
    cellImages=nan(maxCellsPerSq*nSqHorz*nSqVert,numel(DFOF(:,:,1)));
    scaledPhi=nan(maxCellsPerSq*nSqHorz*nSqVert,size(DFOF,3));
    nCellsSoFar=0;
    options.iterTrack=cell(nSqVert,nSqHorz);
    for vertInd=1:nSqVert
        for horzInd=1:nSqHorz

            % update progress square and output
            if ~options.suppressProgressFig
                progressSquare(vertInd,horzInd)=0.5;
                options.progressSquare=progressSquare;
            end

            % get chunk of data
            yLims=(vertInd-1)*sqIncY+(1:sqSizeY);
            yLims(yLims>size(DFOF,1))=[];
            xLims=(horzInd-1)*sqIncX+(1:sqSizeX);
            xLims(xLims>size(DFOF,2))=[];
            imgs=single(DFOF(yLims,xLims,:));

            if ~isempty(imgs)
                %try
                    % if initializing with ICA, restrict to local IC images/traces
                    if strcmp(options.initMethod,'input')
                        if isfield(options, 'initImages') && isfield(options, 'initTraces') %#ok<*UNRCH>
                            [options.localInitImages,options.localInitTraces] = getLocalICs(options.initImages,options.initTraces,xLims,yLims);
                        else
                            disp('Warning: Must input images and traces for initialization if you choose the INPUT method. Initializing with grid...')
                            options.initMethod='grid';
                        end
                    elseif strcmp(options.initMethod,'ica')
                        nFrames=min(5000,size(imgs,3));
                        nPCs=round(2.5*length(xLims)*length(yLims)/(options.gridSpacing^2));
                        icaOptions.TermTolICs=10^(-3);
                        [options.localICimgs,options.localICtraces] = pcaIca(imgs(:,:,1:nFrames),nPCs,nPCs,0.1,'options',icaOptions);
                    end

                    % perform EM
                    [paramStore, changeParams, options] = EM_genFilt(imgs, 'options', options);
                    options.iterTrack{vertInd,horzInd}=zeros(2,length(changeParams));
                    options.iterTrack{vertInd,horzInd}(1,:)=changeParams;
                    options.iterTrack{vertInd,horzInd}(2,:)=cellfun(@(x) size(x,1), paramStore(1:length(changeParams)));

                    % find very small cells and remove
                    estCentroids=paramStore{end}(:,1:2);
                    estImages=paramStore{end}(:,2+(1:numel(imgs(:,:,1))));
                    cellsToDelete = findSmallImages(estImages, size(imgs(:,:,1)), options.sizeThresh);
                    estCentroids(cellsToDelete,:)=[];
                    estImages(cellsToDelete,:)=[];

                    % find cells very close to border and remove
                    closeCells=false(size(estCentroids,1),1);
                    if vertInd>1
                        closeCells(estCentroids(:,2)<=options.borderRemoveBuffer)=1;
                    end
                    if horzInd>1
                        closeCells(estCentroids(:,1)<=options.borderRemoveBuffer)=1;
                    end
                    if vertInd<nSqVert
                        closeCells(estCentroids(:,2)>=(yLims(end)-yLims(1)-options.borderRemoveBuffer))=1;
                    end
                    if horzInd<nSqHorz
                        closeCells(estCentroids(:,1)>=(xLims(end)-xLims(1)-options.borderRemoveBuffer))=1;
                    end
                    estCentroids(closeCells,:)=[];
                    estImages(closeCells,:)=[];
                    if options.useScaledPhi
                        estTraces=paramStore{end}(:,2+numel(imgs(:,:,1))+(1:size(imgs,3)));
                        estTraces(cellsToDelete,:)=[];
                        estTraces(closeCells,:)=[];
                    end


                    if ~isempty(estImages)
                        % adjust params for square location and store
                        nCells=size(estCentroids,1);
                        estCentroids(:,1)=estCentroids(:,1)+xLims(1)-1;
                        estCentroids(:,2)=estCentroids(:,2)+yLims(1)-1;
                        [xLims, yLims]=meshgrid(xLims,yLims);
                        linInds=sub2ind(size(DFOF(:,:,1)), yLims, xLims);
                        theseImages=zeros(nCells,numel(DFOF(:,:,1)));
                        theseImages(:,linInds)=estImages;
                        cellImages(nCellsSoFar+(1:nCells),:)=theseImages;
                        cellParams(nCellsSoFar+(1:nCells),:)=estCentroids;
                        if options.useScaledPhi
                            scaledPhi(nCellsSoFar+(1:nCells),:)=estTraces;
                        end
                        nCellsSoFar=nCellsSoFar+nCells;
                    end


        %        catch err
        %            disp(['Error thrown on vert block ' num2str(vertInd) ', horz block ' num2str(horzInd) ', message: ' err.message])
        %        end
            end

            progressSquare(vertInd,horzInd)=1;
        end
    end


    % update the progress figure, if using
    if ~options.suppressProgressFig
        try %#ok<TRYNC>
            figure(777); hold off; imagesc(progressSquare); title('Done with shape estimation, resolving conflicts and calculating final traces...'); drawnow
        end
    end


    % chop off the unused portion of param storage
    cellParams(isnan(cellParams(:,1)),:)=[];
    cellImages(isnan(cellImages(:,1)),:)=[];
    if options.useScaledPhi
        scaledPhi(isnan(scaledPhi(:,1)),:)=[];
    end

    % reshape and permute cell images
    cellImages=reshape(cellImages,[size(cellParams,1), size(DFOF(:,:,1))]);
    cellImages=permute(cellImages, [2 3 1]);

    % remove the images for the squares that were used as constant background
    if options.removeZeroVarImages
        [cellImages, goodInds] = removeZeroVarImages(cellImages);
        cellParams=cellParams(goodInds,:);
        if options.useScaledPhi
            scaledPhi=scaledPhi(goodInds,:);
        end
    end

    % remove any discontiguous regions from the images
    if options.removeDiscontig
        cellImages=removeDiscontig(cellImages);
    end

    % if we ran EM on more than one square, resolve conflicts and reculate
    % traces
    if nSqVert>1 || nSqHorz>1
        [cellImages,cellParams,~,goodCellInds]=resolveBorderConflicts(cellParams,...
            options.areaOverlapThresh,DFOF,0,cellImages,options);
        if options.useScaledPhi
            scaledPhi=scaledPhi(goodCellInds,:);
        end
    end
    if options.recalculateFinalTraces
        cellTraces = calculateTraces(cellImages, DFOF, 'removeBelowThreshPixelsForRecalc',options.removeBelowThreshPixelsForRecalc, 'numSigmasThresh',options.numSigmasThresh, 'noiseSigma',options.noiseSigma);
    end


    % eliminate any cells that were set as redundant during the svd check in
    % the trace calculation (slightly different from resolve conflicts fnc)
    if options.recalculateFinalTraces
        redunCells=sum(cellTraces,2)==0;
        cellParams(redunCells,:)=[];
        cellImages(:,:,redunCells)=[];
        if ~isempty(cellTraces)
            cellTraces(redunCells,:)=[];
        end
        if options.useScaledPhi
            scaledPhi(redunCells,:)=[];
        end
    else
        cellTraces=[];
    end

    % do final event detection, if options set to do so
    if options.doEventDetect
        filtTraces=calculateFilteredTraces(DFOF, cellImages, 'oneCentered', options.oneCentered);
        eventTimes=detectEventsOnPhi(double(filtTraces), double(scaledPhi), cellImages, cellParams, DFOF, 'numSigmasThresh',options.optionsED.numSigmasThresh,'useImageCorrelation',options.optionsED.useImageCorrelation);
    else
        eventTimes=[];
    end
    eventTrigImages=[];

    % update the progress figure, if using
    if ~options.suppressProgressFig
        try %#ok<TRYNC>
            close(777)
        end
    end
end