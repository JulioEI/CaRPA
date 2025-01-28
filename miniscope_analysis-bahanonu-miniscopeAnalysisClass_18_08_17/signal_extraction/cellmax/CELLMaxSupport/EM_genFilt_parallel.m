function [cellImages, cellTraces, cellParams, noiseSigma, eventTimes, eventTrigImages, scaledPhi, optionsOut, profileObj] =...
    EM_genFilt_parallel(DFOF, varargin)

    % Written by Lacey Kitch in 2013-2014

    % Excutes EM Cell-finding method on input movie. Default options should
    % work fine for most mini-microscope movies. Main options to consider
    % changing are options.initMethod (initialization method),
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

    %mpiprofile on

    %%%%%%%%%% Set Options %%%%%%%%%%%%
    optionsOut = getDefaultCELLMaxOptions();

    %%%% replace default options with input options
    optionsOut = getOptions(optionsOut, varargin);

    %%%%%%%%%% EM %%%%%%%%%%%%

    % check movie class
    if optionsOut.readMovieChunks==0
        if ~isa(DFOF, 'double') && ~isa(DFOF, 'single')
            warning('Movie is not single or double class, unknown behavior could result')
        end
        if isa(DFOF, 'double')
            DFOF=single(DFOF);
        end
        % if present in options, restrict spatial extent of movie
        if ~isempty(optionsOut.xLims) && ~isempty(optionsOut.yLims)
            if length(optionsOut.xLims)==2
                optionsOut.xLims=optionsOut.xLims(1):optionsOut.xLims(2);
            end
            if length(optionsOut.yLims)==2
                optionsOut.yLims=optionsOut.yLims(1):optionsOut.yLims(2);
            end
            DFOF=DFOF(optionsOut.yLims,optionsOut.xLims,:);
        end
    end

    % pre-computer values needed in loop
    if optionsOut.readMovieChunks==0
        yDim = size(DFOF,1);
        xDim = size(DFOF,2);
        zDim = size(DFOF,3);
        movieSize = [yDim xDim zDim];
        % xyDim = size(DFOF(:,:,1));
        nElementsDim = numel(DFOF(:,:,1));
        frameSize = size(DFOF(:,:,1));
    else
        movieDims = loadMovieList(optionsOut.movieFilename,'getMovieDims',1,'inputDatasetName',optionsOut.movieDatasetName);
        yDim = movieDims.one;
        xDim = movieDims.two;
        zDim = movieDims.three;
        frameSize = [yDim xDim];
        movieSize = [yDim xDim zDim];
        nElementsDim = yDim*xDim;
    end

    % Obtain a list of frames to subsample, place here instead of EM_genFilt so that all movie chunks get identical treatment.
    [frameVector frameVectorDecision frameVectorReset] = getSubsampleFrameMatrix(zDim,optionsOut.maxIters+2,optionsOut.percentFramesPerIteration,'subsampleMethod',optionsOut.subsampleMethod,'displayPlots',0,'percentRemainingSubsample',optionsOut.percentRemainingSubsample);
    optionsOut.subsampleFrameMatrix = frameVectorDecision;

    % calculate the std dev of the noise (used for discretizing the movie)
    optionsTmp.readMovieChunks = optionsOut.readMovieChunks;
    optionsTmp.movieFilename = optionsOut.movieFilename;
    optionsTmp.movieDatasetName = optionsOut.movieDatasetName;
    if optionsOut.readMovieChunks==0
        [noiseSigma, movieMean, maxVal] = fitNoiseSigma(DFOF,'options',optionsTmp);
        testFrame=DFOF(:,:,1);
    else
        [noiseSigma, movieMean, maxVal] = fitNoiseSigma(movieSize,'options',optionsTmp);
        testFrame = loadMovieList(optionsOut.movieFilename,'frameList',1:2,'inputDatasetName',optionsOut.movieDatasetName);
        testFrame(isnan(testFrame)) = 0;
        testFrame = testFrame(:,:,1);
    end
    if abs(movieMean)<abs(movieMean-1) || min(testFrame(:))<0
        optionsOut.oneCentered=0;
        disp('Detected that movie is centered around 0. Finding cells...')
    else
        optionsOut.oneCentered=1;
        disp('Detected that movie is centered around 1. Finding cells...')
    end
    optionsOut.noiseSigma=noiseSigma;

    % determine size for field of view chunking, set progress square to not
    % used since it doesn't make sense with parallel squares
    [sqSizeX, sqSizeY, nSqHorz, nSqVert, sqIncX, sqIncY] = findOptimalSquareSize(optionsOut.maxSqSize,optionsOut.sqOverlap,[yDim xDim zDim]);
    optionsOut.suppressProgressFig=1;

    % store all results together
    nSqTotal=nSqVert*nSqHorz;
    cellParams=cell(nSqTotal,1);
    cellImages=cell(nSqTotal,1);
    scaledPhi=cell(nSqTotal,1);

    % =======================================
    % create a cell array that can be used to reduce memory bandwidth usage with Matlab workers
    if optionsOut.readMovieChunks==0
        sprintf('converting DFOF to cell array with %d chunks',nSqTotal);
        for sqInd=1:nSqTotal
                vertInd=mod(sqInd-1,nSqVert)+1;
                horzInd=floor((sqInd-1)/nSqVert)+1;

                % get chunk of data
                yLims=(vertInd-1)*sqIncY+(1:sqSizeY);
                yLims(yLims>size(DFOF,1))=[];
                xLims=(horzInd-1)*sqIncX+(1:sqSizeX);
                xLims(xLims>size(DFOF,2))=[];
                imgsCell{sqInd}=DFOF(yLims,xLims,:);
                thisMovie = imgsCell{sqInd};
                j = whos('thisMovie'); j.bytes = j.bytes*9.53674e-7; display(['movie size: ' num2str(j.bytes) 'Mb | ' num2str(j.size) ' | ' j.class]);
                clear thisMovie;
        end
        j = whos('imgsCell'); j.bytes = j.bytes*9.53674e-7; display(['movie size: ' num2str(j.bytes) 'Mb | ' num2str(j.size) ' | ' j.class]);
        disp('done.')
    else
        imgsCell = cell([nSqTotal 1]);
    end

    % yDim = size(DFOF,1);
    % xDim = size(DFOF,2);
    % zDim = size(DFOF,3);
    % xyDim = size(DFOF(:,:,1));
    % nElementsDim = numel(DFOF(:,:,1));
    % frameSize = size(DFOF(:,:,1));
    % =======================================
    % if using parallel, open parpool and set workers, if there is only a single tile, turn off parallel since less efficient
    if nSqTotal==1
        display('Single tile, disabling parallel operation.');
        nWorkers=0;
    elseif optionsOut.useParallel==1&optionsOut.nParallelWorkers~=0
        manageParallelWorkers('setNumCores',optionsOut.nParallelWorkers);
        nWorkers=Inf;
    else
        nWorkers=0;
    end

    display(repmat('=',1,7))
    display('starting parallel CELLMax runs...')

    % initialize loop variables
    initMethod = optionsOut.initMethod;
    estTraces = cell([1 nSqTotal]);
    localICimgs = cell([1 nSqTotal]);
    localICtraces = cell([1 nSqTotal]);
    nCells = cell([1 nSqTotal]);

    parfor (sqInd=1:nSqTotal,nWorkers)

        % get copy of options to modify if needed, for parfor
        thisOptions=optionsOut;

        if thisOptions.readMovieChunks==1
            % theseImages
            [parentPath,~,~] = fileparts(thisOptions.movieFilename);
            newPath = [parentPath filesep 'tmpImages' filesep];
            if ~exist(newPath,'dir');mkdir(newPath);end
            % tmpFilePath = [newPath filesep 'sqInd_' num2str(sqInd) '.h5'];
            tmpFilePath = [newPath filesep 'sqInd_' num2str(sqInd) '.mat'];

            % don't re-process if file still exists
            if exist(tmpFilePath,'file')&thisOptions.loadPreviousChunks==1
                fprintf('Already ran chunk %d/%d\n',sqInd,nSqTotal);
                continue;
            end
        end

        fprintf('running chunk %d/%d\n',sqInd,nSqTotal);

        vertInd=mod(sqInd-1,nSqVert)+1;
        horzInd=floor((sqInd-1)/nSqVert)+1;


        % get chunk of data
        if thisOptions.readMovieChunks==1
            yLims=(vertInd-1)*sqIncY+(1:sqSizeY);
            yLims(yLims>yDim)=[];
            xLims=(horzInd-1)*sqIncX+(1:sqSizeX);
            xLims(xLims>xDim)=[];
            offset = [min(yLims)-1 min(xLims)-1 0];
            block = [length(yLims) length(xLims) zDim];
            display(['loading chunk | offset: ' num2str(offset) ' | block: ' num2str(block)]);
            [imgs] = readHDF5Subset(thisOptions.movieFilename, offset, block,'datasetName',thisOptions.movieDatasetName);
            imgs(isnan(imgs)) = 0;
        else
            display('using pre-loaded chunk');
            yLims=(vertInd-1)*sqIncY+(1:sqSizeY);
            yLims(yLims>yDim)=[];
            xLims=(horzInd-1)*sqIncX+(1:sqSizeX);
            xLims(xLims>xDim)=[];
            % imgs=DFOF(yLims,xLims,:);
            imgs = imgsCell{sqInd};
        end

        % ensure no NaNs in movies
        display('removing NaNs...');drawnow
        imgs(isnan(imgs)) = 0;

        drawnow('update')

        % j = whos('imgs'); j.bytes=j.bytes*9.53674e-7; display(['movie size: ' num2str(j.bytes) 'Mb | ' num2str(j.size) ' | ' j.class]);

        if ~isempty(imgs)
            % try
                % if initializing with ICA, restrict to local IC images/traces
                % if strcmp(thisOptions.initMethod,'ica') && isfield(thisOptions, 'icImgs') && isfield(thisOptions, 'icTraces') %#ok<*UNRCH>
                %     [thisOptions.localICimgs,thisOptions.localICtraces] = getLocalICs(thisOptions.icImgs,thisOptions.icTraces,xLims,yLims);
                % end

                % if initializing with ICA, restrict to local IC images/traces
                if strcmp(initMethod,'input')
                    fprintf('Using input images/traces at xLims=%d, yLims=%d frames.\n',xLims,yLims);
                    if isfield(thisOptions, 'initImages') && isfield(thisOptions, 'initTraces') %#ok<*UNRCH>
                        [localICimgs{sqInd},localICtraces{sqInd}] = getLocalICs(thisOptions.initImages,thisOptions.initTraces,xLims,yLims);
                    else
                        disp('Warning: Must input images and traces for initialization if you choose the INPUT method. Initializing with grid...')
                        continue;
                        % initMethod = 'grid';
                    end
                elseif strcmp(initMethod,'ica')
                    %
                    nFrames=min(5000,size(imgs,3));
                    nPCs=round(2.5*length(xLims)*length(yLims)/(thisOptions.gridSpacing^2));
                    TermTolICs=10^(-3);
                    fprintf('Running PCA-ICA on %d frames.\n',nFrames);
                    [localICimgs{sqInd},localICtraces{sqInd}] = pcaIca(imgs(:,:,1:nFrames),nPCs,nPCs,0.1,'TermTolICs',TermTolICs);
                else
                    localICimgs{sqInd} = [];
                    localICtraces{sqInd} = [];
                    % initMethod = 'grid';
                end

                % perform EM
                [paramStore, ~, ~] = EM_genFilt(imgs, 'options', thisOptions,'initMethod',initMethod,'localICimgs',localICimgs{sqInd},'localICtraces',localICtraces{sqInd});

                % find very small cells and remove
                estCentroids=paramStore{end}(:,1:2);
                estImages=paramStore{end}(:,2+(1:numel(imgs(:,:,1))));
                cellsToDelete = findSmallImages(estImages, size(imgs(:,:,1)), thisOptions.sizeThresh);
                estCentroids(cellsToDelete,:)=[];
                estImages(cellsToDelete,:)=[];

                % find cells very close to border and remove
                closeCells=false(size(estCentroids,1),1);
                if vertInd>1
                    closeCells(estCentroids(:,2)<=thisOptions.borderRemoveBuffer)=1;
                end
                if horzInd>1
                    closeCells(estCentroids(:,1)<=thisOptions.borderRemoveBuffer)=1;
                end
                if vertInd<nSqVert
                    closeCells(estCentroids(:,2)>=(yLims(end)-yLims(1)-thisOptions.borderRemoveBuffer))=1;
                end
                if horzInd<nSqHorz
                    closeCells(estCentroids(:,1)>=(xLims(end)-xLims(1)-thisOptions.borderRemoveBuffer))=1;
                end
                estCentroids(closeCells,:)=[];
                estImages(closeCells,:)=[];
                if thisOptions.useScaledPhi
                    tmpVar = paramStore{end}(:,2+numel(imgs(:,:,1))+(1:size(imgs,3)));
                    tmpVar(cellsToDelete,:) = [];
                    tmpVar(closeCells,:) = [];
                    estTraces{sqInd} = tmpVar;
                else
                    estTraces{sqInd} = [];
                end

                if ~isempty(estImages)
                    % adjust params for square location and store
                    nCells{sqInd}=size(estCentroids,1);
                    estCentroids(:,1)=estCentroids(:,1)+xLims(1)-1;
                    estCentroids(:,2)=estCentroids(:,2)+yLims(1)-1;
                    [xLims, yLims]=meshgrid(xLims,yLims);
                    linInds=sub2ind(frameSize, yLims, xLims);
                    theseImages=zeros(nCells{sqInd},nElementsDim);
                    theseImages(:,linInds)=estImages;
                    if thisOptions.readMovieChunks==1
                        if exist(tmpFilePath,'file');delete(tmpFilePath);end
                        % [success] = writeHDF5Data(theseImages,tmpFilePath,'datasetname',thisOptions.movieDatasetName);
                        fprintf('saving: %s \n',tmpFilePath);
                        subFxnParsave(tmpFilePath, theseImages,estCentroids,estTraces{sqInd});
                        % save(tmpFilePath,'theseImages','estCentroids','estTracesTmp','-v7.3')
                    else
                        cellImages{sqInd}=theseImages;
                        cellParams{sqInd}=estCentroids;
                        if thisOptions.useScaledPhi
                            scaledPhi{sqInd}=estTraces{sqInd};
                        end

                    end
                end
           % catch err
           %     disp(['Error thrown on vert block ' num2str(vertInd) ', horz block ' num2str(horzInd) ', message: ' err.message])
           % end
        end
    end

    clear imgsCell;

    % =======================================
    % close and re-open to clear memory
    if optionsOut.useParallel==1&optionsOut.nParallelWorkers~=0&optionsOut.closeOpenWorkers==1
        manageParallelWorkers('openCloseParallelPool','close');
        manageParallelWorkers('setNumCores',optionsOut.nParallelWorkers);
        nWorkers=Inf;
    else
        nWorkers=0;
    end
    % =======================================

    % convert params, shapes, and scaled probabilities to arrays
    if optionsOut.readMovieChunks==1
        % cellImages = {};
        cellImages=zeros(sum([nCells{:}]),nElementsDim);
        cellParams=cell(nSqTotal,1);
        scaledPhi=cell(nSqTotal,1);
        %
        for sqInd=1:nSqTotal
            [parentPath,NAME,EXT] = fileparts(optionsOut.movieFilename);
            newPath = [parentPath filesep 'tmpImages' filesep];
            tmpFilePath = [newPath filesep 'sqInd_' num2str(sqInd) '.mat'];
            tmpFilePointer = matfile(tmpFilePath);
            nCells{sqInd}=size(tmpFilePointer.estCentroids,1);
        end
        % size(cellImages)
        for sqInd=1:nSqTotal
            [parentPath,NAME,EXT] = fileparts(optionsOut.movieFilename);
            newPath = [parentPath filesep 'tmpImages' filesep];
            % if ~exist(newPath,'dir');mkdir(newPath);end
            % tmpFilePath = [newPath filesep 'sqInd_' num2str(sqInd) '.h5'];
            tmpFilePath = [newPath filesep 'sqInd_' num2str(sqInd) '.mat'];
            fprintf('loading chunk %d/%d: %s\n',sqInd,nSqTotal,tmpFilePath);
            offset = [0 0];
            block = [nCells{sqInd} nElementsDim];
            tmpFilePointer = matfile(tmpFilePath);
            if sqInd==1
                % cellImages(1:nCells{sqInd},:) = readHDF5Subset(tmpFilePath, offset, block,'datasetName',optionsOut.movieDatasetName);
                cellImages(1:nCells{sqInd},:) = tmpFilePointer.theseImages;
            else
                startIdx = sum([nCells{1:sqInd-1}])+1;
                endIdx = sum([nCells{1:sqInd}]);
                % cellImages(startIdx:endIdx,:) = readHDF5Subset(tmpFilePath, offset, block,'datasetName',optionsOut.movieDatasetName);
                cellImages(startIdx:endIdx,:) = tmpFilePointer.theseImages;
            end

            cellParams{sqInd} = tmpFilePointer.estCentroids;
            if optionsOut.useScaledPhi
                scaledPhi{sqInd} = tmpFilePointer.estTracesTmp;
            end

        end
        cellParams=cell2mat(cellParams);
        scaledPhi=cell2mat(scaledPhi);
        % [parentPath,NAME,EXT] = fileparts(optionsOut.movieFilename);
        % newPath = [parentPath filesep 'tmpImages' filesep];
        % rmdir(newPath)
    else
        cellImages=cell2mat(cellImages);
        cellParams=cell2mat(cellParams);
        scaledPhi=cell2mat(scaledPhi);
    end

    fprintf('size images: [%d %d]\n',size(cellImages,1),size(cellImages,2))

    % chop off the unused portion of param storage
    fprintf('\n=======\nRemoving unused images/traces\n')
    cellParams(isnan(cellParams(:,1)),:)=[];
    cellImages(isnan(cellImages(:,1)),:)=[];
    if optionsOut.useScaledPhi==1
        scaledPhi(isnan(scaledPhi(:,1)),:)=[];
    end

    % reshape and permute cell images
    fprintf('\n=======\nReshaping and permuting cell images\n')
    % cellImages=reshape(cellImages,[size(cellParams,1), size(DFOF(:,:,1))]);
    cellImages=reshape(cellImages,[size(cellParams,1), frameSize]);
    cellImages=permute(cellImages, [2 3 1]);

    % remove the images for the squares that were used as constant background
    if optionsOut.removeZeroVarImages==1
        fprintf('\n=======\nRemoving constant background images\n')
        [cellImages, goodInds] = removeZeroVarImages(cellImages);
        cellParams=cellParams(goodInds,:);
        if optionsOut.useScaledPhi
            scaledPhi=scaledPhi(goodInds,:);
        end
    end

    % remove any discontiguous regions from the images
    if optionsOut.removeDiscontig==1
        fprintf('\n=======\nRemoving discontiguous regions in images\n')
        cellImages=removeDiscontig(cellImages);
    end

    % if we ran EM on more than one square, resolve conflicts and recalculate traces
    if nSqVert>1 || nSqHorz>1
        fprintf('\n=======\nResolve square conflicts and recalculate traces\n')
        if optionsOut.readMovieChunks==0
            [cellImages,cellParams,~,goodCellInds]=resolveBorderConflicts(cellParams,...
                optionsOut.areaOverlapThresh,DFOF,0,cellImages,optionsOut);
        else
            optionsTmp.readMovieChunks = optionsOut.readMovieChunks;
            optionsTmp.movieFilename = optionsOut.movieFilename;
            optionsTmp.movieDatasetName = optionsOut.movieDatasetName;

            [cellImages,cellParams,~,goodCellInds]=resolveBorderConflicts(cellParams,...
                optionsOut.areaOverlapThresh,movieSize,0,cellImages,'options',optionsTmp);
        end

        if optionsOut.useScaledPhi
            scaledPhi=scaledPhi(goodCellInds,:);
        end
    end

    if optionsOut.recalculateFinalTraces==1
        fprintf('\n=======\nRecalculating final traces\n')
        if optionsOut.readMovieChunks==0
            cellTraces = calculateTraces(cellImages, DFOF, 'removeBelowThreshPixelsForRecalc',optionsOut.removeBelowThreshPixelsForRecalc, 'numSigmasThresh',optionsOut.numSigmasThresh, 'noiseSigma',optionsOut.noiseSigma);
        else
            optionsTmp.readMovieChunks = optionsOut.readMovieChunks;
            optionsTmp.movieFilename = optionsOut.movieFilename;
            optionsTmp.movieDatasetName = optionsOut.movieDatasetName;
            optionsTmp.removeBelowThreshPixelsForRecalc = optionsOut.removeBelowThreshPixelsForRecalc;
            optionsTmp.numSigmasThresh = optionsOut.numSigmasThresh;
            optionsTmp.noiseSigma = optionsOut.noiseSigma;
            optionsTmp.movieMaxVal = maxVal;
            cellTraces = calculateTraces(cellImages, movieSize, 'options',optionsTmp);
        end
    else
        cellTraces=[];
    end


    % eliminate any cells that were set as redundant during the svd check in
    % the trace calculation (slightly different from resolve conflicts fnc)
    fprintf('\n=======\nRemove redundant cells\n')
    if optionsOut.recalculateFinalTraces==1
        redunCells=sum(cellTraces,2)==0;
        cellParams(redunCells,:)=[];
        cellImages(:,:,redunCells)=[];
        if ~isempty(cellTraces)
            cellTraces(redunCells,:)=[];
        end
        if optionsOut.useScaledPhi
            scaledPhi(redunCells,:)=[];
        end
    end

    % do final event detection, if options set to do so
    if optionsOut.doEventDetect==1
        fprintf('\n=======\nRun final event detection\n')
        if optionsOut.readMovieChunks==1
            optionsTmp = optionsOut.optionsED;
            optionsTmp.readMovieChunks = optionsOut.readMovieChunks;
            optionsTmp.movieFilename = optionsOut.movieFilename;
            optionsTmp.movieDatasetName = optionsOut.movieDatasetName;
            optionsTmp.oneCentered = optionsOut.oneCentered;
            output.filtTraces=calculateFilteredTraces([], cellImages, 'options', optionsTmp);
            eventTimes=detectEventsOnPhi(double(filtTraces), double(scaledPhi), cellImages, cellParams, movieSize, 'options', optionsTmp);
        else
            filtTraces=calculateFilteredTraces(DFOF, cellImages, 'options', optionsOut);
            eventTimes=detectEventsOnPhi(double(filtTraces), double(scaledPhi), cellImages, cellParams, DFOF, 'options', optionsOut.optionsED);
        end

    else
        eventTimes=[];
    end
    eventTrigImages=[];


    % get profile information
    % mpiprofile off
    % profileObj=mpiprofile('info');
    % mpiprofile clear reset
    profileObj=[];
end

function subFxnParsave(tmpFilePath, theseImages,estCentroids,estTracesTmp)
    save(tmpFilePath,'theseImages','estCentroids','estTracesTmp','-v7.3');
end