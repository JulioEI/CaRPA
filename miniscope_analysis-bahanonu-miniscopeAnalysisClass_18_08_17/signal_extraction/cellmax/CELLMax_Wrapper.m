function [output, dsImgs] = CELLMax_Wrapper(movieFilename, varargin)

    % CELLMax : Wrapper function
    % Written by Lacey Kitch in 2013-2015

    % ------------------------------------------------------------------------
    % Inputs
    % movieFilename : filename, including path, of the temporally downsampled
    %   movie. The movie should be normalized framewise (ie dividing by
    %   lowpass filter) and DF/F.
    %   Anything in the range 3Hz-10Hz will work. Higher framerate is fine but
    %   will be slower.
    %
    % options : options structure. all inputs optional.
    %   options.movieDatasetName : if using an hdf5 file, dataset name of
    %       the downsampled movie in the hdf5 file
    %   options.CELLMaxOptions : options structure for main CELLMax function. Defaults
    %       should be fine for most purposes, but see EM_genFilt_main for
    %       description of possible options.
    %   options.eventOptions : options structure for event detection. See
    %       detectEventsOnPhi.m for parameters. Main options are
    %       burstDuration (the length of time that an event must remain above
    %       threshold on average, in seconds) and framerate (in Hz)


    % Outputs
    % output : Output structure
    %      output.cellImages : images representing sources found (candidate cells). not all
    %        will be cells. Size is [x y numCells]
    %      output.centroids : centroids of each cell image, x (horizontal) and
    %        then y (vertical). Size is [numCells 2]
    %      output.cellTraces : fluorescence traces for each cell
    %        Size is [numCells numFrames] for numFrames of full
    %        movie
    %      output.scaledProbability : a scaled probability trace for each cell
    %        Can be used as a denoised fluorescence trace.
    %      output.eventTimes : event timings.
    %        This function calculates event times based on the scaled
    %        probability traces output by EM.
    %      output.filtTraces : traces calculated by applying the CELLMax output
    %        images as filters to the movie. Will be noisier than output.cellTraces.
    %      output.CELLMaxoptions : options that CELLMax was run with. Good to keep for
    %        recordkeeping purposes.
    %      output.eventOptions : options that event detection was run with.
    %        Good to keep for recordkeeping purposes.
    % dsImgs : The movie as a single MATLAB array. Warning - this
    %   occupies a lot of RAM!
    %   ****Delete this variable if you do not need it.
    % ------------------------------------------------------------------------

    % get options
    options.movieDatasetName='/Data/Images';
    options.useParallel=0;
    options.useImageCorrEventDetect=0;

    % THK: Are the next two lines needed?
    options.eventOptions.framerate=10;
    options.eventOptions.numSigmasThresh=3;

    % CellMax options
    options.CELLMaxoptions=getDefaultCELLMaxOptions();

    % mandate some options
    options.CELLMaxoptions.recalculateFinalTraces=1;
    options.CELLMaxoptions.useScaledPhi=1;
    options.CELLMaxoptions.doEventDetect=0;

    % THK: Presumably these lines are not needed, since 'doEventDetect=0'?
    % options.CELLMaxoptions.optionsED.numSigmasThresh = 3;
    % options.CELLMaxoptions.optionsED.useImageCorrelation = 0;
    % options.CELLMaxoptions.optionsED.diffThresh = 0;

    % replace default options with input options
    options = getOptions(options, varargin);

    options.CELLMaxoptions.gridSpacing=9;%7; % spacing of initialization grid, in pixels % CUSTOM EDIT!!!!!!
    options.CELLMaxoptions.gridWidth=4;%8; % CUSTOM EDIT!!!!!! 

    % update CELLMax options
    options.CELLMaxoptions.useParallel = options.useParallel;
    options.CELLMaxoptions.movieDatasetName = options.movieDatasetName;
    if ischar(movieFilename)
        options.CELLMaxoptions.movieFilename = movieFilename;
    end
    % check proper configuration for parallel + reading movie chunk in
    if options.CELLMaxoptions.useParallel==1 && options.CELLMaxoptions.readMovieChunks==1 && isempty(options.CELLMaxoptions.movieFilename)
        display('improper parameter configuration')
        output.runCompletedNoErrors = 0;
        dsImgs = [];
        return
    end

    % =======================================
    % if using parallel, open parpool and set workers
    if options.CELLMaxoptions.useParallel==1&options.CELLMaxoptions.nParallelWorkers~=0
        manageParallelWorkers('setNumCores',options.CELLMaxoptions.nParallelWorkers);
        nWorkers=Inf;
    else
        nWorkers=0;
    end

    % version number for this build of CELLMax. Format is YYYYMMDD.majorVersion.minorVersion.bugFix
    output.versionCellmax = '20160503.0.1.0';
    % store filenames in output
    output.movieFilename=movieFilename;

    % whether run was successful
    output.runCompletedNoErrors = 0;

    % pre-open progress figure
    try
        % Don't open progress figure in parallel mode, workers have no graphical display
        if options.useParallel==1
            options.CELLMaxoptions.suppressProgressFig = 1;
        end

        if options.CELLMaxoptions.suppressProgressFig==0
            figNo = 777;
            if ishandle(figNo)
                set(0,'CurrentFigure',figNo)
                figHandle = figNo;
            else
                figHandle = figure(figNo);
            end
        end
    catch err
        display(repmat('@',1,7))
        disp(getReport(err,'extended','hyperlinks','on'));
        display(repmat('@',1,7))
    end

    try
        % load temporally downsampled data
        if options.CELLMaxoptions.readMovieChunks==1
            display('Using movie chunks, don''t pre-load movie')
            dsImgs = [];
            movieDims = loadMovieList(options.CELLMaxoptions.movieFilename,'getMovieDims',1,'inputDatasetName',options.CELLMaxoptions.movieDatasetName);
            movieDims.one = movieDims.x; %CUSTOM_CODE!!!
            movieDims.two = movieDims.y; %CUSTOM_CODE!!!
            movieDims.three = movieDims.z; %CUSTOM_CODE!!!
            display(['Movie dimensions: ' num2str([movieDims.one movieDims.two movieDims.three])])
        else
            if ischar(movieFilename)
                disp('Loading data...')
                if strcmp(movieFilename(end-2:end), 'tif')
                    dsImgs=loadTifSlow(movieFilename);
                else
                    dsImgs = h5read(movieFilename, options.movieDatasetName);
                end
                % store filenames in output
                output.movieFilename=movieFilename;
            else
                disp('Using input matrix...')
                dsImgs = movieFilename;
                clear movieFilename;
                % store filenames in output
                output.movieFilename='raw matrix input';
            end

            % replace any NaNs with zero to avoid problems later
            display('removing NaNs...');drawnow
            dsImgs(isnan(dsImgs)) = 0;
            disp('Done loading data. Running CELLMax...')
            display(['Movie dimensions: ' num2str(size(dsImgs))])
        end

        % run CELLMax
        tic
        if options.useParallel||options.CELLMaxoptions.readMovieChunks==1
            [output.cellImages, output.cellTraces, output.centroids, ~, ~, ~, output.scaledProbability, output.CELLMaxoptions] =...
                EM_genFilt_parallel(dsImgs, 'options', options.CELLMaxoptions);
        else
            [output.cellImages, output.cellTraces, output.centroids, ~, ~, ~, output.scaledProbability, output.CELLMaxoptions] =...
                EM_genFilt_main(dsImgs, 'options', options.CELLMaxoptions);
        end

        if options.CELLMaxoptions.calculateFilteredTraces==1
            if options.CELLMaxoptions.readMovieChunks==1
                optionsTmp.readMovieChunks = options.CELLMaxoptions.readMovieChunks;
                optionsTmp.movieFilename = options.CELLMaxoptions.movieFilename;
                optionsTmp.movieDatasetName = options.CELLMaxoptions.movieDatasetName;
                optionsTmp.oneCentered = options.CELLMaxoptions.oneCentered;
                output.filtTraces=calculateFilteredTraces([], output.cellImages, 'options', optionsTmp);
            else
                output.filtTraces=calculateFilteredTraces(dsImgs, output.cellImages, 'oneCentered', options.CELLMaxoptions.oneCentered);
            end
        else
            output.filtTraces = [];
        end
        t = toc;
        displayString=sprintf('Done running CELLMax. This run took %.3f hours. Detecting events...',t/60/60);
        disp(displayString)


        if options.useImageCorrEventDetect
            if options.CELLMaxoptions.readMovieChunks==1
                optionsTmp = options.eventOptions;
                optionsTmp.readMovieChunks = options.CELLMaxoptions.readMovieChunks;
                optionsTmp.movieFilename = options.CELLMaxoptions.movieFilename;
                optionsTmp.movieDatasetName = options.CELLMaxoptions.movieDatasetName;
                optionsTmp.oneCentered = options.CELLMaxoptions.oneCentered;
                [output.eventTimes, ~, output.eventOptions] = detectEventsOnPhi(output.filtTraces, output.scaledProbability,...
                    output.cellImages, output.centroids, [movieDims.one movieDims.two movieDims.three], 'options', optionsTmp);
            else
                [output.eventTimes, ~, output.eventOptions] = detectEventsOnPhi(output.filtTraces, output.scaledProbability,...
                    output.cellImages, output.centroids, dsImgs, 'options', options.eventOptions);
            end
        else
            [output.eventTimes, ~, ~, ~, ~, output.eventOptions] = detectEvents(double(output.scaledProbability), 'options', options.eventOptions);
        end
        % DFOF=dsImgs;
        totalTime=round(toc);
        displayString=sprintf('Done with all, ready for manual or automated cell classification. This run took %.1f hours.',totalTime/60/60);
        disp(displayString)
        output.runtime = totalTime;
        output.options = options;
        % run completed successfully
        output.runCompletedNoErrors = 1;
    catch err
        display(repmat('@',1,7))
        disp(getReport(err,'extended','hyperlinks','on'));
        display(repmat('@',1,7))
    end
end