function [scaledPhi, filtTraces, eventTimes] = recalcPhiAndDetectEvents(DFOF, cellImages, dsScaledPhi, CELLMaxOptions, varargin)
	% example function with outline for necessary components
	% Lacey Kitsch & biafra ahanonu
	% started: xxxx.xx.xx
	% inputs
	    %
	% outputs
	    %

	% changelog
	    % 2016.09.06 - altered to allow changing of noiseSigma and other values.
	    % 2016.12.06 - added identical movie feature to allow fast copy when movies are the same but varying some other variable
	% TODO
	    %

    %========================
    % number of iterations to run new scaledPhi estimate
    options.nIterations = 20;
    % whether movie is 1 centered, default classic dfof
    options.oneCentered = 0;
    % whether to run event detection
    options.runEventDetection = 1;
    % whether to read a chunk of the movie and the size of the chunked movie
    options.readMovieChunks = 0;
    options.movieFilename = '';
    options.movieDatasetName = '/Data/Images';
    % change EDoptions
    options.EDoptions = getDefaultEventDetectionOptions();
    % get options
    options = getOptions(options,varargin);
    % display(options)
    % unpack options into current workspace
    % fn=fieldnames(options);
    % for i=1:length(fn)
    % 	eval([fn{i} '=options.' fn{i} ';']);
    % end
    %========================

    try
        foptions.oneCentered = options.oneCentered;
        foptions.readMovieChunks = options.readMovieChunks;
        foptions.movieFilename = options.movieFilename;
        foptions.movieDatasetName = options.movieDatasetName;
    	filtTraces = calculateFilteredTraces(DFOF, cellImages,'options',foptions);

		scaledPhi = calculateScaledPhi(cellImages, dsScaledPhi, DFOF, CELLMaxOptions,'nIterations',options.nIterations);

        if options.runEventDetection==1
		  % options = getDefaultEventDetectionOptions()
		  % EDoptions=getOptions(EDoptions, varargin);
		  eventTimes = detectEventsOnPhi(filtTraces, scaledPhi, 'options', options.EDoptions);
        else
            eventTimes = [];
        end
    catch err
    	display(repmat('@',1,7))
    	disp(getReport(err,'extended','hyperlinks','on'));
    	display(repmat('@',1,7))
    end

end