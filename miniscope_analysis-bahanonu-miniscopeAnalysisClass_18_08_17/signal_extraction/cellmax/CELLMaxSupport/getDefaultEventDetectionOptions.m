function options = getDefaultEventDetectionOptions()
	% put event detection options in a single place for later re-use
	% biafra ahanonu
	% started: 2017.04.22
	% inputs
		%
	% outputs
		%

	% changelog
		%
	% TODO
		%

	%========================
	% options.exampleOption = '';
	% get options
	% options = getOptions(options,varargin);
	% display(options)
	% unpack options into current workspace
	% fn=fieldnames(options);
	% for i=1:length(fn)
	% 	eval([fn{i} '=options.' fn{i} ';']);
	% end
	%========================

	try
		options.framerate=4;
		options.burstDuration=1;
		options.refractoryPeriod=1;
		options.useImageCorrelation=1;
		options.useImageCorrelationTrace=0;
		options.nBaselineFrames=500;
		options.correlationWindow=10;
		options.prctileCorrThresh = 90;
		options.prctileCorr = 95;
		options.prctileDistance = 5;
		options.minThresh = 0.75;
		options.numSigmasThresh=5;
		options.diffThreshold = 0;
		options.phiBaselineThresh = 0.1;
		options.smoothFilter = gaussian2(1,7);

		options.readMovieChunks = 0;
		options.movieFilename = '';
		options.movieDatasetName = '/Data/Images';
		options.displayInfo = 0;
	catch err
		display(repmat('@',1,7))
		disp(getReport(err,'extended','hyperlinks','on'));
		display(repmat('@',1,7))
	end
end