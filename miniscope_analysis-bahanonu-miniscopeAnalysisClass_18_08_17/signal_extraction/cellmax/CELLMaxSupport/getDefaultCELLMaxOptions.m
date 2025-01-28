function options = getDefaultCELLMaxOptions()

% display
options.suppressProgressFig=0;
options.plotDeltaParams = 0;

% outputs
options.recalculateFinalTraces=1;
options.calculateFilteredTraces=0;
options.doEventDetect=0;

% movie chunking and restriction
options.maxSqSize=101;
options.xLims=[];
options.yLims=[];
options.sqOverlap=16;

% Initialization
options.initMethod='ica';
options.inputSizeManual=0;
options.gridSpacing=9;%7; % spacing of initialization grid, in pixels
options.gridWidth=4;%8;
options.useInitialElimHeuristic=0;
options.initProbsWithFiltTraces=1;

% Convergence
options.maxDeltaParams=10^(-2.1);
options.maxDeltaDeltaParams=0.0001;
options.minIters=500;%51;
options.maxIters=1175;%120;
options.saveParamConvergeMovie=0;

% Removal of overlapping and small images
options.borderRemoveBuffer=8;
options.sizeThresh=12;
options.removeZeroVarImages=1;
options.removeDiscontig=1;
options.areaOverlapThresh=.65;%.2%.4%0.65;

% scaledPhi-based removal section
options.removeCorrProbs=1;      % toggle for turning this on
options.scaledPhiCorrThresh=0.5;   % correlation threshold
options.distanceThresh=5;   % distance threshold between centroids (pixels)
options.corrRemovalAreaOverlapThresh=0.1;

% Core algorithm parameters, selecting random frame subset
% random, resampleRemaining [default]
options.subsampleMethod = 'resampleRemaining';
% [maxIters nMovieFrames]
options.subsampleFrameMatrix = [];
% [1 nMovieFrames] - vector of frames to use in a movie
options.subsampleFrameVector = [];
options.selectRandomFrames=1;
options.numFramesRandom=2000;
% 0 to 1, percentage of frames per iteration to select
options.percentFramesPerIteration = 0.5;
% subsampleMethod = 'resampleRemaining', fr
options.percentRemainingSubsample = 0.75;
% 1 = randperm and other generators will use a new seed each time, 0 = use seed as specified in randNumGenSeed, will produce deterministic outputs
options.generateNovelSeed = 1;
% Seed to use for rng()
options.randNumGenSeed = 1;

% Core algorithm parameters
options.elimOnScaledPhi=0;
options.numSigmasThresh=0;
options.numPhotonsPerSigma=5;
options.numSigmasThreshInitial=options.numSigmasThresh;
options.removeBelowThreshPixelsForRecalc=1;
options.useMuChangeOnly=0;
options.useScaledPhi=1;
options.threshForElim=0.005;
options.useConstantBG=1;

% Event detect
options.optionsED.numSigmasThresh=3;
options.optionsED.useImageCorrelation = 0;
options.optionsED.diffThresh = 0;

% preset storage fields
options.oneCentered = 0;
options.noiseSigma = 0.004;
options.iterTrack = {};
options.progressSquare = [];
options.localICimgs = [];
options.localICtraces = [];
options.iterOptions.lastIter = 0;
options.iterOptions.notConverged = 1;

% presets for parallel processing
options.useParallel=0;
options.nParallelWorkers = java.lang.Runtime.getRuntime().availableProcessors;
options.closeOpenWorkers = 0;

% load movie or read chunks
options.readMovieChunks = 0;
options.loadPreviousChunks = 0;
options.movieDatasetName='/Data/Images';
options.movieFilename = '';