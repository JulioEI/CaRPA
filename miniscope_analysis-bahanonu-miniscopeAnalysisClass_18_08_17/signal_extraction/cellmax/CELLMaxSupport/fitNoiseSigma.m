function [noiseSigma, movieMean, maxVal] = fitNoiseSigma(imgs,varargin)

% Written by Lacey Kitch in 2013
% Updated by Maggie Carr Larkin to ignore pixels with zero variation in time (black edges, dead pixels) in miniscope data
% 2017.01.18 - updated Biafra, added support for movie chunking

% vlm=0;
% if ~isempty(varargin)
%     options=varargin{1};
%     if isfield(options, 'vlm')
%         vlm=1;
%     end
% end
% ========================
options.vlm = 0;
options.readMovieChunks = 0;
options.movieFilename = '';
options.movieDatasetName = '/Data/Images';

%%%% replace default options with input options
options = getOptions(options, varargin);

options

if options.readMovieChunks==0
    nFrames=size(imgs,3);
else
    nFrames = imgs(3);
end

if nFrames>1000
    framesForEst=randperm(nFrames,1000);
else
    framesForEst=1:nFrames;
end

if ~options.vlm
    if options.readMovieChunks==0
        % fInc=(max(max(max(imgs)))-min(min(min(imgs))))/1000;
        maxVal = max(imgs(:));
        minVal = min(imgs(:));
        fInc=(maxVal-minVal)/1000;
        % fInc=(max(imgs(:))-min(imgs(:)))/1000;
    else
        % reverseStr = '';
        maxVal = [];
        minVal = [];
        frameIncrement = 200;
        frameVector = 1:frameIncrement:nFrames;
        nFramesVector = length(frameVector);
        movieFilename = options.movieFilename;
        movieDatasetName = options.movieDatasetName;
        imgDim1 = imgs(1);
        imgDim2 = imgs(2);
        parfor frameNo = 1:nFramesVector
            if frameNo==nFramesVector
                thisFrame = readHDF5Subset(movieFilename,[0 0 frameVector(frameNo)-1],[imgDim1 imgDim2 nFrames-frameVector(frameNo)+1],'datasetName',movieDatasetName,'displayInfo',0); %NOTE: ADDED THE +1 TO ACCOUNT FOR MOVIES THAT GO ONE FRAME AVOVE A MULTIPLE OF 200 (19/01/2023)
            else
                thisFrame = readHDF5Subset(movieFilename,[0 0 frameVector(frameNo)-1],[imgDim1 imgDim2 frameIncrement],'datasetName',movieDatasetName,'displayInfo',0);
            end
            thisFrame(isnan(thisFrame)) = 0;
            maxVal(frameNo) = max(thisFrame(:));
            minVal(frameNo) = min(thisFrame(:));
            fprintf('Calculating fInc, done: %d/%d\n',frameNo,nFramesVector);
            % reverseStr = cmdWaitbar(frameNo,nFramesVector,reverseStr,'inputStr','calculating fInc','waitbarOn',1,'displayEvery',25);
        end
        maxVal = max(maxVal(:));
        minVal = min(minVal(:));
        fInc=(maxVal-minVal)/1000;
    end

    try
        if options.readMovieChunks==0
            testFrames=imgs(:,:,1:min(50,size(imgs,3)));
        else
            testFrames = loadMovieList(options.movieFilename,'frameList',1:min(50,nFrames),'inputDatasetName',options.movieDatasetName);
            testFrames(isnan(testFrames)) = 0;
        end
        %figure; hist(testFrames(:)); waitforbuttonpress()
        testFrames=round(testFrames(:)/fInc);
        startMu=mode(testFrames)*fInc;
        startSigma=0.004;
        fVals=(startMu-1000*fInc):fInc:(startMu+1000*fInc);
    catch err
        display(repmat('@',1,7))
        disp(getReport(err,'extended','hyperlinks','on'));
        display(repmat('@',1,7))

        fVals=0.95:fInc:1.2;
        startMu=[];
    end

    if options.readMovieChunks==0
        variation_in_time = std(imgs(:,:,framesForEst),[],3);
    else
        imgsTmp = loadMovieList(options.movieFilename,'frameList',framesForEst,'inputDatasetName',options.movieDatasetName);
        imgsTmp(isnan(imgsTmp)) = 0;
        variation_in_time = std(imgsTmp,[],3);
    end
    valid_xy = variation_in_time~=0; %Real data has some variation in time, black edges do not
else
    fInc=0.001;
    fVals=-0.2:fInc:0.2;
end

dist=zeros(size(fVals));
reverseStr = '';
nFramesEst = length(framesForEst(:));
for frameNo = 1:nFramesEst
    fr = framesForEst(frameNo);
    if options.readMovieChunks==0
        thisFrame=imgs(:,:,fr);
    else
        thisFrame = readHDF5Subset(options.movieFilename,[0 0 fr-1],[imgs(1) imgs(2) 1],'datasetName',options.movieDatasetName,'displayInfo',0);
        thisFrame(isnan(thisFrame)) = 0;
    end
    if ~options.vlm
        dist=dist+hist(thisFrame(valid_xy),fVals);
    else
        dist=dist+hist(thisFrame(:),fVals);
    end
    reverseStr = cmdWaitbar(frameNo,nFramesEst,reverseStr,'inputStr','calculating dist','waitbarOn',1,'displayEvery',25);
end

if ~options.vlm
    % then only take the bottom part of the distribution and reflect it
    % this gets rid of the heavy tail at the top, but it does tend to
    % underestimate the noise std dev
    if isempty(startMu)
        [~,zeroInd]=min(abs(fVals-1));
    else
        [~,zeroInd]=min(abs(fVals-startMu));
    end
    dist=[dist(1:zeroInd), fliplr(dist(1:zeroInd-1))];
    fVals=fVals(1:length(dist));
else
    [~,zeroInd]=min(abs(fVals));
    dist(zeroInd)=0;
    dist=dist(2:end-1);
    fVals=fVals(2:end-1);
end

dist=dist/(sum(dist(:))*fInc);

if ~options.vlm
    if isempty(startMu)
        startMu=1;
        startSigma=0.004;
    end
else
    startMu=0;
    startSigma=0.05;
end
try
    [noiseSigma,movieMean]=gaussfit(fVals,dist,startSigma,startMu);
catch
    noiseSigma=0.8*std(fVals);
    if ~isempty(startMu)
        movieMean=mode(fVals);
    else
        movieMean=startMu;
    end
end


% from file fitNoiseSigma on 12/11/13 at 5:59pm
