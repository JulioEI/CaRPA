function [position,velocity,score] = getMouseTrajectory(moviePath,varargin)% This is the original written by David
% function [position,velocity,score] = getMouseTrajectory(Thresh,moviePath,varargin)%% This one I wrote to check whci backgroundTreshold works better
%Load movie
mov = VideoReader(moviePath);

%Process options
options.bgLength = ceil(mov.NumberOfFrames/10);
options.mouseColorT = 0.2;
options.nightMode = false;
options.backgroundTreshold = -0.8; % David's Original (Works for Global Remapping)
% options.backgroundTreshold = -0.9 % Thresh; % Pasha (Better for Novel Object Locations)
options.smoothingWindowSize = 10;
options.lonlyValLenTresh = 5;
options.doPlot = 0;
options.verbose = 1;

for k = 1:2:length(varargin)
    options.(varargin{k}) = varargin{k+1};
end  

%Compute background
if options.verbose;disp('Computing background...');end

if ~options.nightMode
    movieBG = zeros([mov.Height,mov.Width,options.bgLength]);
    randIdx = randperm(mov.NumberOfFrames);
    for k = 1:options.bgLength
        bgFrame =  rgb2gray(double(mov.read(randIdx(k)))/255);
        bgFrame(bgFrame>options.mouseColorT) = 1;
        movieBG(:,:,k) = bgFrame;
    end
    bg = median(movieBG,3);
else
    bg = zeros([mov.Height,mov.Width,options.bgLength]);
end

if options.verbose;disp('Computing trajectories...');end

traceLog = nan([mov.NumberOfFrames,2]);
parfor k = 1:mov.NumberOfFrames
    %Threshold movie
    frame =  rgb2gray(double(mov.read(k))/255);
    if options.nightMode
        frameBin = frame>.5*max(frame(:));
    else    
        frame(frame>options.mouseColorT) = 1;
        %Substract background and trehshold again
        frameNoBg = frame - bg;
        frameNoBg(frameNoBg>options.backgroundTreshold) = 0;
        frameNoBg = abs(frameNoBg);
        frameBin = double(frameNoBg>0);
        %Open and erode bin movie
        frameBin = bwareaopen(frameBin,10); %get rid of small regions
        frameBin = bwmorph(bwmorph(frameBin,'erode'),'erode'); %erode the binary image to clean up cable
    end
    %Find centroid of largest object
    regProps = regionprops(frameBin,'Centroid','Area');
    [~,idx] = max([regProps.Area]);
    if ~isempty(regProps)
        traceLog(k,:) = regProps(idx).Centroid;
    end    
end

%%
%Get velocity, remove outliers
if options.verbose;disp('Processing trajectories and computing velocities...');end

[position,velocity,score] = postProcessTracesTracker(traceLog,options.smoothingWindowSize,options.lonlyValLenTresh,options.doPlot);

%%
if options.doPlot
    plotTracesOverMovie(position,velocity,mov)
disp('-----------------------------------------------------------------------------')
disp(['Current threshold = ',num2str(Thresh)])
disp('-----------------------------------------------------------------------------')
end



