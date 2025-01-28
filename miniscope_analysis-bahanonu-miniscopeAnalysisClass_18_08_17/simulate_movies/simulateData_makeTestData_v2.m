function [imgs,realCellTraces,params,spikeTimes,noiseSigma]=simulateData_makeTestData_v2(nCellsDesired,...
    imgSize, minDist, eventRate, totTime,varargin)

% nCellsDesired : number of cells to simulate
% imgSize : size of frame in pixels, ie [250 300]
% minDist : minimum distance, in pixels, to enforce between cells
% eventRate : mean event rate, in events per frame
% totTime : total time of movie to simulate, in frames

% example call :
% nCellsDesired=800;
% imgSize=[250 300];
% minDist=8;
% eventRate=1/100;  % in events/frame
% totTime=1000;

% example call for a small snippet of movie :
% [imgs,realCellTraces,params]=simulateData_makeTestData_v2(35, [40 40], 3, 6/1000, 500);


options.minAmpStdDevs = 4;
options.maxAmpStdDevs = 6;
options.minTimeConstFrames = 7;
options.maxTimeConstFrames = 15;
options.noiseSigma = 0.06;
options.minCellWidthPix = 2;
options.minCellWidthPix = 3;
options.theNoise=[];
options.C=[];
options.displayTraces=0;
options.displayImgs=0;
options.maxJitter=0;
options.extraBorder=0;
options.jitterVec=[];
options = getOptions(options,varargin);


% select random positions for cells, guaranteeing a minimum spacing
noiseSigma = options.noiseSigma;
nCells=1;
xpos=zeros(nCellsDesired,1);
ypos=zeros(nCellsDesired,1);
xpos(1)=rand*imgSize(2);
ypos(1)=rand*imgSize(1);

nLoops=0;
nLoopsMax=10^6;
while nCells<nCellsDesired && nLoops<nLoopsMax
    
    nLoops=nLoops+1;
    
    candX=rand*imgSize(2);
    candY=rand*imgSize(1); 

    dists=squareform(pdist([[xpos(1:nCells);candX] [ypos(1:nCells);candY]]))+...
        diag(100*ones(nCells+1,1));
    
    if min(dists(:))>minDist
        xpos(nCells+1)=candX;
        ypos(nCells+1)=candY;
        nCells=nCells+1;
    end
end
nCells=nCells-1;
xpos=xpos(1:nCells);
ypos=ypos(1:nCells);

% these parameters set the positions (x,y), std devs (x,y), and angle of
% gaussians. the cells will be gaussians with these params
params=zeros(nCells,5);
params(:,1)=xpos;
params(:,2)=ypos;
params(:,3)=options.minCellWidthPix+...
    (options.minCellWidthPix - options.minCellWidthPix)*rand(nCells,1);
params(:,4)=options.minCellWidthPix+...
    (options.minCellWidthPix - options.minCellWidthPix)*rand(nCells,1);
params(:,5)=pi/2*rand(nCells,1);


% generate poisson spike times and random amplitudes and decay constants,
% based on the specified inputs
spikeTimes=cell(1,nCells);
amps=cell(1,nCells);
for i=1:nCells
    spikeYesNo=rand([1,totTime]);       % generate poisson spikes (yes/no) for each bin
    spikeTimes{i}=find(spikeYesNo<=eventRate);
    amps{i}=options.minAmpStdDevs*noiseSigma+...
        (options.maxAmpStdDevs - options.minAmpStdDevs)*options.noiseSigma*rand(1,length(spikeTimes{i}));
end
taus=options.minTimeConstFrames+...
    (options.maxTimeConstFrames - options.minTimeConstFrames)*rand(1,nCells);

% simulate the movie
[imgs,realCellTraces,params]=simulateData_FromSpikeTimes(spikeTimes,amps,taus,...
    totTime,params,imgSize,noiseSigma,'options', options);