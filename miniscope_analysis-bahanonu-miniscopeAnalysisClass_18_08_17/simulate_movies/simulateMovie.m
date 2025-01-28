function [imgs,realCellTraces,params,spikeTimes,options]=simulateMovie(nCellsDesired,...
    imgSize, totTime, varargin)

%%%% INPUTS
% nCellsDesired : number of cells to simulate
% imgSize : size of frame in pixels, ie [250 300]
% totTime : total time of movie to simulate, in frames

%%%% OPTIONS
options.minDist=8;  % minimum distance between cell centroids, in pixels. cells are set to be 2D gaussians with std dev of 2-3 pixels (so each cell is about 10 pixels wide)
options.eventRate=1/100;    % event rate in events/frame (timings create 20hz movie)
options.noiseSigma=0.008;
options.minSpikeAmplitude=0.032;
options.spikeAmplitudeAddlRange=0.08;
options.minTimeConstant=7;
options.timeConstAddlRange=8;
options=getOptions(options,varargin);


% position the cells randomly
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
    
    if min(dists(:))>options.minDist
        xpos(nCells+1)=candX;
        ypos(nCells+1)=candY;
        nCells=nCells+1;
    end
end
nCells=nCells-1;
xpos=xpos(1:nCells);
ypos=ypos(1:nCells);


% set cell sizes randomly
params=zeros(nCells,5);
params(:,1)=xpos;
params(:,2)=ypos;
params(:,3)=2+rand(nCells,1);
params(:,4)=2+rand(nCells,1);
params(:,5)=pi/2*rand(nCells,1);


% assign spike times, amplitudes, and decay times randomly
spikeTimes=cell(1,nCells);
amps=cell(1,nCells);
for i=1:nCells
    spikeYesNo=rand([1,totTime]);       % generate poisson spikes (yes/no) for each bin
    spikeTimes{i}=find(spikeYesNo<=options.eventRate);
    amps{i}=options.minSpikeAmplitude+options.spikeAmplitudeAddlRange*rand(1,length(spikeTimes{i}));
end
taus=options.minTimeConstant+options.timeConstAddlRange*rand(1,nCells);

[imgs,realCellTraces,params]=simulateData_FromSpikeTimes(spikeTimes,amps,taus,totTime,params,imgSize,options.noiseSigma,options);