function [imgs,cellTraces,params]=simulateData_FromSpikeTimes(spikeTimes, amps, taus, totTime,...
    params, imgSize, noiseSigma, varargin)

% spikeTimes: nCells x 1 cell, each element 1 x nSpikes_cInd
% amps: nCells x 1 cell, each element 1 x nSpikes_cInd
% taus
theNoise=[];
C=[];
displayTraces=0;
displayImgs=0;
maxJitter=0;
extraBorder=0;
jitterVec=[];
if ~isempty(varargin)
    options=varargin{1};
    if isfield(options, 'C')
        C=options.C;
    end
    if isfield(options, 'theNoise')
        theNoise=options.theNoise;
    end
    if isfield(options, 'displayImgs')
        displayImgs=options.displayImgs;
    end
    if isfield(options, 'maxJitter')
        maxJitter=options.maxJitter;
    end
    if isfield(options, 'extraBorder')
        extraBorder=options.extraBorder;
    end
    if isfield(options, 'jitterVec')
        jitterVec=options.jitterVec;
        maxJitter=100;
    end
end

nCells=length(spikeTimes);
cellTraces=zeros(nCells,totTime);
for cInd=1:nCells
    theseSpikeTimes=spikeTimes{cInd};
    theseAmps=amps{cInd};
    if length(theseAmps)~=length(theseSpikeTimes)
        error('Amplitudes and spike times must be the same length')
    end
    for sInd=1:length(theseSpikeTimes)
        sTime=theseSpikeTimes(sInd);
        ampl=theseAmps(sInd);
        tInds=sTime:totTime;
        cellTraces(cInd,tInds)=cellTraces(cInd,tInds)+...
            ampl*exp(-(tInds-sTime)/taus(cInd));
        if displayTraces
            figure(13); plot(1:length(cellTraces(cInd,:)), cellTraces(cInd,:))
            pause(0.5)
        end
    end
end

if maxJitter>0
    if isempty(jitterVec)
        jitterSigma=maxJitter/3;
        jitter=randn(2,totTime)*jitterSigma^2;
        jitter(abs(jitter)>maxJitter)=maxJitter;
        padAmount=ceil(maxJitter);
    else
        jitter=jitterVec';
        if size(jitter,1)>2
            jitter=jitter';
        end
        padAmount=ceil(max(jitterVec(:)));
    end
    imgSize=imgSize+2*padAmount;
    params(:,1:2)=params(:,1:2)+padAmount;
end
if extraBorder>0
    disp('putting extra space!')
    params(:,1:2)=params(:,1:2)+extraBorder;
    imgSize=imgSize+2*extraBorder;
end

cellImgs=calcCellImgs(params, imgSize);
cellImgs=reshape(cellImgs, [imgSize(1)*imgSize(2), nCells]);

imgs=cellImgs*cellTraces;
imgs=reshape(imgs, [imgSize totTime])+1;


if isempty(theNoise) && isempty(C)
    imgs=imgs+noiseSigma*randn(size(imgs));
elseif ~isempty(theNoise)
    for fr=1:size(imgs,3)
        imgs(:,:,fr)=imgs(:,:,fr)+theNoise(:,:,fr);
    end
else
    for fr=1:size(imgs,3)
        imgs(:,:,fr)=imgs(:,:,fr)+reshape(mvnrnd(zeros(1,size(C,1)),C),size(imgs(:,:,1)));
    end
end


if maxJitter>0
    [xTarg,yTarg]=meshgrid(1:imgSize(2),1:imgSize(1));
    for fr=1:size(imgs,3)
        xSrc=xTarg-jitter(1,fr);
        ySrc=yTarg-jitter(2,fr);
        imgs(:,:,fr)=interp2(xTarg,yTarg,imgs(:,:,fr),xSrc,ySrc);
    end
    imgs=imgs((padAmount+1):(end-padAmount),(padAmount+1):(end-padAmount),:);
    params(:,1:2)=params(:,1:2)-padAmount;
end

    
if displayImgs
    % colors for plotting centroids and traces
    colors=[1 0 0;
            0 1 0;
            0 0 0;
            0 1 1;
            1 0 1;
            1 1 0];
    nColors=size(colors,1);

    for t=1:totTime
        figure(14); imagesc(imgs(:,:,t))
        set(gca, 'CLim', [0.98, 1.08])
        hold all
        for cInd=1:nCells
            plot(params(cInd,1),params(cInd,2),'*',...
                'Color', colors(mod(cInd,nColors)+1,:))
        end
        hold off
        
        pause(1/19.3)
    end
end