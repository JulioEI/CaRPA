function [initImages, initProbs] = initializeEM(imgs, varargin)

options=getDefaultCELLMaxOptions();
options=getOptions(options,varargin);

initImages=[];
initProbs=[];
switch lower(options.initMethod)
    case 'input'
%         initImages=options.localInitImages;
%         initImages=permute(initImages, [3 1 2]);
%         initImages=reshape(initImages, [size(initImages,1), numel(imgs(:,:,1))]);
        [initImages, ~] = initImagesFromICs(options.localInitImages, options.localInitTraces, 'options', options);
    case 'ica'
        if isfield(options, 'localICimgs')
            [initImages, ~] = initImagesFromICs(options.localICimgs, options.localICtraces, 'options', options);%*********************************************
        else
            disp('Warning: to initialize with ICA results, input IC images. Initializing with grid...')
        end
    case 'vertstripes'
        initImages=initPurkinjeParams(imgs, 4);%*********************************************
    case 'random'
        initParams=gridInitialize(imgs, options.gridSpacing);
        initImages=rand(size(initParams,1),numel(imgs(:,:,1)));%*********************************************
        initImages=bsxfun(@rdivide,initImages,sum(initImages,2));%*********************************************
    case 'spikes'
        nBlobs=size(imgs,1);
        initParams=[round(size(imgs,2)/2)*ones(nBlobs,1),(1:nBlobs)',50*ones(nBlobs,1),10*ones(nBlobs,1),zeros(nBlobs,1)];
        initImages=calcCellImgs(initParams,size(imgs(:,:,1)));
        initImages=permute(initImages, [3 1 2]);
        initImages=reshape(initImages, [size(initParams,1), numel(imgs(:,:,1))]);%*********************************************
        initImages=bsxfun(@rdivide, initImages, sum(initImages,2));
end
if isempty(initImages)
    [initParams, ~]=gridInitialize(imgs, 'gridSpacing',options.gridSpacing,'gridWidth',options.gridWidth);
    initImages=calcCellImgs(initParams,size(imgs(:,:,1)));
    % figure(101010101);imagesc(nanmax(initImages,[],3));drawnow
    initImages=permute(initImages, [3 1 2]); %*********************************************
    initImages=reshape(initImages, [size(initParams,1), numel(imgs(:,:,1))]);%*********************************************
    % normalize, but causes edge images to be given greater weight
    initImages=bsxfun(@rdivide, initImages, sum(initImages,2));
end

nCells=size(initImages,1);%*********************************************
if options.useConstantBG
    initImages=[initImages; 1/(size(initImages,2))*ones(1,size(initImages,2))];%*********************************************
    nCells=nCells+1;
end
% tmpIterImg = initImages(1:end-1,:);
% tmpIterImg=reshape(tmpIterImg,[size(tmpIterImg,1), size(imgs(:,:,1))]);
% tmpIterImg=permute(tmpIterImg, [2 3 1]);
% % playMovie(tmpIterImg);
% figure(101010102);imagesc(nanmax(tmpIterImg,[],3));drawnow

if isempty(initProbs)
    % initialize all cells as equally probable on all frames at first
    % this could be changed to use relative brightness etc.
    initProbs=(1/nCells)*ones(nCells,size(imgs,3));
    if options.initProbsWithFiltTraces && size(initProbs,1)>1
        tempImages=reshape(initImages', [size(imgs(:,:,1)) nCells]);
        initProbs(1:nCells-1,:)=calculateFilteredTraces(imgs,tempImages(:,:,1:nCells-1),'oneCentered',options.oneCentered); clear tempImages
        initProbs(initProbs<0)=0;
        initProbs(nCells,:)=mean(initProbs(1:nCells-1,:),1);
        nonZeroInds=sum(initProbs,1)>0;
        initProbs(:,nonZeroInds)=bsxfun(@rdivide,initProbs(:,nonZeroInds),sum(initProbs(:,nonZeroInds),1));
        clear nonZeroInds
    end
elseif size(initProbs,1)<size(initImages,1)%*********************************************
    probRowsToAdd=size(initImages,1)-size(initProbs,1);%*********************************************
    initProbs=[initProbs; 1/nCells*(ones(probRowsToAdd,size(initProbs,2)))];
end

% consider adding this here - or initializing the background with the mean
% of the other cells, or something similar
% oh, or maybe the background is already being included in the filt trace
% calculation so don't need this
% in that case we shouldn't even need the else clause
% nonZeroInds=sum(initProbs,1)>0;
% initProbs(:,nonZeroInds)=bsxfun(@rdivide,initProbs(:,nonZeroInds),sum(initProbs(:,nonZeroInds),1));
% clear nonZeroInds