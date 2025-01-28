function allCellTraces = calculateTraces(cellImgs, imgs, varargin)

    % Written by Lacey Kitch in 2013
    % 2017.01.18 - updated to include support for movie chunks

    % haveBG=0;
    % if ~isempty(varargin)
    % 	options=varargin{1};
    %     if isfield(options, 'suppressOutput')
    %         suppressOutput=options.suppressOutput;
    %     else
    %         suppressOutput=0;
    %     end
    %     if isfield(options, 'bg')
    %         haveBG=1;
    %         bg=options.bg;
    %     end
    % else
    % 	suppressOutput=0;
    % end

    options.suppressOutput = 0;
    options.bg = [];
    options.removeBelowThreshPixelsForRecalc=0;
    options.numSigmasThresh=0;
    options.noiseSigma=0.004;

    options.datasetName = '/1';
    options.displayInfo = 1;
    options.movieFilename = '';
    options.movieDatasetName = '/Data/Images';
    options.readMovieChunks = 0;
    options.movieMaxVal = [];

    options=getOptions(options,varargin);

    if isempty(options.bg)
        haveBG = 0;
    else
        haveBG = 1;
        bg = options.bg;
    end

    if options.readMovieChunks==0
        if numel(size(cellImgs))==3
            % if cellImgs and imgs are still 3D and have not been reshaped
            nCells=size(cellImgs,3);
            imgSize=size(imgs(:,:,1));
            cellImgs=reshape(cellImgs,[imgSize(1)*imgSize(2),nCells]);
        else
            % if cellImgs and imgs are 2D and so have been reshaped
            nCells=size(cellImgs,2);
        end

        if numel(size(imgs))==2
            nFrames=size(imgs,2);
            reshapeImgs=0;
        else
            imgSize=size(imgs(:,:,1));
            nFrames=size(imgs,3);
            reshapeImgs=1;
        end
    else
        if numel(size(cellImgs))==3
            % if cellImgs and imgs are still 3D and have not been reshaped
            nCells=size(cellImgs,3);
            imgSize=[imgs(1) imgs(2)];
            cellImgs=reshape(cellImgs,[imgSize(1)*imgSize(2),nCells]);
        else
            % if cellImgs and imgs are 2D and so have been reshaped
            nCells=size(cellImgs,2);
        end

        imgSize=[imgs(1) imgs(2)];
        nFrames=imgs(3);
        reshapeImgs=1;

        % if numel(size(imgs))==2
        %     nFrames=size(imgs,2);
        %     reshapeImgs=0;
        % else
        %     imgSize=size(imgs(:,:,1));
        %     nFrames=size(imgs,3);
        %     reshapeImgs=1;
        % end
    end


    if haveBG
        if numel(size(bg))==2
            bg=reshape(bg, [imgSize(1)*imgSize(2) 1]);
        end
    end

    % check for repetitive cellImgs, which make cellImgs singular, and
    % eliminate them from the computation
    [~,s,v]=svd(cellImgs'*cellImgs);
    s=diag(s);
    lowSVs=find(s<(10*eps*max(s)))';
    goodCellInds=1:nCells;
    badCellInds=zeros(1,nCells);
    for lowSVind=lowSVs
        [~,badCellInd]=max(v(:,lowSVind));
        badCellInds(badCellInd)=1;
        v(badCellInd,:)=-10000;
    end
    if ~isempty(lowSVs) && ~options.suppressOutput
        disp(['Found ' num2str(length(lowSVs)) ' redundant cell images']);
    end
    goodCellInds(logical(badCellInds))=[];
    cellImgs(:,logical(badCellInds))=[];

    if isempty(options.movieMaxVal)
        options.movieMaxVal = max(imgs(:));
    end

    if options.movieMaxVal>1
        subOne=1;
    else
        subOne=0;
    end

    % cycle through chunks of time and of cells (for RAM preservation) and calculate traces
    allCellTraces=zeros(nCells, nFrames);
    cImgPrefactor=(cellImgs'*cellImgs)\cellImgs';
    clear cellImgs
    reverseStr = '';
    for fr=1:200:nFrames

        if mod(fr,200)==0
            reverseStr = cmdWaitbar(fr,nFrames,reverseStr,'inputStr','calculating traces frame-by-frame','waitbarOn',1,'displayEvery',200);
        end

        tLims=fr:min(fr+199,nFrames);
        if reshapeImgs
            if options.readMovieChunks==0
                tLimsImgs = imgs(:,:,tLims);
            else
                tLimsImgs = loadMovieList(options.movieFilename,'frameList',tLims,'inputDatasetName',options.movieDatasetName,'displayInfo',0);
                tLimsImgs(isnan(tLimsImgs)) = 0;
            end
            if haveBG
                thisImgs=reshape(tLimsImgs,[imgSize(1)*imgSize(2),length(tLims)])-...
                    repmat(bg, [1 length(tLims)]);
            else
                thisImgs=reshape(tLimsImgs,[imgSize(1)*imgSize(2),length(tLims)]);
                if subOne
                    thisImgs=thisImgs-1;
                end
            end
        else
            if haveBG
                thisImgs=imgs(:,tLims)-repmat(bg, [1 length(tLims)]);
            elseif max(imgs(:,tLims))>1
                thisImgs=imgs(:,tLims)-1;
            else
                thisImgs=imgs(:,tLims);
            end
        end
        if options.removeBelowThreshPixelsForRecalc
            thisImgs(thisImgs<options.numSigmasThresh*options.noiseSigma)=0;
        end
        allCellTraces(goodCellInds,tLims)=cImgPrefactor*thisImgs;
    end
end