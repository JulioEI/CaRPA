function [cEstNew, phiEstNew, toDelete, scaledPhi, options] = EM_genFilt_oneIteration(imgs, cEst, phiEst, varargin) %noiseSigma, useConstantBG, numSigmasThresh)

    cEst=single(cEst);
    phiEst=single(phiEst);

    %%%%%%%%%% Set Options %%%%%%%%%%%%
    options = getDefaultCELLMaxOptions();

    %%%% replace default options with input options
    options=getOptions(options,varargin);

    nFrames=size(imgs,3);
    imgSize=size(imgs(:,:,1));
    nCells=size(cEst,1);
    nPixX=imgSize(2);
    nPixY=imgSize(1);
    nPix=nPixX*nPixY;

    thresh=options.numSigmasThresh*options.noiseSigma;
    fInc=options.noiseSigma/options.numPhotonsPerSigma;

    phiEstNew=phiEst;

    pXgivenZ=cEst;

    totCounts=zeros(1,nFrames,'single');
    labelProbabilities=zeros(nPix,nCells,'single'); % same size as pZgivenX

    %pX=pXgivenZ*phiEst;
    pX=phiEst'*pXgivenZ;
    pX(pX==0)=10e-8;
    pX=pX';
    pXgivenZ=pXgivenZ';

    %thisTimeLabelProbs=zeros(nPix,nCells,'single'); % same size as pZgivenX

    if options.selectRandomFrames && options.numFramesRandom<nFrames && ~options.iterOptions.lastIter
        % determine whether to use a novel seed or a predetermined one (for deterministic output)
        framesToProcess = options.subsampleFrameVector;
        % switch options.generateNovelSeed
        %     case 0
        %         rng(options.randNumGenSeed);
        %     case 1
        %         % nothing
        %     otherwise
        %         % body
        % end
        % framesToProcess=randperm(nFrames, options.numFramesRandom);
    else
        framesToProcess=1:nFrames;
    end
    if size(framesToProcess,1)>1
        framesToProcess=framesToProcess';
    end

    for t=framesToProcess

        % get photon counts
        counts=floor((imgs(:,:,t)-thresh)/fInc);
        counts(imgs(:,:,t)<thresh)=0;
        counts=counts(:);
        totCounts(t)=sum(counts);

        % update estimates from counts
        nonZeroInds=counts>0;
        if sum(nonZeroInds)>0

            countsOverPx=counts./pX(:,t);
            nonInfInds=~isinf(countsOverPx);
            countsOverPx(isinf(countsOverPx))=0;
            nonZeroInds=and(nonZeroInds,nonInfInds);

            if sum(nonZeroInds)>0
    %             %%% version with non-zero inds indexing:
    %              thisTimeLabelProbs=...
    %                  bsxfun(@times,pXgivenZ(nonZeroInds,:),countsOverPx(nonZeroInds)); % these are the indices where neither counts nor pX is zero
    %              thisTimeLabelProbs=...
    %                  bsxfun(@times, thisTimeLabelProbs, phiEst(:,t)'); % p(z|x) =
    %
    %             phiEstNew(:,t)=sum(thisTimeLabelProbs,1); % don't divide by total number of counts yet, since need to divide mu est by phi est before
    %
    %             %%%% reminder for size: labelProbabilities=zeros(nPix,nCells); % same size as pZgivenX
    %             labelProbabilities(nonZeroInds,:)=labelProbabilities(nonZeroInds,:)+thisTimeLabelProbs;



                %%% version WITHOUT non-zero inds indexing:
                  thisTimeLabelProbs=...
                      bsxfun(@times,pXgivenZ,countsOverPx); % these are the indices where neither counts nor pX is zero
                  thisTimeLabelProbs=...
                      bsxfun(@times, thisTimeLabelProbs, phiEst(:,t)'); % p(z|x) =


                phiEstNew(:,t)=sum(thisTimeLabelProbs,1); % don't divide by total number of counts yet, since need to divide mu est by phi est before

                %%%% reminder for size: labelProbabilities=zeros(nPix,nCells); % same size as pZgivenX
                labelProbabilities=labelProbabilities+thisTimeLabelProbs;
            end
        end
    end



    %%%%% calc cEst by dividing labelProbabilities by the total weight for that
    %%%%% cell
    nonZeroInds=sum(phiEstNew,2)>0;
    cEstNew=cEst;
    cEstNew(nonZeroInds,:)=bsxfun(@rdivide,labelProbabilities(:,nonZeroInds)',sum(phiEstNew(nonZeroInds,:),2));
    if sum(isnan(cEstNew(:)))>0
        error('nans!!!!!!!!!!!!!!!!!!!!!')
    end


    scaledPhi=phiEstNew*fInc;
    for t=1:nFrames
        if totCounts(t)>0
            phiEstNew(:,t)=phiEstNew(:,t)/totCounts(t);
        end
    end

    toDelete=zeros(1,nCells);
    if options.elimOnScaledPhi
        for cInd=1:nCells
            if max(scaledPhi(cInd,:))<options.threshForElim
                toDelete(cInd)=1;
            end
        end
    else
        for cInd=1:nCells
            if max(phiEstNew(cInd,:))<options.threshForElim
                toDelete(cInd)=1;
            end
        end
    end

    if sum(toDelete)==nCells
        toDelete(end)=0;
    end

    if options.useConstantBG && ~isempty(cEstNew)
        cEstNew(end,:)=1/(size(cEstNew,2))*ones(1,size(cEstNew,2));
        toDelete(end)=0;
    end

    toDelete=logical(toDelete);
end