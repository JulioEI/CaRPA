function [paramStore, changeParams, options, phiEst] = EM_genFilt(imgs, varargin)


    %%%%%%%%%% Set Options %%%%%%%%%%%%
    options = getDefaultCELLMaxOptions();

    %%%% Replace default with input options
    options=getOptions(options, varargin);

    % set progress square settings
    if isfield(options, 'progressSquare')
        progressSquare=options.progressSquare;
    else
        progressSquare=[];
    end

    % initialize.
    [initImages, initProbs] = initializeEM(imgs, 'options', options);
    nCells=size(initProbs,1);
    nFrames=size(imgs,3);

    % if imgs centered on 1, subtract 1 to center on 0
    % note that this breaks for other movie normalizations
    if options.oneCentered
        imgs=imgs-1;
    end
    imgs(imgs<0)=0;


    % iterate through EM cycles
    options.threshForElim=min(1/(3*nCells), options.threshForElim);
    numIters=1;
    if ~isempty(progressSquare)
        figNo = 777;
        if ishandle(figNo)
            set(0,'CurrentFigure',figNo)
            figHandle = figNo;
        else
            figHandle = figure(figNo);
        end
        % figure(777);
        hold off; imagesc(progressSquare); colormap(gray); drawnow;
    end
    paramStore=cell(options.maxIters+1,1);
    changeParams=zeros(options.maxIters+1,1);
    cEst=initImages;
    phiEst=initProbs;
    muEst = getImageCentroids(cEst,size(imgs(:,:,1)));%*********************************************
    if options.useScaledPhi
        paramStore{numIters}=[muEst, cEst, zeros(size(phiEst))];%*********************************************
    else
        paramStore{numIters}=[muEst, cEst];%*********************************************
    end

    % initialize param converge movie if asked for
    if options.saveParamConvergeMovie
        writerObj=initAVIwriter(options.paramConvergePathname, 4);
        paramxlims=[1 size(imgs,2)]; paramylims=[1 size(imgs,1)];
        saveParamConvergeFrame(muEst, options.trueCentroids, options.doppCentroids, writerObj, paramxlims, paramylims)
        saveParamConvergeFrame(muEst, options.trueCentroids, options.doppCentroids, writerObj, paramxlims, paramylims)
        saveParamConvergeFrame(muEst, options.trueCentroids, options.doppCentroids, writerObj, paramxlims, paramylims)
    end

    % set iteration options and begin iterating to convergence
    options.iterOptions.notConverged=1;
    options.iterOptions.lastIter=0;
    while options.iterOptions.notConverged

    %     cellImages=reshape(cEst,[size(cEst,1), size(imgs(:,:,1))]);
    %     cellImages=permute(cellImages, [2 3 1]);
    %     figure(876); imagesc(max(cellImages(:,:,1:end-1),[],3));
    %     waitforbuttonpress
    %     for cInd=1:size(cellImages,3)
    %         imagesc(cellImages(:,:,cInd))
    %         waitforbuttonpress;
    %     end

        if ~isempty(options.subsampleFrameMatrix)
            options.subsampleFrameVector = find(options.subsampleFrameMatrix(numIters,:));
        end

        % do one EM iteration
        [cEstNew, phiEstNew, toDelete, scaledPhi, options] = EM_genFilt_oneIteration(imgs, cEst, phiEst,'options', options);

        % calculate the change in images/centroids as specified
        muEstNew = getImageCentroids(cEstNew,size(imgs(:,:,1)));%*********************************************
        if options.useMuChangeOnly
            muEst = getImageCentroids(cEst,size(imgs(:,:,1)));%*********************************************
            thisChangeParams=max(abs(muEst(:)-muEstNew(:)));
        else
            goodInds=zeros(size(cEstNew));
            for cInd=1:size(cEstNew,1)
                goodInds(cInd,:)=cEstNew(cInd,:)>0.01*max(cEstNew(cInd,:));
            end
            goodInds=logical(goodInds);
            thisChangeParams=mean(abs(cEst(goodInds)-cEstNew(goodInds))./cEst(goodInds));
        end
        thisChangeParams(isnan(thisChangeParams))=1000;


        % set(0,'DefaultFigureWindowStyle' , 'normal')
        % openFigure(numIters);
            % subplot(options.maxIters,3,(numIters-1)*3+1);
            % subplot(1,3,1)
                % tmpIterImg = initImages(1:end-1,:);
                % tmpIterImg=reshape(tmpIterImg,[size(tmpIterImg,1), size(imgs(:,:,1))]);
                % tmpIterImg=permute(tmpIterImg, [2 3 1]);
                % playMovie(tmpIterImg);
                % imagesc(nanmax(tmpIterImg,[],3));
                % title('initial grid')
                % iterationImages.init = nanmax(tmpIterImg,[],3);
            % subplot(options.maxIters,3,(numIters-1)*3+2);
            % subplot(1,3,2)
                tmpIterImg=reshape(cEst,[size(cEst,1), size(imgs(:,:,1))]);
                tmpIterImg=permute(tmpIterImg, [2 3 1]);
                % imagesc(nanmax(tmpIterImg,[],3));
                % title(sprintf('iteration #%d',numIters))
                iterationImages.iterImg{numIters} = nanmax(tmpIterImg,[],3);
            % subplot(options.maxIters,3,3)
            % subplot(1,3,3)
                % imagesc(nanmax(imgs,[],3))
                % title('movie max')
                % drawnow
                % iterationImages.maxMovie = nanmax(tmpIterImg,[],3);
                % save('cellmax_iter_images.mat','iterationImages','-v7.3');

        % remove cells with scaledPhi that is too correlated
        if options.removeCorrProbs && ((mod(numIters-1, 5)==0 && numIters-1>min(15, options.minIters-1)) || ...
                ~((thisChangeParams>options.maxDeltaParams || numIters<options.minIters) && numIters<options.maxIters+1))
            C=corr(scaledPhi');
            for cInd1=1:size(scaledPhi,1)
                for cInd2=cInd1+1:size(scaledPhi,1)
                    if  C(cInd1,cInd2)>options.scaledPhiCorrThresh && ...
                            norm(muEstNew(cInd1,:)-muEstNew(cInd2,:))<options.distanceThresh && ...
                            ~toDelete(cInd2) && ~toDelete(cInd1)
                        image1=reshape(cEst(cInd1,:),size(imgs(:,:,1)));
                        image2=reshape(cEst(cInd2,:),size(imgs(:,:,1)));
                        areaOverlap = calcAreaOverlap(image1, image2);
                        if areaOverlap>options.corrRemovalAreaOverlapThresh
                            cEstNew(cInd1,:)=(cEstNew(cInd1,:)+cEstNew(cInd2,:))/2;
                            scaledPhi(cInd1,:)=(scaledPhi(cInd1,:)+scaledPhi(cInd2,:))/2;
                            toDelete(cInd2)=1;
        %                     %%% To visually inspect which images are getting
        %                       combined at this step...
        %                     filtTraces=calculateFilteredTraces(imgs, cEst([cInd1,cInd2],:)', 'options', options);
        %                     eventTimes=detectEventsOnPhi(double(filtTraces), double(scaledPhi([cInd1,cInd2],:)));
        %                     [v,~] = manualCellChecker(imgs, scaledPhi([cInd1,cInd2],:), single(reshape(cEst([cInd1,cInd2],:)', [size(imgs(:,:,1)) 2])),...
        %                         eventTimes);
                        end
                    end
                end
            end
        end

        % delete cells that were marked for deletion in the iteration or
        %   scaledPhi corr removal
        phiEstNew(toDelete,:)=[];
        cEstNew(toDelete,:)=[]; %*********************************************
        muEstNew(toDelete,:)=[];
        scaledPhi(toDelete,:)=[];
        % ignore first several iterations for removing cells to avoid eliminating good cells too quickly
        % if options.selectRandomFrames && options.numFramesRandom<nFrames && ~options.iterOptions.lastIter
        %     if ~isempty(options.subsampleFrameMatrix)
        %         skipIterDelete = round(1/nanmean(options.subsampleFrameMatrix(numIters,:)));
        %     else
        %         skipIterDelete = round(1/(options.numFramesRandom/nFrames));
        %     end
        %     if numIters>skipIterDelete
        %         phiEstNew(toDelete,:)=[];
        %         cEstNew(toDelete,:)=[]; %*********************************************
        %         muEstNew(toDelete,:)=[];
        %         scaledPhi(toDelete,:)=[];
        %     else
        %         fprintf('iter %d, ignore delete until %d\n',numIters,skipIterDelete);
        %     end
        % else
        %     phiEstNew(toDelete,:)=[];
        %     cEstNew(toDelete,:)=[]; %*********************************************
        %     muEstNew(toDelete,:)=[];
        %     scaledPhi(toDelete,:)=[];
        % end

        % write to param converge movie if asked for
        if options.saveParamConvergeMovie
            saveParamConvergeFrame(muEstNew, options.trueCentroids, options.doppCentroids, writerObj, paramxlims, paramylims)
        end

        % store value for comparison next iteration
        cEst=cEstNew;
        phiEst=phiEstNew;

        % store change params, params etc.
        changeParams(numIters)=thisChangeParams;
        numIters=numIters+1;
        if options.useScaledPhi
            paramStore{numIters}=[muEstNew, cEst, scaledPhi];%*********************************************
        else
            paramStore{numIters}=[muEstNew, cEst];%*********************************************
        end


         % if progress square is input, display it, with current convergence
         if ~isempty(progressSquare)
           xMin=1;
           xMax=options.maxIters;
           xRange=xMax-xMin;
           xMinCenter=xMin+xRange/(2*size(progressSquare,2));
           xMaxCenter=xMax-xRange/(2*size(progressSquare,2));
           if xMax==xMin
               xMax=xMin+1;
           end

           yMin=log10(min(options.maxDeltaParams, min(changeParams(1:numIters-1))));
           yMin(or(isnan(yMin),isinf(yMin)))=log10(options.maxDeltaParams);
           yMax=max(log10(changeParams(1:numIters-1)));
           yMax(or(isnan(yMax),isinf(yMax)))=max(0,yMin+1);
           yRange=yMax-yMin;
           yMinCenter=yMin+yRange/(2*size(progressSquare,1));
           yMaxCenter=yMax-yRange/(2*size(progressSquare,1));
           if yMax<=yMin
               yMax=yMin+1;
           end

           plot(log10(changeParams(1:numIters-1))); plot(log10(changeParams(1:numIters-1)), '.'); hold on;
           imagesc(linspace(xMinCenter, xMaxCenter, size(progressSquare,2)), linspace(yMinCenter, yMaxCenter, size(progressSquare,1)), flipud(progressSquare));
           xlim([xMin xMax]); ylim([yMin yMax]);
           colormap(gray); set(gca, 'CLim', [0 1]);
           hold on; plot(log10(changeParams(1:numIters-1))); plot(log10(changeParams(1:numIters-1)), '.');
           title('Progress in movie space...'); drawnow; hold off
         elseif options.plotDeltaParams
             hold off; plot(log10(changeParams(1:numIters-1))); ...
                 hold on; plot(log10(changeParams(1:numIters-1)), '.'); ...
                 ylabel('log mean change in params'); xlabel('Iteration'); drawnow
         end

         paramsChangeIndicator = thisChangeParams<=options.maxDeltaParams; % if this is 1, converged
         minItersIndicator = numIters>options.minIters; % if this is 1, reached minimum number of iterations
         maxItersIndcator = numIters>options.maxIters; % if this is 1, forced to finish bc of max iters
         changeChangeIndicator = mean(diff(changeParams(max(numIters-10,1):numIters-1)))<=options.maxDeltaDeltaParams;

         if maxItersIndcator || ((paramsChangeIndicator || changeChangeIndicator) && minItersIndicator)
             if options.iterOptions.lastIter
                options.iterOptions.notConverged=0;
                options.iterOptions.lastIter=0;
             else
                options.iterOptions.lastIter=1;
             end
         end
    end
    if numIters<length(paramStore)
        paramStore=paramStore(1:numIters);
        changeParams=changeParams(1:numIters);
    end

    if options.saveParamConvergeMovie
        close(writerObj)
    end

    %
    % tmpIterImg = initImages(1:end-1,:);
    % tmpIterImg=reshape(tmpIterImg,[size(tmpIterImg,1), size(imgs(:,:,1))]);
    % tmpIterImg=permute(tmpIterImg, [2 3 1]);
    % iterationImages.init = nanmax(tmpIterImg,[],3);
    % iterationImages.maxMovie = nanmax(tmpIterImg,[],3);
    % currentDateTimeStr = datestr(now,'yyyy_mm_dd_HHMMSS','local');
    % save(['cellmax_iter_images_' currentDateTimeStr '.mat'],'iterationImages','-v7.3');
end