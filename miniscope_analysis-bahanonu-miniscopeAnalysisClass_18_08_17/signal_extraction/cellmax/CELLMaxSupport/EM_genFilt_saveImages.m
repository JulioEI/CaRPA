function [paramStore, changeParams, options] = EM_genFilt_saveImages(imgs, saveDir, varargin)

  % make directory if it does not exist already
  mkdir(saveDir)

  % input movie properties
  options.oneCentered=1;

  % display
  options.suppressProgressFig=0;

  % outputs
  options.recalculateFinalTraces=1;
  options.doEventDetect=1;

  % movie chunking and restriction
  options.maxSqSize=60;
  options.xLims=[];
  options.yLims=[];
  options.sqOverlap=12;

  % Initialization
  options.initMethod='grid';
  options.gridSpacing=6; % spacing of initialization grid, in pixels
  options.gridWidth=8;

  % Convergence
  options.maxDeltaParams=10^(-5);
  options.minIters=25;
  options.maxIters=300;

  % Removal of overlapping and small images
  options.borderRemoveBuffer=5;
  options.sizeThresh=12;
  options.removeZeroVarImages=1;
  options.removeDiscontig=1;
  options.areaOverlapThresh=0.65;

  % scaledPhi-based removal section
  options.removeCorrProbs=1;      % toggle for turning this on
  options.scaledPhiCorrThresh=0.7;   % correlation threshold
  options.distanceThresh=5;   % distance threshold between centroids (pixels)

  % Core algorithm parameters
  options.numSigmasThresh=0;
  options.numPhotonsPerSigma=5;
  options.numSigmasThreshInitial=options.numSigmasThresh;
  options.removeBelowThreshPixelsForRecalc=1;
  options.useMuChangeOnly=0;
  options.useScaledPhi=1;
  options.threshForElim=0.005;
  options.useConstantBG=1;

  % replace default options with input options
  options.plotDeltaParams=0;
  options=getOptions(options, varargin);

  % set progress square settings
  if isfield(options, 'progressSquare')
      progressSquare=options.progressSquare;
  else
      progressSquare=[];
  end

  % if imgs centered on 1, subtract 1 to center on 0
  % note that this breaks for other movie normalizations
  if options.oneCentered
      imgs=imgs-1;
  end
  imgs(imgs<0)=0;


  % initialize.
  %%%% make initializeEM method!
  [initImages, initProbs] = initializeEM(imgs, 'options', options);



  % iterate through EM cycles
  options.threshForElim=min(options.noiseSigma, options.threshForElim);
  thisChangeParams=10000000;
  numIters=1;
  if ~isempty(progressSquare)
      figure(777); hold off; imagesc(progressSquare); colormap(gray); drawnow;
  end
  paramStore=cell(options.maxIters,1);
  changeParams=zeros(options.maxIters,1);
  cEst=initImages;
  phiEst=initProbs;
  muEst = getImageCentroids(cEst,size(imgs(:,:,1)));%*********************************************
  if options.useScaledPhi
      paramStore{numIters}=[muEst, cEst, zeros(size(phiEst))];%*********************************************
  else
      paramStore{numIters}=[muEst, cEst];%*********************************************
  end

  numSigmasThreshFinal=options.numSigmasThresh;
  numSigmasThreshInitial=options.numSigmasThreshInitial;
  changeSigmaThreshDelta=10^(-5.5);
  pastThreshold=0;
  while (thisChangeParams>options.maxDeltaParams || numIters<options.minIters) && numIters<options.maxIters+1

      if thisChangeParams>changeSigmaThreshDelta && ~pastThreshold
          options.numSigmasThresh=numSigmasThreshInitial;
      else
          options.numSigmasThresh=numSigmasThreshFinal;
          pastThreshold=1;
      end

      % save max image of cells so far
      cellImages=reshape(cEst', [size(imgs(:,:,1)) size(cEst,1)]);
      %cellImages=normalizeImages(cellImages);
      maxImg=max(cellImages(20:end-20,20:end-20,1:end-1),[],3);
      figh=figure; imagesc(maxImg); colormap(gray);
  %     if numIters==1
  %         cLims=get(gca, 'CLim');
  %     else
  %         set(gca, 'CLim', cLims);
  %     end
      cLims(1)=min(maxImg(:));
      cLims(2)=prctile(maxImg(:),99.8);
      set(gca, 'CLim', cLims)
      set(gca, 'XTick', [], 'YTick', [])
      print([saveDir filesep 'CellMap_iter' num2str(numIters) '.eps'], '-depsc')
      close(figh); drawnow

      % do one EM iteration
      [cEstNew, phiEstNew, toDelete, scaledPhi] = EM_genFilt_oneIteration(imgs, cEst, phiEst,'options', options);

      % calculate the change in images/centroids as specified
      muEstNew = getImageCentroids(cEstNew,size(imgs(:,:,1)));%*********************************************
      if options.useMuChangeOnly
          muEst = getImageCentroids(cEst,size(imgs(:,:,1)));%*********************************************
          thisChangeParams=max(abs(muEst(:)-muEstNew(:)));
      else
          thisChangeParams=mean(abs(cEst(:)-cEstNew(:)));
      end
      thisChangeParams(isnan(thisChangeParams))=1000;

      % remove cells with scaledPhi that is too correlated
      if options.removeCorrProbs && mod(numIters-1, 10)==0 && numIters-1>25
          C=corr(scaledPhi');
          for cInd1=1:size(scaledPhi,1)
              for cInd2=cInd1+1:size(scaledPhi,1)
                  if  C(cInd1,cInd2)>options.scaledPhiCorrThresh && ...
                          norm(muEstNew(cInd1,:)-muEstNew(cInd2,:))<options.distanceThresh && ... %%% the centroid is stored in cEst(cellInd,1:2). the rest of that cEst row is the image
                          ~toDelete(cInd2) && ~toDelete(cInd1)
                      cEstNew(cInd1,:)=(cEstNew(cInd1,:)+cEstNew(cInd2,:))/2;
                      scaledPhi(cInd1,:)=(scaledPhi(cInd1,:)+scaledPhi(cInd2,:))/2;
                      toDelete(cInd2)=1;
  %                     %%% To visually inspect which images are getting
  %                       combined at this step...
  %                     filtTraces=calculateFilteredTraces(imgs, cEst([cInd1,cInd2],:)','options',options);
  %                     eventTimes=detectEventsOnPhi(double(filtTraces), double(scaledPhi([cInd1,cInd2],:)));
  %                     [v,~] = manualCellChecker(imgs, scaledPhi([cInd1,cInd2],:), single(reshape(cEst([cInd1,cInd2],:)', [size(imgs(:,:,1)) 2])),...
  %                         eventTimes);
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
         yMax=max(log10(changeParams(1:numIters-1)));
         yRange=yMax-yMin;
         yMinCenter=yMin+yRange/(2*size(progressSquare,1));
         yMaxCenter=yMax-yRange/(2*size(progressSquare,1));
         if yMax==yMin
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
  end
  paramStore=paramStore(1:numIters);
end