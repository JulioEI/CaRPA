function [cellImgs,finalCellParams,finalCellTraces,goodCellInds] = resolveBorderConflicts(cellFitParams,...
    areaOverlapThresh,imgs,calcTraces,cellImgs,varargin)

    % Written by Lacey Kitch in 2013
    % 2017.01.18 - updated to include support for movie chunks - biafra

    if ~isempty(varargin)
    	options=varargin{1};
    else
        if isa(cellImgs, 'struct')
            options=cellImgs;
            cellImgs=[];
        else
            options=[];
        end
    end

    % old way of saving, only temporary until full switch
    optionsTwo.datasetName = '/1';
    optionsTwo.displayInfo = 1;
    optionsTwo.movieFilename = '';
    optionsTwo.movieDatasetName = '/Data/Images';
    optionsTwo.readMovieChunks = 0;
    % get options
    optionsTwo = getOptions(optionsTwo,varargin);

    if ~ (calcTraces==1 || calcTraces==0)
        calcTraces %#ok<NOPRT>
        error('ResolveBorderConflicts being called incorrectly');
    end
    if optionsTwo.readMovieChunks==0
        imgSize=size(imgs(:,:,1));
    else
        imgSize = [imgs(1) imgs(2)];
    end

    nCells=size(cellFitParams,1);
    if size(cellFitParams,2)>2
        cellImgs=calcCellImgs(cellFitParams,imgSize);
    end
    cellImgsCopy=cellImgs;
    for cInd=1:nCells
        thisImg=cellImgs(:,:,cInd);
        thisImg(thisImg<0.2*max(thisImg(:)))=0;
        cellImgs(:,:,cInd)=thisImg;
    end
    cellImgs(cellImgs>0)=1;

    finalCellParams=cellFitParams;

    cellsToDelete=nan(nCells,1);
    nCellsBad=0;

    reverseStr = '';
    for cInd=1:nCells
        if mod(cInd,25)==0
            reverseStr = cmdWaitbar(cInd,nCells,reverseStr,'inputStr','resolving conflicts','waitbarOn',1,'displayEvery',25);
        end

        c1Pixels=find(cellImgs(:,:,cInd)>0);
        if ~ismember(cInd,cellsToDelete)
            for matchInd=(cInd+1):nCells
                if ~ismember(cInd,cellsToDelete) && ~ismember(matchInd,cellsToDelete)
                    if abs(cellFitParams(cInd,1)-cellFitParams(matchInd,1))<20 && abs(cellFitParams(cInd,2)-cellFitParams(matchInd,2))<20
                        c2Pixels=find(cellImgs(:,:,matchInd)>0);
                        thisOverlap=length(intersect(c1Pixels,c2Pixels));
                        thisOverlap1=thisOverlap/length(c1Pixels);
                        thisOverlap2=thisOverlap/length(c2Pixels);

                        if thisOverlap1>areaOverlapThresh && thisOverlap2>areaOverlapThresh
                            nCellsBad=nCellsBad+1;
                            if length(c1Pixels)<length(c2Pixels)
                                cellsToDelete(nCellsBad)=cInd;
                            else
                                cellsToDelete(nCellsBad)=matchInd;
                            end
                        end
                    end
                end
            end
        end
    end
    cellsToDelete=cellsToDelete(1:nCellsBad);
    finalCellParams(cellsToDelete,:)=[];
    cellImgs=cellImgsCopy;
    cellImgs(:,:,cellsToDelete)=[];
    goodCellInds=1:nCells;
    goodCellInds(cellsToDelete)=[];

    if optionsTwo.readMovieChunks==0
        if calcTraces
            finalCellTraces = calculateTraces(cellImgs, imgs, options);
        else
            finalCellTraces=[];
        end
    else
        if calcTraces
            finalCellTraces = calculateTraces(cellImgs, imgs, 'options',optionsTwo);
        else
            finalCellTraces=[];
        end
    end

end