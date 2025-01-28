function [icImgs,icTraces]=pcaIca(imgs,nPCs,nICs,muVal,varargin)

    options.selectPCs=0;
    options.selectICs=0;
    options.displayIcImgs=0;
    options.TermTolICs=10^(-5);
    options=getOptions(options,varargin);


    [pcImgs, pcTraces]=my_pca(imgs, nPCs+5);

    pcImgs=pcImgs(:,:,1:min(nPCs,size(pcImgs,3)));
    pcTraces=pcTraces(1:min(nPCs,size(pcImgs,3)),:);
    nPCs=size(pcImgs,3);

    if options.selectPCs
        pcsToKeep=1:nPCs;
        h=figure;
        for pcInd=1:25
            figure(h)
            imagesc(pcImgs(:,:,pcInd))
            colormap(gray)
            title(['PC ' num2str(pcInd) 'Mouse Click to keep, space bar to throw out'])
            k=waitforbuttonpress;
            if k==1
                pcsToKeep(pcInd)=0;
            end
        end
        pcsToKeep(pcsToKeep==0)=[];
        pcImgs=pcImgs(:,:,pcsToKeep);
        pcTraces=pcTraces(pcsToKeep,:);
        close(h)
    end

    [icImgs,icTraces]=my_ica(imgs, pcImgs, pcTraces, nICs, muVal, 'TermTolICs', options.TermTolICs);
    icImgs=permute(icImgs,[2 3 1]);

    if options.displayIcImgs
        h=figure;
        for icNum=1:nICs
            figure(h)
            subplot(4,1,1:3)
            imagesc(icImgs(:,:,icNum))
            subplot(4,1,4)
            plot(icTraces(icNum,:))
            pause(1)
        end
        close(h)
    end

    nICs=size(icImgs,3);
    if options.selectICs
        icsToKeep=1:nICs;
        h=figure;
        for icInd=1:nICs
            figure(h)
            imagesc(icImgs(:,:,icInd))
            colormap(gray)
            title(['IC ' num2str(icInd) 'Mouse Click to keep, space bar to throw out'])
            k=waitforbuttonpress;
            if k==1
                icsToKeep(icInd)=0;
            end
        end
        icsToKeep(icsToKeep==0)=[];
        icImgs=icImgs(:,:,icsToKeep);
        icTraces=icTraces(icsToKeep,:);
        close(h)
    end
end