function cellImgs  = calcCellImgs(params, imgSize)

% Written by Lacey Kitch in 2013

% params = [MuX, MuY, SigX, SigY, Theta]

nCells=size(params,1);
[xdata, ydata]=meshgrid(1:imgSize(2), 1:imgSize(1));

if size(params,2)==5

    mux=params(:,1);
    muy=params(:,2);
    sigx=params(:,3);
    sigy=params(:,4);
    theta=params(:,5);

    a=(cos(theta).^2)./(2*sigx.^2)+(sin(theta).^2)./(2.*sigy.^2);
    b=-sin(2.*theta)./(4.*sigx.^2)+sin(2.*theta)./(4.*sigy.^2);
    c=(sin(theta).^2)./(2.*sigx.^2)+(cos(theta).^2)./(2.*sigy.^2);

    cellImgs=zeros([imgSize, nCells]);

    for cInd=1:nCells
        cellImgs(:,:,cInd)=exp(-(a(cInd).*(xdata-mux(cInd)).^2+...
            2.*b(cInd).*(xdata-mux(cInd)).*(ydata-muy(cInd))+...
            c(cInd).*(ydata-muy(cInd)).^2));
    end

    % from file calcCellImgs on 12/11/13 at 6:04pm

else

    muEst=params(:,1:2);
    sigEst=params(:,3:6);

    cellImgs=zeros([imgSize, nCells]);

    xdata=xdata(:);
    ydata=ydata(:);

    for cInd=1:nCells
        if sum(sigEst(cInd,:))>0
            P = mvnpdf([xdata, ydata],muEst(cInd,:),reshape(sigEst(cInd,:),[2 2]));
            %
            if sum(P(:))>0
                cellImgs(:,:,cInd)=reshape(P,imgSize)/max(P(:));
            else
                cellImgs(:,:,cInd)=0;
            end
        end
    end
end