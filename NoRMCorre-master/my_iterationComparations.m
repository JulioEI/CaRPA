%% read data and convert to double
name = 'E:\Processing\newTurboreg\Mouse2028\20150305-icx\concat_recording_20150305_132231.h5';
%addpath(genpath('../../NoRMCorre'));
Yf = read_file(name);
Yf = single(Yf);
Yf = Yf(:,:,1:1000);
[d1,d2,T] = size(Yf);
%% perform some sort of deblurring/high pass filtering

if (0)    
    hLarge = fspecial('average', 40);
    hSmall = fspecial('average', 2); 
    for t = 1:T
        Y(:,:,t) = filter2(hSmall,Yf(:,:,t)) - filter2(hLarge, Yf(:,:,t));
    end
    %Ypc = Yf - Y;
    bound = size(hLarge,1);
else
    gSig = 7; 
    gSiz = 17; 
    psf = fspecial('gaussian', round(gSiz), gSig);
    ind_nonzero = (psf(:)>=max(psf(:,1)));
    psf = psf-mean(psf(ind_nonzero));
    psf(~ind_nonzero) = 0;   % only use pixels within the center disk
    Y = imfilter(Yf,psf,'same');
    %Ypc = Yf - Y;
    bound = 2*ceil(gSiz/2);
end

itVect = [1,2,4,6,8];
cM2Log = [];
mM2Log = [];
vM2Log = [];
cM2fLog = [];
mM2fLog = [];
vM2fLog = [];
for k = 1:length(itVect)
    %% now apply non-rigid motion correction
    % non-rigid motion correction is likely to produce very similar results
    % since there is no raster scanning effect in wide field imaging

    % %% register data and apply shifts to removed percentile
    % %% first try out rigid motion correction
    %     % exclude boundaries due to high pass filtering effects
    options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',50,'max_shift',20,'iter',itVect(k),'correct_bidir',false);

    tic; [M1,shifts1,template1] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r); toc % register filtered data
        % exclude boundaries due to high pass filtering effects
    % tic; Mr = apply_shifts(Yf,shifts1,options_r,bound/2,bound/2); toc % apply shifts to full dataset
        % apply shifts on the whole movie


    options_nr = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',50, ...
        'grid_size',[128,128],'mot_uf',4,'correct_bidir',false, ...
        'overlap_pre',32,'overlap_post',32,'max_shift',20,'iter',itVect(k));

    tic; [M2,shifts2,template2] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_nr,template1); toc % register filtered data
    movie = double(Yf);
    myDFOF = movie./mean(movie,3) - 1;
    movie = single(myDFOF);
    tic; Mpr = apply_shifts(movie,shifts2,options_nr,bound/2,bound/2); toc % apply the shifts to the removed percentile
    

    %% compute metrics



    [cM2,mM2,vM2] = motion_metrics(M2,options_nr.max_shift);
    cM2Log = [cM2Log,cM2];
    mM2Log = cat(3,mM2Log,mM2);
    vM2Log = [vM2Log,vM2];
    [cM2f,mM2f,vM2f] = motion_metrics(Mpr,options_nr.max_shift);
    cM2fLog = [cM2fLog,cM2f];
    mM2fLog = cat(3,mM2fLog,mM2f);
    vM2fLog = [vM2fLog,vM2f];
end

[cY,mY,vY] = motion_metrics(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r.max_shift);
[cYf,mYf,vYf] = motion_metrics(Yf,options_r.max_shift);

%%
figure;
iosr.statistics.boxPlot([cY,cM2Log],'OutlierSize',0);
% shadedErrorBar([0,itVect],[mean(cY),mean(cM2Log)],[std(cY),std(cM2Log)],'-o');xticks([0,itVect]);axis square;xlabel('iterations');ylabel('mean correlation')
% shadedErrorBar([0,itVect],[mean(cYf),mean(cM2fLog)],[std(cYf),std(cM2fLog)],'-o');xticks([0,itVect]);axis square;xlabel('iterations');ylabel('mean correlation')
title('mean correlation coefficient of each frame with the mean','fontsize',14,'fontweight','bold')

figure;
plot([0,itVect],[vY,vM2Log],'-o');xticks([0,itVect]);axis square;xlabel('iterations');ylabel('gradient');
% plot([0,itVect],[vYf,vM2fLog],'-o');xticks([0,itVect]);axis square;xlabel('iterations');ylabel('gradient');
title('norm of gradient of mean image','fontsize',14,'fontweight','bold')

%%
caxisY = [quantile(Y(:),0.0005),quantile(Y(:),1-0.0005)];
caxisYf = [quantile(Yf(:),0.0005),quantile(Yf(:),1-0.0005)];
caxisM2 = [quantile(M2(:),0.0005),quantile(M2(:),1-0.0005)];
caxisMpr = [quantile(Mpr(:),0.0005),quantile(Mpr(:),1-0.0005)];
%%
figure;
for k = 1:size(Y,3)
    subplot(2,2,1)
    imagesc(Y(:,:,k));caxis(caxisY);
    title('original highpass')
    subplot(2,2,2)
    imagesc(Yf(:,:,k));caxis(caxisYf);
    title('original')
    subplot(2,2,3)
    imagesc(M2(:,:,k));caxis(caxisM2);
    title('registered highpass')
    subplot(2,2,4)
    imagesc(Mpr(:,:,k));caxis(caxisMpr);
    title('registered')
    drawnow;%pause(.1)
end