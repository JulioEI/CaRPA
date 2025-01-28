orig = 'E:\Processing\Mouse2028 - 200itEM\20150303-icx\concat_recording_20150303_103622.h5';
orig = read_file(orig);
orig = single(orig);
orig = orig(:,:,1:1000);
gSig = 7; 
gSiz = 17; 
psf = fspecial('gaussian', round(gSiz), gSig);
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;   % only use pixels within the center disk
Yorig = imfilter(orig,psf,'same');
bound = 2*ceil(gSiz/2);

reg1 = 'E:\Processing\Mouse2028 - 200itEM\20150303-icx\2015_03_03_p000_mouse2028_NULL000_turboreg_1.h5';
reg1 = read_file(reg1);
reg1 = single(reg1);
reg1 = reg1(:,:,1:1000);
Yreg1 = imfilter(reg1,psf,'same');

reg2 = 'E:\Processing\newTurboreg\Mouse2028 - reg - linear\20150303-icx\2015_03_03_p000_mouse2028_NULL000_turboreg_1.h5';
reg2 = read_file(reg2);
reg2 = single(reg2);
reg2 = reg2(:,:,1:1000);
Yreg2 = imfilter(reg2,psf,'same');

options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',50,'max_shift',20,'iter',2,'correct_bidir',false);
[cYorig,mYorig,vYorig] = motion_metrics(Yorig(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r.max_shift);
[cYreg1,mYreg1,vYreg1] = motion_metrics(Yreg1(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r.max_shift);
[cYreg2,mYreg2,vYreg2] = motion_metrics(Yreg2(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r.max_shift);

figure;plot(cYorig);hold on;plot(cYreg1);plot(cYreg2);legend('Original','Turboreg','NoRMCorre');title('Correlation between frames and mean');xlabel('Frames');ylabel('p');
set(gcf,'color','white')

figure;
% boxplot([cYorig,cYreg1,cYreg2])
iosr.statistics.boxPlot([cYorig,cYreg1,cYreg2],'OutlierSize',0);
% h = barwitherr([std(cYorig),std(cYreg1),std(cYreg2)], [mean(cYorig),mean(cYreg1),mean(cYreg2)]);
set(gca,'XTickLabel',{'Original','Turboreg','NoRMCorre'})
% shadedErrorBar([0,itVect],[mean(cYorig),mean(cM2Log)],[std(cY),std(cM2Log)],'-o');xticks([0,itVect]);axis square;xlabel('iterations');ylabel('mean correlation')
% shadedErrorBar([0,itVect],[mean(cYf),mean(cM2fLog)],[std(cYf),std(cM2fLog)],'-o');xticks([0,itVect]);axis square;xlabel('iterations');ylabel('mean correlation')
title('mean correlation coefficient of each frame with the mean','fontsize',14,'fontweight','bold')

figure;
h = bar([vYorig,vYreg1,vYreg2]);
set(gca,'XTickLabel',{'Original','Turboreg','NoRMCorre'})
set(h(1),'FaceColor',[52/255,152/255,219/255]);
% plot([0,itVect],[vY,vM2Log],'-o');xticks([0,itVect]);axis square;xlabel('iterations');ylabel('gradient');
% plot([0,itVect],[vYf,vM2fLog],'-o');xticks([0,itVect]);axis square;xlabel('iterations');ylabel('gradient');
title('norm of gradient of mean image','fontsize',14,'fontweight','bold')

%% read data and convert to double
%addpath(genpath('../../NoRMCorre'));
Yf = read_file(name);
Yf = single(Yf);
Yf = Yf(:,:,1:1000);
[d1,d2,T] = size(Yf);
%% perform some sort of deblurring/high pass filtering

gSig = 7; 
gSiz = 17; 
psf = fspecial('gaussian', round(gSiz), gSig);
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;   % only use pixels within the center disk
Y = imfilter(Yf,psf,'same');
%Ypc = Yf - Y;
bound = 2*ceil(gSiz/2);

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
shadedErrorBar([0,itVect],[mean(cY),mean(cM2Log)],[std(cY),std(cM2Log)],'-o');xticks([0,itVect]);axis square;xlabel('iterations');ylabel('mean correlation')
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