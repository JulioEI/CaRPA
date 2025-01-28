function [MprCrop] = register1p(Yf)
    iter = 2;
    Yf = single(Yf);
    [d1,d2,~] = size(Yf);
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

    %% now apply non-rigid motion correction
    % non-rigid motion correction is likely to produce very similar results
    % since there is no raster scanning effect in wide field imaging

    % %% register data and apply shifts to removed percentile
    % %% first try out rigid motion correction
    %     % exclude boundaries due to high pass filtering effects
    options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',50,'max_shift',20,'iter',iter,'correct_bidir',false,'shifts_method','linear');

    tic; [test1,~,template1] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r); toc % register filtered data
        % exclude boundaries due to high pass filtering effects
    % tic; Mr = apply_shifts(Yf,shifts1,options_r,bound/2,bound/2); toc % apply shifts to full dataset
        % apply shifts on the whole movie


    options_nr = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',50, ...
        'grid_size',[128,128],'mot_uf',4,'correct_bidir',false, ...
        'overlap_pre',32,'overlap_post',32,'max_shift',20,'iter',iter,'use_parallel',1,'boundary','NaN','shifts_method','linear');

%     options_nr = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'use_parallel',1,'boundary','NaN');
    
    tic; [mTemp,shifts2,~] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_nr,template1); toc % register filtered data
    tic; MprFull = apply_shifts(Yf(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),shifts2,options_nr,bound/2,bound/2); toc % apply the shifts to the removed percentile
    MprCrop = MprFull(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:);

%     c1Axis = [min(M_final(:)),max(M_final(:))];
%     c2Axis = [min(Mpr(:)),max(Mpr(:))];
%     c3Axis = [min(Yf(:)),max(Yf(:))];
%     figure;
%     for k = 100:size(Mpr,3)
%         subplot(1,3,1);imagesc(M_final(:,:,k));
%         caxis(c1Axis)
%         subplot(1,3,2);imagesc(Mpr(:,:,k));
%         caxis(c2Axis)
%         subplot(1,3,3);imagesc(Yf(:,:,k));
%         caxis(c3Axis)
%         drawnow;
%     end
%     
%     
%     thisClip = mTemp;%MprFull./mean(MprFull,3)-1;
%     thisClip2 = MprFull./mean(MprFull,3)-1;
%     figure;
%     myCAxis = [min(thisClip(:)),max(thisClip(:))];
%     myCAxis2 = [min(thisClip2(:)),max(thisClip2(:))];
%     for k = 1:size(thisClip,3)
%         subplot(1,2,1)
%         imagesc(thisClip(:,:,k))
%         caxis(myCAxis)
%         subplot(1,2,2)
%         imagesc(thisClip2(:,:,k))
%         caxis(myCAxis2)
%         pause(.2)
%     end

end
