function [psfStack_ideal, psfStack_aberration] = ...
    generating_system_psf(M,NA,n,lambda,zmin,zmax,zspacing,pixelPitch,pixelSize,opticalResolution)
%%%%%%%%Generating System PSF Stack with Fresnel Propagation%%%%%%%%%%%%%%%

fov = pixelSize*pixelPitch;   %the size of field of view for the PSF
k = 2*pi/lambda;
dz = zmin:zspacing:zmax;

%%%%%%create the coordinates in Fourier plane using sinalpha as unit%%%%%%%
sinalpha_max = NA / n / M;
fx_sinalpha = lambda / 2 / pixelPitch / lambda;
fx_step = lambda / fov / lambda;
fx_max = fx_sinalpha ;
fx= -fx_max : fx_step : fx_max-fx_step;
[fxcoor fycoor] = meshgrid( fx , fx );
fx2coor=fxcoor.*fxcoor;
fy2coor=fycoor.*fycoor;

%%%%%%%%%create the mask for aperture size%%%%%%%%%%%%
aperture_mask=((fx2coor+fy2coor)<=((sinalpha_max/lambda*n).^2));
%%%%%%%%%%%%first step detection PSF%%%%%%%%%%%%%%%%%%
psfWAVE_f=ones(pixelSize,pixelSize).*aperture_mask;
psfWAVE=fftshift(ifft2(ifftshift(psfWAVE_f)));

%%%%%%%%%%%%%%%%%%calculate defocus aberration matrix%%%%%%%%%%%%%%%%%%%%%%
psfWAVE_fAFTER_STACK=zeros(pixelSize,pixelSize,length(dz));
psfWAVE_AFTER_STACK=zeros(pixelSize,pixelSize,length(dz));
psfStack_ideal=zeros(pixelSize,pixelSize,length(dz));

for i=1:length(dz)
    tempZ=dz(i);
    tempP=k*n*tempZ*realsqrt((1-(fxcoor.*lambda./n.*M).^2-(fycoor.*lambda./n.*M).^2).*aperture_mask);
    psfWAVE_fAFTER_STACK(:,:,i)=psfWAVE_f.*exp(1j*tempP);
    psfWAVE_AFTER_STACK(:,:,i)=fftshift(ifft2(ifftshift(squeeze(psfWAVE_fAFTER_STACK(:,:,i)))));
    psfStack_ideal(:,:,i)=abs(squeeze(psfWAVE_AFTER_STACK(:,:,i)));
    psfStack_ideal(:,:,i)=psfStack_ideal(:,:,i)./sum(sum(squeeze(psfStack_ideal(:,:,i))));   
end

%%%%%%%%%%%%%%%%%%%%introducing optical aberrations%%%%%%%%%%%%%%%%%%%%%%%%
pixelPitch_objectSide=pixelPitch./M;
psfAberration_pixelNumber=opticalResolution./pixelPitch_objectSide;

count=1;
step=0.1;
for psfAberration_sigma=step:step:psfAberration_pixelNumber
    psfAberration = fspecial('gaussian', pixelSize,psfAberration_sigma);
    psfAberration_focus=imfilter(psfStack_ideal(:,:,floor(-zmin/zspacing)+1),psfAberration,'replicate');
    psfAberration_focus_mask = psfAberration_focus >= max(max(psfAberration_focus))./2;  %FWHM Resolution
    psfAberration_focus_size(count)=length(find(psfAberration_focus_mask(pixelSize./2,:)~=0));
    count=count+1;
end

psfAberration_focus_size=abs(psfAberration_focus_size-psfAberration_pixelNumber);
index=find(psfAberration_focus_size==min(psfAberration_focus_size));
psfAberration = fspecial('gaussian', pixelSize,mean(index)*step);  % Aberration PSF
psfAberration_focus=imfilter(psfStack_ideal(:,:,floor(-zmin/zspacing)+1),psfAberration,'replicate');
psfAberration_focus_mask = psfAberration_focus >= max(max(psfAberration_focus))./2;
psfAberration_focus_size=length(find(psfAberration_focus_mask(pixelSize./2,:)~=0));

%%%%%%%%%%%%%convolving ideal PSF with optical aberrations%%%%%%%%%%%%%%%%%
for i=1:length(dz)
    psfStack_aberration(:,:,i)=imfilter(psfStack_ideal(:,:,i),psfAberration,'replicate');
%     psfStack_Chromatic_Aberration(:,:,i)=conv2(psfStack_Chromatic(:,:,i),psfAberration,'same');
    psfStack_aberration(:,:,i)=psfStack_aberration(:,:,i)./sum(sum(psfStack_aberration(:,:,i)));
end
end
