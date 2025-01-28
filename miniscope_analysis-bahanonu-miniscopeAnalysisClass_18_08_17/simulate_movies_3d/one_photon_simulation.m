%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Use two-photon volumetric movies as the ground truth to simulate the one-photon movie 
%   by modeling the tissue scattering, optical system and camera sampling.
%
%   Xing Lin [xinglin@stanford.edu]
%   March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Setup System Parameters%%%%%%%%%%%%%%%%%%%%%%%%
NA = 0.5;  % numerical aperture of the objective lens
n = 1.33;  % the refractive index of the sample
lambda = 525*1e-9;  % wavelength
opticalResolution = 2.5*1e-6;  % optical resolution
samplingResolution = 2.5*1e-6;  % pixel sampling resolution
camera_pixelPitch = 6.5*1e-6;  % camera pixel size
M = camera_pixelPitch/samplingResolution;  % total system magnification
pixelSize=100;  % the number of the pixels for each PSF

zmin = -14*10*1e-6;  % minimum axial position of 2p volume
zmax = 13*10*1e-6;  % maximum axial position of 2p volume
zspacing = 10*1e-6;  % axial increment of 2p volume
pixelPitch = 460/512*M*1e-6;  % pixel size of 2p volume
dz = zmin:zspacing:zmax;

enable_scattering = 1;
enable_noise = 1;

fprintf('Simulating One-photon Image...\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Generating System PSF%%%%%%%%%%%%%%%%%%%%%%%%%%
[psfStack_ideal, psfStack_aberration] = ...
    generating_system_psf(M,NA,n,lambda,zmin,zmax,zspacing,pixelPitch,pixelSize,opticalResolution);

psfStack_aberration = psfStack_aberration./sum(sum(sum(psfStack_aberration)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Volumetric Covolution%%%%%%%%%%%%%%%%%%%%%%
if enable_scattering==1
psfScattering_amount=7; % scattering amount, arbitrary unit
h_1 = fspecial_3d('gaussian',psfScattering_amount);
end

frame=1;

for i=2:length(dz)+1
    twoPhoton_stack(:,:,i-1)=double(imread(['one_frame/TSeries-10312016-1609-',num2str(sprintf('%03d',i))...
        ,'_Cycle001_CurrentSettings_Ch2_',num2str(sprintf('%06d',frame)),'.tif']));
%     twoPhoton_stack(:,:,i-1)=imresize(twoPhoton_stack(:,:,i-1),1,'bicubic');
end

if enable_scattering==1
twoPhoton_stack=imfilter(twoPhoton_stack,h_1);
end

for i=2:length(dz)+1
    conv_stack(:,:,i-1)=imfilter(twoPhoton_stack(:,:,i-1),psfStack_aberration(:,:,i-1),'replicate');
end
onePhoton_org=sum(conv_stack,3);
onePhoton=single(imresize(onePhoton_org,...
    [round( size(onePhoton_org,1)*(pixelPitch./M)./(camera_pixelPitch./M) ),...
     round( size(onePhoton_org,2)*(pixelPitch./M)./(camera_pixelPitch./M) )],'bicubic'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Adding Poisson Noise%%%%%%%%%%%%%%%%%%%%%%%
if enable_noise==1
noise_amount=35; % PSNR dB
onePhoton = adding_poisson_noise(onePhoton,noise_amount);
end

%%%%%%%%%%%%%%%%%%%%%%Save Simulated One-photon Data%%%%%%%%%%%%%%%%%%%%%%%
t = Tiff(['./one_photon_',num2str(sprintf('%06d',frame)),'.tif'], 'w');
tagstruct.ImageLength = size(onePhoton,1);
tagstruct.ImageWidth = size(onePhoton,2);
tagstruct.Photometric = 1;
tagstruct.BitsPerSample = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB';
tagstruct.SampleFormat = 3;
t.setTag(tagstruct);
t.write(onePhoton);  
t.close();

fprintf('Simulated Image Saved!\n\n')

