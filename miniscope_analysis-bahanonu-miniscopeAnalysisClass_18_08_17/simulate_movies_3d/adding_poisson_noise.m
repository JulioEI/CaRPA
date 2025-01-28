function imgdata_noise = adding_poisson_noise(imgdata,noise_amount)
%%%%%%%%%%%%Adding Poisson Noise to Obtain the Certain PSNR%%%%%%%%%%%%%%%%

count=1;
delta=0.1;
photon_convert_min=1;
photon_convert_max=5;

for photon_convert=photon_convert_min:delta:photon_convert_max  %arbitrary unit
imgdata_noise_org(:,:,count)=uint16(imgdata*photon_convert);
imgdata_noise(:,:,count)=single(imnoise(imgdata_noise_org(:,:,count),'poisson'))/photon_convert;
[peaksnr(count), snr(count)] = psnr(mat2gray(imgdata_noise(:,:,count)), mat2gray(imgdata));
count=count+1;
end

peaksnr=abs(peaksnr-noise_amount);
index=find(peaksnr==min(peaksnr));
imgdata_noise=imgdata_noise(:,:,index);

[peaksnr, snr] = psnr(mat2gray(imgdata_noise), mat2gray(imgdata));
end
