%DPRIME
%Load data
load('D:\Storage\Processing\Mouse2024\Mouse-2024-20150317-linear-track\Mouse-2024-20150317_202708-linear-track-TracesAndEvents.mat')
%%
%Set up parameters
param.dtCamera = 0.05; %period of miniscope (seconds/frame)
param.dT = 0.05; %Integration parameter (seconds)
param.binN = [20,1]; %Number of bins in X Y
param.minVel = [4, 0]; %Frames under this speed will be discarded (cm/s)
param.pxPerCm = [5.62, 5.62]; %How many px a centimeter is (px/cm)

x = tracesEvents.position;%svmT.fw_ks;%
r = tracesEvents.rawProb;%svmT.fw_X;%
%Compute dprime
disp('-----------------------')
disp('Computing dPrime of raw')
shuffleTrials = 0;
[dVectRawMy,dPrimeDistRawMy,otherRawMy] = computeDPrime(x,r,param,shuffleTrials);
disp('-----------------------')
disp('Computing dPrime of shuffle')
shuffleTrials = 1;
[dVectShuffleMy,dPrimeDistShuffleMy,otherShuffleMy] = computeDPrime(x,r,param,shuffleTrials);
%%
%Plot results
figure;
h = shadedErrorBar(dVectRawMy(2:end),dPrimeDistRawMy.mean(2:end),1.96*dPrimeDistRawMy.ste(2:end),'lineprops',{'-o','color',[178, 83, 0]/255,'linewidth',1},'patchSaturation',0.2);h.mainLine.DisplayName = 'Raw';
h = shadedErrorBar(dVectShuffleMy(2:end),dPrimeDistShuffleMy.mean(2:end),1.96*dPrimeDistShuffleMy.ste(2:end),'lineprops',{'-o','color',[70, 160, 178]/255,'linewidth',1},'patchSaturation',0.2);h.mainLine.DisplayName = 'Shuffle';
xlabel('Distance (bins)');ylabel('d''');axis square;grid on;
legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'),'Location','best');
figText = ['m2024 SVM ','rawProb ',num2str(1000*0.05),'msec'];
annotation('textbox',[0.3,0,1,1],'String',figText,'FitBoxToText','on','FontSize',12);

