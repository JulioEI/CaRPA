%Load data
load('C:\Users\csc\Desktop\caImagg\Data\Mouse2024\Mouse-2024-20150311-linear-track\Mouse-2024-20150311_073912-linear-track-TracesAndEvents.mat')
movie = 'C:\Users\csc\Desktop\caImagg\Data\Mouse2024\Mouse-2024-20150311-linear-track\Mouse-2024-20150311_073912-linear-track-behavior.avi';
%%
%Get position and responses
param.dtCamera = 0.05; %period of miniscope (seconds/frame)
param.dT = 0.05; %Integration parameter (seconds)
param.binN = [20,1]; %Number of bins in X Y
param.minVel = [4, 0]; %Frames under this speed will be discarded (cm/s)
param.pxPerCm = [6.5, 6.5]; %How many px a centimeter is (px/cm)

x = tracesEvents.position;
r = zscore(tracesEvents.rawProb);
[x,r] = dataAnalysis.curateXandR(x,r,param); %Bin, remove slow frames
%r = decoderAnalysis.shuffleKeepingTrials1D(r,x(:,1)); %Shuffle sample times in the same bin
%%
%Do pca on the full time
coeff = pca(r);
meanR = mean(r);
T = (r-meanR)*coeff + meanR; %T is pca space
T12 = T(:,1:2);%Grab the first 2 components
T13 = T(:,1:3);%Grab the first 2 components
%%
%Normalize position
xNorm = (x(:,1)-min(x(:,1)))/(max(x(:,1))-min(x(:,1)));
%%
% %Show pc activity dependent on the position (joining the dots)
% a = T12(:,1)';
% b = T12(:,2)';
% c = zeros(size(a));
% col = xNorm';
% colSat = col;
% %colSat(col > quantile(col,.95)) = quantile(col,.95);
% %colSat(col < quantile(col,.05)) = quantile(col,.05);
% figure;surface([a;a],[b;b],[c;c],[colSat;colSat],...
% 'facecol','no',...
% 'edgecol','interp',...
% 'linew',1);
% colorbar
% axis square;
%%
%Show pc activity dependent on the position (without joining the dots)
figure;
cMap = parula(length(xNorm)+1);
plot(T12(:,1),T12(:,2),'--','color',[0.7,0.7,0.7]);
hold on;
for t = 1:size(T12,1)
    plot(T12(t,1),T12(t,2),'.','color',cMap(1+floor(xNorm(t)*length(xNorm)),:),'markersize',10);
end
colorbar
%%
%Show 3pc (without joining the dots)
figure;
cMap = parula(length(xNorm)+1);
plot3(T13(:,1),T13(:,2),T13(:,3),'--','color',[0.7,0.7,0.7]);
hold on;
for t = 1:size(T13,1)
    plot3(T13(t,1),T13(t,2),T13(t,3),'.','color',cMap(1+floor(xNorm(t)*length(xNorm)),:),'markersize',10);
end
%plot3(T13(1,1),T13(1,2),T13(1,3),'r.','markersize',30);%start
colorbar