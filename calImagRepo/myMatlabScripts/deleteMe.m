close all;
home;
%%
%Generate data

neurons = 2;
time = 100;
trials = 100;
covMat = [20,10;10,30];

data = zeros([neurons,time,trials]);
for k = 1:trials
    for t = 1:time
        meanMat = [t;2*t];
        data(:,t,k) = mvnrnd(meanMat,covMat);
    end
end

%%
%Mean across trials

dataAvgTr = mean(data,3);
figure;
subplot(1,2,1)
plot(dataAvgTr','x');
legend({'Neuron1','Neuron2'});xlabel('Time');ylabel('Amplitude')
subplot(1,2,2)
%mdl = fitlm(dataAvgTr(1,:),dataAvgTr(2,:));
%plot(dataAvgTr(1,:),dataAvgTr(2,:),'-x');xlabel('Neuron1');ylabel('Neuron2');axis square;title(num2str(mdl.Rsquared.Ordinary))
plot2Dpca(dataAvgTr')

%%
%Mean across time

dataAvgTi = squeeze(mean(data,2));
figure;
subplot(1,2,1)
plot(dataAvgTi','x');
legend({'Neuron1','Neuron2'});xlabel('Trial');ylabel('Amplitude')
subplot(1,2,2)
%mdl = fitlm(dataAvgTi(1,:),dataAvgTi(2,:));
%plot(dataAvgTi(1,:),dataAvgTi(2,:),'x');xlabel('Neuron1');ylabel('Neuron2');axis square;title(num2str(mdl.Rsquared.Ordinary))
plot2Dpca(dataAvgTi')


%%


