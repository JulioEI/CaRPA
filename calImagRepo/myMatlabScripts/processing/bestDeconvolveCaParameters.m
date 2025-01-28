myTraces = emAnalysisOutput.cellTraces(find(validCellMax),:)';
% normalizedUT = (tracesToGetEvents'./repmat(max(tracesToGetEvents)',[1,size(tracesToGetEvents,1)]));

kernelTypes = {'exp2','thresholded';'ar2','thresholded';'ar1','foopsi'};

%%
%Look at a single trace
celli = 10;
k = 1;
figure;
for celli = 10:100
clf;
[myC,myS] = deconvolveCa(myTraces(:,celli),kernelTypes{k,1},kernelTypes{k,2});
mySFind = find(myS>0.10);mySProb = myS(myS>0.10);
subplot(1,2,1);histogram(mySProb(mySProb<.5),100);title(['Spike histogram of ',kernelTypes{k,1},' for cell ',num2str(celli)]);
subplot(1,2,2);plot(0.05*(1:length(myC)),myC,'k','lineWidth',2);hold on;plot(0.05*(1:length(myTraces)),myTraces(:,celli));grid on;title([num2str(length(mySFind)),' spikes'])
plot(0.05*[find(myS>0),find(myS>0)]',([zeros([1,length(myS(myS>0))])',myS(myS>0)])'+repmat([0,0],[length(myS(myS>0)),1])','k','lineWidth',1.7)
plot(0.05*[mySFind,mySFind]',([zeros([1,length(mySProb)])',mySProb])'+repmat([0,0],[length(mySFind),1])','r','lineWidth',1.7)
pause;
end
% myS20 = mySFind(find(mySProb > 0));
% myS50 = mySFind(find(mySProb > .15));
% myS80 = mySFind(find(mySProb > .5));
% myS100 = mySFind(find(mySProb > 1));
% 
% plot(0.05*[myS20,myS20]',repmat([-2,-1.5],[length(myS20),1])','r')
% plot(0.05*[myS50,myS50]',repmat([-1.5,-1],[length(myS50),1])','r')
% plot(0.05*[myS80,myS80]',repmat([-1,-.5],[length(myS80),1])','r')
% plot(0.05*[myS100,myS100]',repmat([-.5,0],[length(myS100),1])','r')

%%
%Compare spike trains
celli = 10;
figure;plot(0.05*(1:length(myTraces)),myTraces(:,celli));hold on;grid on
myColormap = lines;
k = 1;
[myC,myS] = deconvolveCa(myTraces(:,celli),kernelTypes{k,1},kernelTypes{k,2});
mySFind = find(myS);mySProb = myS(myS~=0);
plot((k-1)*.005+0.05*[mySFind,mySFind]',([zeros([1,length(mySProb)])',mySProb])'+repmat([0,0],[length(mySFind),1])','lineWidth',1.7,'Color',myColormap(k+1,:))
k = 2;
[myC,myS] = deconvolveCa(myTraces(:,celli),kernelTypes{k,1},kernelTypes{k,2});
mySFind = find(myS);mySProb = myS(myS~=0);
plot((k-1)*.005+0.05*[mySFind,mySFind]',([zeros([1,length(mySProb)])',mySProb])'+repmat([0,0],[length(mySFind),1])','lineWidth',1.7,'Color',myColormap(k+1,:))
k = 3;
[myC,myS] = deconvolveCa(myTraces(:,celli),kernelTypes{k,1},kernelTypes{k,2});
mySFind = find(myS);mySProb = myS(myS~=0);
plot((k-1)*.005+0.05*[mySFind,mySFind]',([zeros([1,length(mySProb)])',mySProb])'+repmat([0,0],[length(mySFind),1])','lineWidth',1.7,'Color',myColormap(k+1,:))

h = [plot(nan,nan,'Color',myColormap(1,:)),plot(nan,nan,'Color',myColormap(2,:)),plot(nan,nan,'Color',myColormap(3,:)),plot(nan,nan,'Color',myColormap(4,:))];
legend(h,'Original','exp2','ar2','ar1')

%%

%Compute spike trains (s) and deconvolved traces(c)
c = nan([size(myTraces),size(kernelTypes,1)]);
s = nan([size(myTraces),size(kernelTypes,1)]);
for k = 1:size(kernelTypes,1)
    disp(k/size(kernelTypes,1))
    for celli = 1:size(myTraces,2)
        try
            [c(:,celli,k), s(:,celli,k)] = deconvolveCa(myTraces(:,celli),kernelTypes{k,1},kernelTypes{k,2});
        catch
            disp(['Unable to perform ',kernelTypes{k,1},' for cell ',num2str(celli)])
        end
    end
end
%%
%Compute mse between trace and deconvolved trace
cVec = zeros([size(myTraces,2),size(kernelTypes,1)]);
for k = 1:size(kernelTypes,1)
    for celli = 1:size(myTraces,2)
        cVec(celli,k) = mean((myTraces(:,celli)-c(:,celli,k)).^2);
    end
end
%%
%Compute mean and std
cMean = nanmean(cVec,1); 
cStd = nanstd(cVec,1);

%Display results
figure;
barwitherr(cStd,cMean)
set(gca,'XTickLabel',kernelTypes(:,1))

% figure;
% plot(normalizedUT(k,:));hold on;plot(c)
% hold on;
% for i = 1:size(sM,2)
%     s = sM(:,i);
%     stem(find(s~=0)+(i-1)*.1,s(s~=0));
% end