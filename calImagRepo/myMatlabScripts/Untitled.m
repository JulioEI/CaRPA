%%
%Create fake surrogate data
rSurrogate = zeros(size(r));
lambda = 5;
treshold = 1.5;
spikes = poissrnd(lambda,size(r)) < treshold;
%%
gammaDist = makedist('Gamma',2,5);
myX = 0:1:100;
gPdf = gammaDist.pdf(myX);
hold on

%%
for celli = 1:size(spikes,2)
    thisConv = conv(double(spikes(:,celli)),gPdf');
    rSurrogate(:,celli) = thisConv(1:length(x));
end
figure;
plot(conv(double(spikes(:,celli)),gPdf'),'linewidth',2)
axis off