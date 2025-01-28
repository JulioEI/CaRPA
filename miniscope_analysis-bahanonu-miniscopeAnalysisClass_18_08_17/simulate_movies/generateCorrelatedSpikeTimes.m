function [spikeTimes,C] = generateCorrelatedSpikeTimes(spikeRate, nCells, meanCorr, totTime)

meanRates=spikeRate*ones(nCells,1);
C=zeros(nCells,nCells);
spikeVar=spikeRate;
for i=1:nCells
    C(i,i)=spikeVar;
    if meanCorr>0
        for j=i+1:nCells
            thisCorr=2*meanCorr*rand;
            thisCorr(thisCorr>0.99)=0.99;
            thisCov=thisCorr*spikeVar;
            C(i,j)=thisCov;
            C(j,i)=thisCov;
        end
    end
end

spikeYesNo = sampleCovPoisson(meanRates,C,totTime);
spikeTimes=cell(nCells,1);
for i=1:nCells
    spikeTimes{i}=find(spikeYesNo(i,:)>0);
end