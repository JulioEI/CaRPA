function [signalCorrelation,trialByTrialCorrelation] = computeCorrelationPairs(x,r)
    
    %Compute the PF of the cells
    PF = cell([max(x),size(r,2)]);
    for t = 1:length(x)
        for celli = 1:size(r,2)
            PF{x(t),celli} = [PF{x(t),celli},r(t,celli)];    
        end
    end 
    PF = cellfun(@mean,PF);

    %Compute the PF crosscorrelation
    PFxCorr = corr(PF,PF);
    
    %%Compute the activity xcorrelation for everything
    %xCorr = corr(r,r);
    %Compute the noise correlations for everything
    xCorr = computeNoiseCorrelations(x,r);

    signalCorrelation = [];
    trialByTrialCorrelation = [];
    for i = 1:size(PF,1)
        for j = i:size(PF,2)
            if i ~= j
                signalCorrelation = [signalCorrelation,PFxCorr(i,j)];
                trialByTrialCorrelation = [trialByTrialCorrelation,xCorr(i,j)];
            end
        end
    end

end