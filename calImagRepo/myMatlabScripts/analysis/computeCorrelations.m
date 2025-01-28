function [meanCorr,steCorr,pfCorr,otherMat] = computeCorrelations(x,r)
%X should be binned and one-dimensional

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

%Take the pairs in quantiles
%myTriu = @(M,k) M(logical(triu(ones(size(M)),k)));
%corrElem = myTriu(PFxCorr,1);

pfCorr = [-.75,-.5,-.25,0,.25,.5,.75,1];

quantilePairs = cell([length(pfCorr),1]);
for i = 1:size(PFxCorr,1)
    for j = i:size(PFxCorr,2)
        if i ~= j
            if ~isnan(PFxCorr(i,j))
                belongsTo = find(PFxCorr(i,j) <= pfCorr,1,'first');
                quantilePairs{belongsTo} = [quantilePairs{belongsTo};[i,j]];
            end
        end
    end
end
%
%%Compute the activity xcorrelation for everything
%xCorr = corr(r,r);
%Compute the noise coreraltions
xCorr = computeNoiseCorrelations(x,r);

%sigma = det(covCor(xCorr));

%For each pair in each group look up the correlation
xCorrG = cell(size(quantilePairs));
for g = 1:length(quantilePairs)
    for pair = 1:size(quantilePairs{g},1)
        xCorrG{g} = [xCorrG{g};xCorr(quantilePairs{g}(pair,1),quantilePairs{g}(pair,2))];
    end
end

%Compute the mean and ste's
meanCorr = cellfun(@nanmean,xCorrG);
steCorr = cellfun(@nanstd,xCorrG)./sqrt(cellfun(@numel,xCorrG));

if nargout > 3
    otherMat.pairs = quantilePairs;
    otherMat.AcorrPerPair = xCorrG;
    
    PFxCorrPairs = cell(size(quantilePairs));
    for g = 1:length(quantilePairs)
        for pair = 1:size(quantilePairs{g},1)
            PFxCorrPairs{g} = [PFxCorrPairs{g};PFxCorr(quantilePairs{g}(pair,1),quantilePairs{g}(pair,2))];
        end
    end
    otherMat.PFcorrPerPair = PFxCorrPairs;
end
end
