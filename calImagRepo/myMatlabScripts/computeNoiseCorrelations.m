function [ncMean,ncMat] = computeNoiseCorrelations(x,r)

if numel(x) ~= prod(size(x))
    error('x should be one dimensional')
end

ncMat = [];
for bin = min(x):max(x)
    rBin = r(x==bin,:);
    ncMat = cat(3,ncMat,corr(rBin,rBin));    
end

%Average over bins

ncMean = nanmean(ncMat,3);

% %Take all the pairs
% 
% myTriu = @(M,k) M(logical(triu(ones(size(M)),k)));
% nc = myTriu(ncMean,1);

end

