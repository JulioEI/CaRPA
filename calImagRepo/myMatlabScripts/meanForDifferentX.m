function [x,y,ste,N] = meanForDifferentX(allX,allY,myBinEdges)
    
    if nargin < 3
        myBinEdges = 'auto';
    end
    
    if ~iscell(allX)
        allX = {allX};
    end
    if ~iscell(allY)
        allY = {allY};
    end
    
    joinX = cat(find(size(allX)==max(size(allX)),1),allX{:});
    joinY = cat(find(size(allY)==max(size(allY)),1),allY{:});
    
    if strcmp(myBinEdges,'auto')
        [N,edges] = histcounts(joinX);
    else
        edges = myBinEdges;
    end
    
    binnedY = cell([1,length(edges)-1]);
    for k = 1:length(binnedY)
        binnedY{k} = joinY(joinX >= edges(k) & joinX < edges(k+1));
    end
    
    N = cellfun(@length,(binnedY));
    
    x = edges(1:end-1) + diff(edges)/2;
    y = cellfun(@mean,binnedY);
    ste = cellfun(@(x) std(x)./sqrt(length(allY)),binnedY);
    
end

