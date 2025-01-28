function [myData,allValidCellMax,wheresNan] = removeNansAnd3FromData(myData,allValidCellMax)
    fieldNames = fields(myData);
    nanIdx = [];
    for k = 1:length(fieldNames)
        nanIdx = [nanIdx;find(isnan(myData.(fieldNames{k})))];
    end
    wheresNan = unique(nanIdx);

    for k = 1:length(fieldNames)
        myData.(fieldNames{k})(wheresNan) = [];
    end
    if nargin > 1
        allValidCellMax(wheresNan) = [];
        allValidCellMax(allValidCellMax==3)=0;
    else
        allValidCellMax = 'WARNING: No truth vector specified';
    end
end

