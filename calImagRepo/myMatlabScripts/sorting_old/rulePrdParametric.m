function [pred,newData] = rulePrdParametric(myData,lBound,uBound)
    fieldNames = fields(myData);
    pred = ones([length(myData.(fieldNames{1})),1]);
    for k = 1:length(fieldNames)
        pred = pred & myData.(fieldNames{k}) > lBound(k) & myData.(fieldNames{k}) < uBound(k);
    end
    
    [~,~,wheresNan] = removeNansAnd3FromData(myData);
    maskVect = ones([length(myData.(fieldNames{1})),1]);
    maskVect(wheresNan) = 0; %Mask to predict as bad the cells with nan in some of their features (this may not be ideal).
    pred = maskVect & pred;
    
    for k = 1:length(fieldNames)
        newData.(fieldNames{k}) = myData.(fieldNames{k})(pred);
    end
end

