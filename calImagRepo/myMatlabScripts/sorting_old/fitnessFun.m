function [ y ] = fitnessFun(x,myData,allValidCellMax)

    fieldNames = fields(myData);
    pred = ones([length(myData.(fieldNames{1})),1]);
    centers = x(1:2:length(x));
    ranges = x(2:2:length(x));
    limits = zeros([1,length(x)]);
    limits(1:2:length(x)) = centers - ranges;
    limits(2:2:length(x)) = centers + ranges;
    for k = 1:2:length(limits)
        pred = pred & myData.(fieldNames{ceil(k/2)}) > limits(k) & myData.(fieldNames{ceil(k/2)}) < limits(k+1);
    end
    
%     FP = sum(pred'==1 & allValidCellMax==0)/sum(allValidCellMax==0);
%     TP = sum(pred'==1 & allValidCellMax==1)/sum(allValidCellMax==1);
%     FN = sum(pred'==0 & allValidCellMax==1)/sum(allValidCellMax==1);
%     TN = sum(pred'==0 & allValidCellMax==0)/sum(allValidCellMax==0);

    y = confusion(allValidCellMax,pred');
%    y = (FP+FN)/(FP+FN+TP+TN);

end

