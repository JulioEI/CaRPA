function [commonCenters] = findcommonCenters(globalIDs,images,i,j)
    commonIn1And2 = find(prod(globalIDs(:,[i,j]),2));
    commonCenters = zeros([length(commonIn1And2),2,2]);
    for k = 1:length(commonIn1And2)
        i1 = images{i}(:,:,globalIDs(commonIn1And2(k),i));
        i2 = images{j}(:,:,globalIDs(commonIn1And2(k),j));
        [cx1,cy1] = find(max((max(i1)))==i1);
        [cx2,cy2] = find(max((max(i2)))==i2);
        commonCenters(k,1,:) = [cx1,cy1];
        commonCenters(k,2,:) = [cx2,cy2];
    end
end

