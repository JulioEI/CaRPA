function [interpVector] = interpolMissingFrames(vector,missing)
    originalLength = size(vector,2) + length(missing);
    samplePoints = setdiff(1:originalLength,missing);
    
    interpVector = zeros([size(vector,1),originalLength]);
    for row = 1:size(vector,1)
        interpVector(row,:) = interp1(samplePoints,vector(row,:),1:originalLength,'pchip');
    end  
%     figure;plot(samplePoints,vector(1,:),'-o');hold on;plot(interpVector(1,:),'-o');
end

