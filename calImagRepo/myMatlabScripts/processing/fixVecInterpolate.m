function xFix = fixVecInterpolate(x,myLen,interpMethod)
    if nargin < 3
        interpMethod = 'pchip';
    end
    
    samplePoints = linspace(1,myLen,length(x));

    xFix = zeros([myLen,size(x,2)]);
    for k = 1:size(xFix,2)
        xFix(:,k) = interp1(samplePoints,x(:,k),1:myLen,interpMethod);
    end  
    
end

