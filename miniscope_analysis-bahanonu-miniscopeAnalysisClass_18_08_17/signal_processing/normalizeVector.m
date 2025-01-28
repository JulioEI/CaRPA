function [outputVector] = normalizeVector(inputVector,varargin)
    % normalizes a vector or matrix (2D or 3D) using the method specified
    % biafra ahanonu
    % started: 2014.01.14 [23:42:58]
    % inputs
    	%
    % outputs
    	%

    % changelog
    	% 2014.02.17 added zero centered normalization
    % TODO
    	%

    %========================
    options.normRange = 'oneToNegativeOne';
    % get options
    options = getOptions(options,varargin);
    % display(options)
    % unpack options into current workspace
    % fn=fieldnames(options);
    % for i=1:length(fn)
    % 	eval([fn{i} '=options.' fn{i} ';']);
    % end
    %========================
    maxVec = nanmax(inputVector(:));
    minVec = nanmin(inputVector(:));
    meanVec = nanmean(inputVector(:));
    switch options.normRange
        case 'oneToNegativeOne'
            outputVector = ((inputVector-minVec)./(maxVec-minVec) - 0.5 ) *2;
        case 'oneToOne'
            outputVector = (inputVector-minVec)./(maxVec-minVec);
        case 'zeroToOne'
            outputVector = (inputVector-minVec)./(maxVec-minVec);
        case 'zeroCentered'
            vectorMean = nanmean(inputVector,2);
            outputVector = bsxfun(@rdivide,inputVector,vectorMean)-1;
        case 'zeroCenteredCorrect'
            outputVector = inputVector-meanVec;
        case 'dfof'
            outputVector = inputVector/meanVec-1;
        otherwise
            % body
    end
end