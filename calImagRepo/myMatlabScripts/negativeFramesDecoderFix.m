% %%%%%%%%%%%%%%%%%%%%%
% %%%Load everything%%%
% %%%%%%%%%%%%%%%%%%%%%
load('E:\Processing\Mouse2021\Mouse-2021-20150321-linear-track\Mouse-2021-20150321_103313-linear-track-emAnalysis.mat')%load('E:\Processing\Mouse2029\Mouse-2029-20150321-linear-track\Mouse-2029-20150321_122709-linear-track-emAnalysis.mat')
load('E:\Processing\Mouse2021\Mouse-2021-20150321-linear-track\Mouse-2021-20150321_103313-linear-track-emAnalysisSorted.mat')%load('E:\Processing\Mouse2029\Mouse-2029-20150321-linear-track\Mouse-2029-20150321_122709-linear-track-emAnalysisSorted.mat')
calciumStr = 'E:\Processing\Mouse2021\Mouse-2021-20150321-linear-track\Mouse-2021-20150321_103313-linear-track-dfof.h5';%'E:\Processing\Mouse2029\Mouse-2029-20150321-linear-track\Mouse-2029-20150321_122709-linear-track-dfof.h5';
behavior = 'E:\Processing\Mouse2021\Mouse-2021-20150321-linear-track\Mouse-2021-20150321_103313-linear-track-behavior.avi';%'E:\Processing\Mouse2029\Mouse-2029-20150321-linear-track\Mouse-2029-20150321_122709-linear-track-behavior.avi';
missingFramesStr = 'E:\Processing\Mouse2021\Mouse-2021-20150321-linear-track\Mouse-2021-20150321_103313-linear-track-log.txt';%'E:\Processing\Mouse2029\Mouse-2029-20150321-linear-track\Mouse-2029-20150321_122709-linear-track-log.txt';
% 
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%Compute resampled traces and mouse position%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%position = getMouseTrajectory(behavior);
% calcium = fileSolution.readHdf5(calciumStr);
% disp('Computing upscaled traces...')
% robustMovie = calcium;
% robustMovie(isnan(calcium)) = min(calcium(:));
% [upScaledPhi, upFilteredTraces] = recalcPhiAndDetectEvents(robustMovie,emAnalysisOutput.cellImages(:,:,logical(validCellMax)),emAnalysisOutput.dsCellTraces(logical(validCellMax),:),emAnalysisOutput.CELLMaxoptions,'runEventDetection',0);
% upScaledPhi = upScaledPhi*10^5;
% 
% %%
% %%%%%%%%%%%%%%%%%%%%
% %%%Get the inputs%%%
% %%%%%%%%%%%%%%%%%%%%
traces = upScaledPhi';
x = dataAnalysis.binPosition(position,20);
missingFrames = readFramesFromLog(missingFramesStr);
clearvars -except traces x missingFrames
%%z`
%%%%%%%%%%%%%%%%%%%
%%%Do the search%%%
%%%%%%%%%%%%%%%%%%%

%Separate positive and negative frames
frameBlocks = {missingFrames(missingFrames>0),missingFrames(missingFrames<=0)};

%Reset first frame idx's to 1
frameBlocks = cellfun(@(x) x-min(x)+1,frameBlocks,'UniformOutput',0);

%Assume the first block comes before the second, compute scores
disp(repmat('%',[3,100]));disp('Computing scores of first block before second');
[bestShiftFrames.AB,bestTrace.AB,bestTraceScore.AB,shiftScore.AB] = ABshiftOptima(traces,x(:,1),frameBlocks{1},frameBlocks{2});

%Assume the second block comes before the first, compute scores
disp(repmat('%',[3,100]));disp('Computing scores of second block before first');
[bestShiftFrames.BA,bestTrace.BA,bestTraceScore.BA,shiftScore.BA] = ABshiftOptima(traces,x(:,1),frameBlocks{1},frameBlocks{2});

%Compare results and get the optimal frames
if bestTraceScore.AB < bestTraceScore.BA
    optimalFrames = [bestShiftFrames.AB.A, bestShiftFrames.AB.B];
else
    optimalFrames = [bestShiftFrames.BA.A, bestShiftFrames.BA.B];
end