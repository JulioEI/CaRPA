function [position,velocity,score] = postProcessTracesTracker(traceLog,smoothingWindowSize,lonlyValLenTresh,doPlot)

if nargin < 4
    doPlot = 0;
    if nargin < 3
        lonlyValLenTresh = 5;
        if nargin < 2
            smoothingWindowSize = 10;
        end
    end
end

frameN = length(traceLog);

velLog = diff(traceLog);
if doPlot;figure;subplot(2,1,1);plot(velLog);legend('x','y');subplot(2,1,2);plot(traceLog);legend('x','y');end

absVel = sqrt(velLog(:,2).^2 +velLog(:,1).^2);
badIdx = (abs(absVel) > nanmean(absVel)+3*nanstd(absVel));

traceLogR = traceLog;
traceLogR(badIdx,:) = nan;

badIdx2 = (abs(traceLogR(:,1)) > nanmean(traceLogR(:,1))+3*nanstd(traceLogR(:,1)) | abs(traceLogR(:,2)) > nanmean(traceLogR(:,2))+3*nanstd(traceLogR(:,2)));
traceLogR(badIdx2,:) = nan;
velLogR = diff(traceLogR);

if doPlot;figure;subplot(2,1,1);plot(traceLogR);legend('x','y');subplot(2,1,2);plot(velLogR);legend('x','y');end

%%
%Remove small isolated trajectories
traceReversed = traceLogR(end:-1:1,:);
nonNanVal = find(~isnan(traceReversed(:,1)) & ~isnan(traceReversed(:,2)));

% figure;plot(traceReversed,'.')
% hold on;
% line([find(isnan(traceReversed(:,1)) | isnan(traceReversed(:,2))) find(isnan(traceReversed(:,1)) | isnan(traceReversed(:,2)))],[0,800],'Color',[.9 .9 .9])

nonNalValDiff = diff(nonNanVal);
lonlyVal = find(1~=nonNalValDiff);
lonlyValLen = diff(lonlyVal);
badLen = find(lonlyValLen<=lonlyValLenTresh);
badIdx3 = [];
for k = 1:length(badLen)
    badIdx3 = [badIdx3,1+nonNanVal(lonlyVal(badLen(k))):nonNanVal(lonlyVal(1+badLen(k)))];
end
badIdx3 = unique(badIdx3);

traceReversed(badIdx3) = nan;

traceLogR = traceReversed(end:-1:1,:);
velLogR = diff(traceLogR);

if doPlot;figure;subplot(2,1,1);plot(traceLogR);legend('x','y');subplot(2,1,2);plot(velLogR);legend('x','y');end
%%
%Interpolate values
score = sum(isnan(traceLogR(:)))/numel(traceLogR);
disp([num2str(100*score),'% of the trace is nan'])
x = find(~isnan(traceLogR(:,1))&~isnan(traceLogR(:,2)));
y = traceLogR(x,:);

traceLogInterP = [interp1(x,y(:,1),1:frameN,'pchip')',interp1(x,y(:,2),1:frameN,'pchip')'];
velLogInterP = diff(traceLogInterP);

if doPlot;figure;subplot(2,1,1);plot(traceLogInterP);legend('x','y');subplot(2,1,2);plot(velLogInterP);legend('x','y');end

%Smooth

% Smooth signal
b = (1/smoothingWindowSize)*ones(1,smoothingWindowSize);
a = 1;
traceLogS = filtfilt(b,a,traceLogInterP);
velLogS = diff(traceLogS);

if doPlot;figure;subplot(2,1,1);plot(traceLogS);legend('x','y');subplot(2,1,2);plot(velLogS);legend('x','y');end

% diffSig = traceLogR - traceLogFilt;
% badIdx2 = abs(diffSig) > nanmedian(diffSig)+3*nanstd(diffSig));
% traceLogR(badIdx2,:) = nan;

traceLogF = traceLogS;
% velLogF = velLogS;

position = traceLogF;
velocity = sqrt(diff(position(:,1)).^2 + diff(position(:,2)).^2);
velocity = [median(velocity(1:50));velocity];
end

