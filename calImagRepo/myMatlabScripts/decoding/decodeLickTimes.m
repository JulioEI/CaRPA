% moviePath = 'C:\Users\csc\Desktop\caImagg\Behavior\Linear-track\Mouse-2028\Behavioral-videos\Mouse2028_20150303_103504.avi';
% mousePath = 'E:\Processing\Mouse2028 - 500itEM\20150303-icx\old\';
% [position,~] = getMouseTrajectory(moviePath);
% load([mousePath,'upFilteredTraces.mat'])
% 
load([mousePath,'upScaledPhi.mat'])
load([mousePath,'validCells.mat'])
scaledPhi = upScaledPhi(logical(validCells),:);
% validCellsIdx = find(validCells);
% spikeTimes = cell([1,length(validCellsIdx)]);
% for k = 1:length(validCellsIdx)
%     [c, s, options] = deconvolveCa(upFilteredTraces(validCellsIdx(k),:)','exp2','foopsi');
%     spikeTimes{k} = find(s)';
% end
% figure;plot(upFilteredTraces(validCellsIdx(1),(1:1000))');hold on;plot(c(1:1000));stem(s(1:1000),'-','marker','none');legend('scaledPhi','denoised','spikeTrain')
mov = VideoReader(moviePath);

sess = 1;

period = 1/20;

if size(position,1) ~= mov.numberOfFrames
    error('Movie and calcium have movie differ in frames')
end
realWidth = 1;
metersPerPixel = realWidth/mov.Width;
Trajectory{sess} = [period*(0:(size(position,1)-1))',position*metersPerPixel];
Trajectory{sess}(:,2) = Trajectory{sess}(:,2) - (max(Trajectory{sess}(:,2)) + min(Trajectory{sess}(:,2)))/2; 
Trajectory{sess}(:,3) = Trajectory{sess}(:,3) - (max(Trajectory{sess}(:,3)) + min(Trajectory{sess}(:,3)))/2; 


load('C:\Users\csc\Desktop\caImagg\Behavior\Linear-track\Mouse-2028\Time-events\Mouse-2028-03-Mar-2015licking_events.mat')
eventsSec = time_events(:,4)*60*60+time_events(:,5)*60+time_events(:,6);
eventsSec = eventsSec-eventsSec(1);

eventsFra = zeros([1,length(eventsSec)]);
for k = 1:length(eventsSec)
    tdiff = (eventsSec(k) - Trajectory{sess}(:,1)).^2;
    [~,ind] = min(tdiff);    
    eventsFra(k) = ind(1);
end
absTraX = abs(Trajectory{sess}(:,2));
figure;plot(absTraX)
absTraXDiff = abs(diff(absTraX));
absTraX(absTraXDiff > median(absTraXDiff)) = nan;
% hold on;plot(absTraX)
traXNaN = find(~isnan(absTraX));
traXNaNDiff = diff(traXNaN);
traXFillIdx = find(traXNaNDiff > 1 & traXNaNDiff < 30); 
traXFill = [traXNaN(traXFillIdx),traXNaN(traXFillIdx)+traXNaNDiff(traXFillIdx)];
for k = 1:size(traXFill,1)
    absTraX(traXFill(k,1):traXFill(k,2)) = abs(Trajectory{sess}(traXFill(k,1):traXFill(k,2),2));
end
absTraX(absTraX<max(absTraX)*.9) = nan;
hold on;plot(absTraX)

idxTillOutflow = length(absTraX) - max(eventsFra);
scoreVect = zeros([1,idxTillOutflow]);
for k = 1:idxTillOutflow
    scoreVect(k) = sum(~isnan(absTraX(eventsFra+k)));
end
bestScore = find(scoreVect==max(scoreVect));

if length(bestScore) > 1
    stdVect = zeros([1,length(bestScore)]);
    for k = 1:length(bestScore)
        stdVect(k) = nanstd(absTraX(eventsFra+bestScore(k)));
    end
    bestScore = bestScore(find(stdVect==min(stdVect),1,'first'));
end

eventsFra = eventsFra + bestScore;

figure;
plot(abs(Trajectory{sess}(:,2)))
hold on
plot([eventsFra;eventsFra],get(gca,'ylim'),'color',[.8,.8,.8])
hold off

% figure;
% for k = 1:mov.numberOfFrames
%     subplot(2,1,1);imagesc(mov.read(k));
%     subplot(2,1,2);
%     plot(Trajectory{sess}(1:k,2))
%     hold on
%     try;plot([eventsFra(eventsFra<=k);eventsFra(eventsFra<=k)],get(gca,'ylim'),'color',[.8,.8,.8]);end
%     hold off
%     drawnow;
% end

timeOfEvents = eventsFra(end:-1:1)*period;
afterAllLicks = Trajectory{sess}(:,1)>timeOfEvents(1);
Trajectory{sess}(afterAllLicks,:) = [];
scaledPhi(:,afterAllLicks) = [];

timeTillGoal = zeros([size(Trajectory{sess},1),1]);
for k = 1:length(timeOfEvents)
    idx = (Trajectory{sess}(:,1)<=timeOfEvents(k));
    timeTillGoal(idx) = Trajectory{sess}(idx,1) - timeOfEvents(k);
end


normalizedSP = scaledPhi./repmat(max(scaledPhi')',[1,size(scaledPhi,2)]);
x = normalizedSP';
y = timeTillGoal;
