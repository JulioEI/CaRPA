function [output] = infoFromData(extractionFile, decisions, calcium, data, missingFrames, behavior)

%Returns struct with traces, events and trajectories
if nargin == 6
    [position,velocity] = getMouseTrajectory(behavior);
    output.position = position;
    output.velocity = velocity;     
end

%Get traces and events
tracesAndEvents = getTracesAndEvents(extractionFile, decisions, calcium, missingFrames);
for fn = fieldnames(tracesAndEvents)'
   output.(fn{1}) = tracesAndEvents.(fn{1});
end

%Return onset of events
disp('Computing event onsets')
NumberCalciumFrames = size(calcium,3);
eventOnsetVector = zeros([1,NumberCalciumFrames]);
cueTime=data.datatimes.cueTime;

if(isfield(data.datatimes,'FrameTimes')==1)
    FrameTimes=data.datatimes.FrameTimes;
    NumberRegistFrames = size(FrameTimes,1);
    
    Indx2=[];Indx2=find(diff(FrameTimes(:,2))>1.71*mode(diff(FrameTimes(:,2))));
    for kk=1:length(Indx2)
        FrameTimes([Indx2(kk)+1:size(FrameTimes)+1],:)=FrameTimes([Indx2(kk):size(FrameTimes)],:);
        FrameTimes(Indx2(kk),2)=mean([FrameTimes(Indx2(kk)-1,2),FrameTimes(Indx2(kk)+1,2)]);
        FrameTimes(Indx2(kk),1)=0;
    end    
    for k = 1:(length(cueTime)-1)
        Indx=[]; Indx = find(0<(FrameTimes(:,2)-cueTime{k}(1)),1,'first');
        destinationFrame = Indx;
        if Indx > NumberCalciumFrames
            warning('Shorter calcium movie')
        end
        eventOnsetVector(destinationFrame) = cueTime{k}(2);
    end
else    
    for k = 1:(length(cueTime)-1)
        EndExperiment=cueTime{end}(1);
        Indx=[]; Indx = round((EndExperiment-cueTime{k}(1))*20);
        destinationFrame = NumberCalciumFrames - Indx + 1 ;
        if Indx > NumberCalciumFrames
            warning('Shorter calcium movie')
        end
        eventOnsetVector(destinationFrame) = cueTime{k}(2);
    end
    
end
output.eventOnsetVector = eventOnsetVector;

end

