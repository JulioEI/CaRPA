function [output] = infoFromData(extractionFile, decisions, calcium, data, missingFrames, behavior)

%Get traces and events
tracesAndEvents = getTracesAndEvents(extractionFile, decisions, calcium, missingFrames,behavior);
for fn = fieldnames(tracesAndEvents)'
   output.(fn{1}) = tracesAndEvents.(fn{1});
end

%Return onset of events
disp('Computing event onsets')
NumberCalciumFrames = size(output.rawTraces,1);
eventOnsetVector = zeros([1,NumberCalciumFrames]);
cueTime=data.datatimes.cueTime;

if(isfield(data.datatimes,'FrameTimes')==1)
    FrameTimes=data.datatimes.FrameTimes;  
    for k = 1:(length(cueTime)-1)        
        Indx=[]; Indx = round((cueTime{k}(1)-FrameTimes(1,2))*20);
        if(0<Indx)
            destinationFrame = Indx;
            if Indx > NumberCalciumFrames
                warning('Shorter calcium movie')
            end
            eventOnsetVector(destinationFrame) = cueTime{k}(2);
        end
    end
else   
    %Early version of the experiment file format
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

