function [bestShiftFrames,bestTrace,bestTraceScore,shiftScore] = ABshiftOptima(traces,x,A,B)

%Assuming frames A come before frames B, find optimal of B
traceWithA = interpolMissingFrames(traces',A)';
minShift = max(A)+1;
maxShift = size(traces,1)-1;
shiftScore.B = zeros([1,(maxShift-minShift)+1]);
textprogressbar('Computing first block shift ');
for shift = 0:(maxShift-minShift)
    textprogressbar(100*shift/(maxShift-minShift));
    shiftFrames = minShift + B + shift;
    shiftFrames(shiftFrames > maxShift) = mod(shiftFrames(shiftFrames > maxShift),maxShift) + minShift;
    traceWithAB = interpolMissingFrames(traceWithA',shiftFrames)';
    shiftScore.B(shift+1) = mean(quickBayesDecoder(traceWithAB,x,10));
end
textprogressbar(' done');

%Find optimal of a A from begining to the optimal B
[bestShiftScore.B,bestShift.B] = min(shiftScore.B); if sum(bestShiftScore.B == shiftScore.B) > 1; error('Did not find unique optima for B');end
bestShiftFrames.B = minShift + B + bestShift.B;
bestShiftFrames.B(bestShiftFrames.B > maxShift) = mod(bestShiftFrames.B(bestShiftFrames.B > maxShift),maxShift) + minShift;

traceWithB = interpolMissingFrames(traces',bestShiftFrames.B)';
minShift = 0;
maxShift = min(shiftFrames);
shiftScore.A = zeros([1,(maxShift-minShift)+1]);
textprogressbar('Computing second block shift ');
for shift = 0:(maxShift-minShift)
    textprogressbar(100*shift/(maxShift-minShift));
    shiftFrames = minShift + A + shift;
    shiftFrames(shiftFrames > maxShift) = mod(shiftFrames(shiftFrames > maxShift),maxShift) + minShift;
    traceWithBA = interpolMissingFrames(traceWithB',shiftFrames)';
    shiftScore.A(shift+1) = mean(quickBayesDecoder(traceWithBA,x,10));
end
textprogressbar(' done');

[bestShiftScore.A,bestShift.A] = min(shiftScore.A); if sum(bestShiftScore.A == shiftScore.A) > 1; error('Did not find unique optima for A');end
bestShiftFrames.A = minShift + A + bestShift.A;bestShiftFrames.A(bestShiftFrames.A > maxShift) = mod(bestShiftFrames.A(bestShiftFrames.A > maxShift),maxShift) + minShift;

bestTrace = interpolMissingFrames(traces',[bestShiftFrames.A,bestShiftFrames.B])';
bestTraceScore = mean(quickBayesDecoder(bestTrace,x,10));
end

