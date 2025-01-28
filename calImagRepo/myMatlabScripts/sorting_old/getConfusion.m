function [errorVal] = getConfusion(targets,outputs)

trueP = sum(targets == 1 & outputs == 1);
falseP = sum(targets == 0 & outputs == 1);
trueN = sum(targets == 0 & outputs == 0);
falseN = sum(targets == 1 & outputs == 0);

totalP = trueP+falseP;
totalN = trueN+falseN;
errorVal = (falseN+falseP)/(totalP+totalN);


disp(['trueP: ',num2str(trueP),' ',num2str(100*trueP/totalP),'%'])
disp(['falseP: ',num2str(falseP),' ',num2str(100*falseP/totalP),'%'])
disp(['trueN: ',num2str(trueN),' ',num2str(100*trueN/totalN),'%'])
disp(['falseN: ',num2str(falseN),' ',num2str(100*falseN/totalN),'%'])
disp(['Error: ',num2str(errorVal),'%'])
end

