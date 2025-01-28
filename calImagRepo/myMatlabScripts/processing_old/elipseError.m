close all

myPath = 'D:\processing\david\Mouse2009\20141120-icx\model_comparison';
confidence = .95;
plotData = true;
myColormap = lines;

methodNameLog = {};
plotLog = {};
allFiles = dir(myPath);
figure()
for i = 3:length(allFiles) %Two first files are windows garbage.
    myColor = myColormap(i,:);
    methodNameFull = allFiles(i).name;
    load([myPath,'\',methodNameFull])
    methodName = methodNameFull(1:end-4);
    methodNameLog = [methodNameLog,methodName];
    p1 = drawElipseError(eval(methodName),myColor,confidence,plotData);
    plotLog = [plotLog,p1];
end

legend(plotLog,methodNameLog)
xlabel('SNR');ylabel('S-ratio');box off;

