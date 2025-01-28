function valid = sortCells( inputSignals,testpeaksArray,cellImages,inputMovie )

%valid = sortCells( rawSignals,signalPeaksArray,rawImages,ioptions.inputMovie );
myData = getCellInfo(inputSignals,testpeaksArray,permute(cellImages,[2,3,1]),inputMovie);
load C:\Users\csc\Desktop\caImagg\calImagRepo\mdl.mat

%Compact the data into a observationsXfeatures matrix
dataFields = fieldnames(myData);
dataFields = dataFields(1:end);

compactData = zeros([size(myData.(dataFields{1}),1),length(dataFields)]);
for k = 1:length(dataFields)
    compactData(:,k) = myData.(dataFields{k}); %No normalitzation
end
    

valid = zeros([1,size(compactData,1)]);
mask = find(sum(isnan(compactData),2)==0);
validCrop = mdl.predict(compactData(mask,:))';

if iscell(validCrop)
    validCrop = cellfun(@str2num,validCrop);
end

valid(mask) = validCrop;
   
end

