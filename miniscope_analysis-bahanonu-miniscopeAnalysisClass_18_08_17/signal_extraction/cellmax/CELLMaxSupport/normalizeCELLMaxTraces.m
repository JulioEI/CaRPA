function cellTraces = normalizeCELLMaxTraces(cellTraces,cellImages)

nCells = size(cellTraces,1);

for cInd = 1:nCells
    thisImage = cellImages(:,:,cInd);
    normFactor = max(thisImage(:));
    cellTraces(cInd,:) = cellTraces(cInd,:)*normFactor;
end