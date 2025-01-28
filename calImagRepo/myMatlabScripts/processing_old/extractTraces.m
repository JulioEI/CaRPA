cImag = emAnalysisOutput.cellImages;

cTraces = zeros([size(cImag,3),size(inputMovie,3)]);
for j = 1:size(cImag,3)
    convMovie = zeros(size(inputMovie));
    for k = 1:size(convMovie,3)
        convMovie(:,:,k) = cImag(:,:,j).*inputMovie(:,:,k);
    end
    cTraces(j,:) = squeeze(nanmean(nanmean(convMovie)));
end
