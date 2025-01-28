function assessAlignmentAcrossDays(cellMaps,globalIDs)
%INPUTS: 
% cellMaps: either the root folder of several days to search cellmaps
% (cellmax) or a cell with all the cellmaps, with separate cells
% globalIDs: matrix with the global idx for each session

%Load cell maps of root folder
if isstr(cellMaps)
    folders = dir(cellMaps);
    allCellMaps = {};
    for folder = {folders.name} %one file per folder!
        allFiles = dir(fullfile(path,folder{1}));
        for file = {allFiles(find(~cellfun(@isempty,(regexp({allFiles.name},'emAnalysis.mat'))))).name};
            disp(file{1});
            decisionRegexp = [file{1}(1:end-4),'Sorted.mat'];
            decisions = load(fullfile(path,folder{1},allFiles(find(~cellfun(@isempty,(regexp({allFiles.name},decisionRegexp))))).name));
            emFile = load(fullfile(path,folder{1},file{1}));
            allCellMaps{k} = emFile.emAnalysisOutput.cellImages(:,:,logical(decisions.validCellMax));
            k = k + 1;
        end
    end
    allCellMaps = cellfun(@(x) thresholdImages(x),allCellMaps,'UniformOutput',0);
end

%Evolution of cells

idxFirstDay = find(globalIDs(:,1));

cellLog = zeros([1,size(globalIDs,2)]);
for k = 1:size(globalIDs,2)
    cellLog(k) = length(find(prod(globalIDs(:,1:k),2)));
end

cellLog2 = zeros([1,size(globalIDs,2)]);
for k = 1:size(globalIDs,2)
    cellLog2(k) = sum(sum(globalIDs(:,setdiff(1:size(globalIDs,2),k)),2) & globalIDs(:,k));
end

cellsInDay = zeros([1,size(globalIDs,2)]);
for k = 1:size(globalIDs,2)
    cellsInDay(k) = length(find(globalIDs(:,k)));
end

cellLog3 = zeros([1,size(globalIDs,2)]);
for k = 1:size(globalIDs,2)
    cellsInOtherSessions = sum(globalIDs(:,setdiff(1:size(globalIDs,2),k))~=0,2);
    cellsInOtherSessions = cellsInOtherSessions(cellsInOtherSessions&globalIDs(:,k));
    cellLog3(k) = sum(cellsInOtherSessions > floor(size(globalIDs,2)/2));
end

cellMat = zeros(size(globalIDs,2));
for j = 1:size(globalIDs,2)
    for k = 1:size(globalIDs,2)
        cellsInOtherSessions = sum(globalIDs(:,setdiff(1:size(globalIDs,2),k))~=0,2);
        cellsInOtherSessions = cellsInOtherSessions(cellsInOtherSessions&globalIDs(:,k));
        cellMat(j,k) = 100*sum(cellsInOtherSessions > floor(size(globalIDs,2)-j))/cellsInDay(k);
    end
end

figure;imagesc(flipud(cellMat));axis square;xlabel('days');ylabel('shared in at least');colorbar;set(gca, 'YTick', 1:size(globalIDs,2));set(gca, 'XTick', 1:size(globalIDs,2));set(gca,'YDir','normal')


figure;bar([100*cellLog./cellsInDay;100*cellLog2./cellsInDay;100*cellLog3./cellsInDay]')
ylabel('%Cells of day')
xlabel('Day')
legend('%Cells aligned to all previous days','%Cells aligned to any other day','%Cells in half the other sessions')


end

