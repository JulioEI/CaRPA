cellStruct = getCellStruct({'C:\Users\csc\Desktop\caImagg\Data\Mouse2029'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2026'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2028'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2010'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2012'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2019'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2023'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2025'
    'E:\Processing\Mouse2021'
    'D:\Storage\Processing\Mouse-2027'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2022'
    'C:\Users\csc\Desktop\caImagg\Data\Mouse2024'
    'E:\Processing\Mouse2020'
    'E:\Processing\Mouse2011'
    },cellStruct);
%%

mouseNames = {cellStruct.mouse};
cellsPerSess = cell([1,length(cellStruct)]);
namePerSess = cell([1,length(cellStruct)]);
for mouseK = 1:length(cellStruct)
        cellsPerSess{mouseK} = cat(2,cellStruct(mouseK).days.cells);
        for c = 1:length(cellStruct(mouseK).days)
            if length(cellStruct(mouseK).days(c).cells) > 1
                for k = 1:length(cellStruct(mouseK).days(c).cells)
                    namePerSess{mouseK} = [namePerSess{mouseK};[cellStruct(mouseK).days(c).name,'_',num2str(k)]];
                end
            else
                namePerSess{mouseK} = [namePerSess{mouseK};[cellStruct(mouseK).days(c).name,'  ']];
            end
        end
end
maxSessLen = max(cellfun(@length,cellsPerSess));
[~,sessOrdering] = sort(cellfun(@length,cellsPerSess),'descend');
namePerSess = namePerSess(sessOrdering);
cellsPerSess = cellsPerSess(sessOrdering);
myTable = nan([length(mouseNames),maxSessLen]);

for mouseK = 1:length(cellStruct)
        myTable(mouseK,1:length(cellsPerSess{mouseK})) = cellsPerSess{mouseK};
end
mouseNames = mouseNames(sessOrdering);

allRows = cell([1,size(myTable,2)]);
sessionNum = cell([1,size(myTable,2)]);
for row = 1:size(myTable,2)
    allRows{row} = myTable(:,row);
    sessionNum{row} = ['S', num2str(row)];
end

T = table(allRows{:},'rowNames',mouseNames,'VariableNames',sessionNum);
%writetable(T,'C:\Users\csc\Desktop\caImagg\Data\cellTable.csv','Delimiter','\t','WriteRowNames',true)
%%
%Find animal with max mean cells
[meanCells,meanSort] = sort(nanmean(myTable,2),'descend');
table(meanCells,maxCells','rowNames',mouseNames(meanSort),'VariableNames',{'MeanCells','MaxCells'})
%%
%Find animal with max cells
[maxCells,maxSort] = sort(max(myTable'),'descend');
sessionMaxCells = cell([1,length(maxCells)]);
for k = 1:length(maxCells)
    sessionMaxCellsK = find(myTable(maxSort(k),:) == maxCells(k));
    sessionMaxCells{k} = namePerSess{maxSort(k)}(sessionMaxCellsK,:);
end
%table(maxCells',meanCells,'rowNames',mouseNames(maxSort),'VariableNames',{'MaxCells','MeanCells'})
table(maxCells',meanCells,sessionMaxCells','rowNames',mouseNames(maxSort),'VariableNames',{'MaxCells','MeanCells','Session'})
