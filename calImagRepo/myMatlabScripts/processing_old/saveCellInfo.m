saveTo = [obj.dataPath{1},'\cellInfo'];
mkdir(saveTo)
inputSignals = rawSignals;
testpeaksArray = signalPeaksArray;
cellImages = rawImages;
inputMovie = ioptions.inputMovie;
validCellMax = valid;

for variable = {'inputSignals','testpeaksArray','cellImages','validCellMax'}
	save([saveTo,'\',variable{1},'.mat'],variable{1})
end

save([saveTo,'\','inputMovie','.mat-v7.3'],'inputMovie','-v7.3')