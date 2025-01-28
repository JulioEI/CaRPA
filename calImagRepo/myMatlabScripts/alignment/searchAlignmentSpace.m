%%%%%%%%%%%%%%%%%%%%%%%
%%%LOAD THE CELLMAPS%%%
%%%%%%%%%%%%%%%%%%%%%%%

path = 'E:\Processing\Mouse2028';
%folders = {'Mouse-2028-20150301-linear-track','Mouse-2028-20150303-linear-track','Mouse-2028-20150305-linear-track','Mouse-2028-20150307-linear-track','Mouse-2028-20150309-linear-track','Mouse-2028-20150311-linear-track','Mouse-2028-20150313-linear-track'};
folders = {'Mouse-2028-20150301-linear-track','Mouse-2028-20150303-linear-track','Mouse-2028-20150305-linear-track','Mouse-2028-20150307-linear-track','Mouse-2028-20150309-linear-track','Mouse-2028-20150311-linear-track','Mouse-2028-20150313-linear-track','Mouse-2028-20150315-linear-track','Mouse-2028-20150317-linear-track','Mouse-2028-20150319-linear-track','Mouse-2028-20150321-linear-track','Mouse-2028-20150323-linear-track','Mouse-2028-20150325-linear-track','Mouse-2028-20150327-linear-track'};
k = 1;
allCellMaps = {};
for folder = folders %one file per folder!
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
allCellMapsMax = cellfun(@(x) max(x,[],3),allCellMaps,'UniformOutput',0);
%%
%%%%%%%%%%%%%%%%%
%%%OPTIMIZE TB%%%
%%%%%%%%%%%%%%%%%
fitFunTB = @(x) TBcostFn(x,allCellMaps,allCellMapsMax);

smoothing = [1,10,40,80,200];
iterations = [1,3,7,10];
perf = zeros([length(smoothing),length(iterations)]);
for a = 1:length(smoothing)
    disp(repmat('%',[3,100]))
    disp([num2str(a/length(smoothing)),'%'])
    for b = 1:length(iterations)
        perf(a,b) = TBcostFn([smoothing(a),iterations(b)],allCellMaps,allCellMapsMax);
    end
end
   %% 
%%%%%%%%%%%%%%%%%
%%%OPTIMIZE NC%%%
%%%%%%%%%%%%%%%%%
gridSize = [.2,.5,.7];
minPatchSize = [1,5,10];
maxShift = [20,200,500];
iterations = [1,2,5,10];
us_fac = [10,20,50];
mot_uf = [1,4,8];

fitFunNC = @(x) NCcostFn(x,d1,d2,allCellMapsMax);

perf = zeros([length(gridSize),length(minPatchSize),length(maxShift),length(iterations),length(us_fac),length(mot_uf)]);
for a = 1:length(gridSize)
    disp(repmat('#',[3,100]))
    disp([num2str(a/length(gridSize)),'%'])
    for b = 1:length(minPatchSize)
        disp(repmat('%',[2,100]))
        disp([num2str(b/length(minPatchSize)),'%'])
        for c = 1:length(maxShift)
            disp(repmat('-',[1,100]))
            disp([num2str(c/length(maxShift)),'%'])
            for d = 1:length(iterations)
                for e = 1:length(us_fac)
                    for f = 1:length(mot_uf)
                        perf(a,b,c,d,e,f) = NCcostFn([gridSize(a),minPatchSize(b),maxShift(c),iterations(d),us_fac(e),mot_uf(f)],d1,d2,allCellMapsMax);                
                    end
                end
            end
        end
    end
end

% [boundTB,fvalTB] = ga(fitFunTB,length(fieldsTB),[],[],[],[],TBLB,TBUB,[],1:length(fieldsTB));
%[bound,fval] = patternsearch(fitFun,x0,[],[],[],[],LB,UB);
% [bound,fval] = particleswarm(fitFun,varN,LB,UB);