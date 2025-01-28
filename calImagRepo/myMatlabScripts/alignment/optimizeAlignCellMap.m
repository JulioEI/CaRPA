% %%%%%%%%%%%%%%%%%%%%%%%
% %%%LOAD THE CELLMAPS%%%
% %%%%%%%%%%%%%%%%%%%%%%%
% 
% path = 'E:\Processing\Mouse2028';
% folders = {'Mouse-2028-20150301-linear-track','Mouse-2028-20150303-linear-track','Mouse-2028-20150305-linear-track','Mouse-2028-20150307-linear-track','Mouse-2028-20150309-linear-track','Mouse-2028-20150311-linear-track','Mouse-2028-20150313-linear-track','Mouse-2028-20150315-linear-track','Mouse-2028-20150317-linear-track','Mouse-2028-20150319-linear-track','Mouse-2028-20150321-linear-track','Mouse-2028-20150323-linear-track','Mouse-2028-20150325-linear-track','Mouse-2028-20150327-linear-track'};
% 
% k = 1;
% allCellMaps = {};
% for folder = folders %one file per folder!
%     allFiles = dir(fullfile(path,folder{1}));
%     for file = {allFiles(find(~cellfun(@isempty,(regexp({allFiles.name},'emAnalysis.mat'))))).name};
%         disp(file{1});
%         decisionRegexp = [file{1}(1:end-4),'Sorted.mat'];
%         decisions = load(fullfile(path,folder{1},allFiles(find(~cellfun(@isempty,(regexp({allFiles.name},decisionRegexp))))).name));
%         emFile = load(fullfile(path,folder{1},file{1}));
%         allCellMaps{k} = emFile.emAnalysisOutput.cellImages(:,:,logical(decisions.validCellMax));
%         k = k + 1;
%     end
% end
% 
% allCellMaps = cellfun(@(x) thresholdImages(x),allCellMaps,'UniformOutput',0);
% allCellMapsMax = cellfun(@(x) max(x,[],3),allCellMaps,'UniformOutput',0);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%SET PARAMETER RANGES TO OPTIMIZE OVER%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TURBOREG
optimizeTB.smoothing = [10,100];
optimizeTB.minGain = [0,1];
optimizeTB.iterations = [1,10];

fieldsTB = fields(optimizeTB);
TBLB = zeros([1,length(fieldsTB)]);
TBUB = zeros([1,length(fieldsTB)]);
for k = 1:length(fieldsTB)
    TBLB(k) = optimizeTB.(fieldsTB{k})(1);
    TBUB(k) = optimizeTB.(fieldsTB{k})(2);
end

%%
%NORMCORR
d1 = size(allCellMapsMax{1},1);
d2 = size(allCellMapsMax{1},2);

optimizeNC.grid_sizeX = [1,d1];
optimizeNC.grid_sizeY = [1,d2];
optimizeNC.min_patch_sizeX = [1,d1];
optimizeNC.min_patch_sizeY = [1,d2];
optimizeNC.max_shift = [1,1000];
optimizeNC.iter = [1,10];
optimizeNC.overlap_preX = [1,d1];
optimizeNC.overlap_preY = [1,d2];
optimizeNC.overlap_postX = [1,d1];
optimizeNC.overlap_postY = [1,d2];
optimizeNC.us_fac = [1,100];
optimizeNC.max_dev = [1,100];
optimizeNC.mot_uf = [1,100];

fieldsNC = fields(optimizeNC);
NCLB = zeros([1,length(fieldsNC)]);
NCUB = zeros([1,length(fieldsNC)]);
for k = 1:length(fieldsNC)
    NCLB(k) = optimizeNC.(fieldsNC{k})(1);
    NCUB(k) = optimizeNC.(fieldsNC{k})(2);
end
%%
%%%%%%%%%%%%%%%%%
%%%OPTIMIZE TB%%%
%%%%%%%%%%%%%%%%%
fitFunTB = @(x) TBcostFn(x,allCellMaps,allCellMapsMax);
[boundTB,fvalTB] = ga(fitFunTB,length(fieldsTB),[],[],[],[],TBLB,TBUB,[],1:length(fieldsTB));
%[bound,fval] = patternsearch(fitFun,x0,[],[],[],[],LB,UB);
% [bound,fval] = particleswarm(fitFun,varN,LB,UB);

%%
%%%%%%%%%%%%%%%%%
%%%OPTIMIZE NC%%%
%%%%%%%%%%%%%%%%%
fitFunNC = @(x) NCcostFn(x,d1,d2,allCellMapsMax);
[boundNC,fvalNC] = ga(fitFunNC,length(fieldsNC),[],[],[],[],NCLB,NCUB,[],1:length(fieldsNC));

