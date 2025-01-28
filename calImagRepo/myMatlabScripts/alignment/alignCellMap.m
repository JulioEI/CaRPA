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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%CREATE FAKE MANUALLY SHIFTED CELLMAP%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test = allCellMapsMax{1};
% shf = 16;
% blank = zeros(size(test));
% blank(1:(end-(shf-1)),:) = test(shf:end,:);
% allCellMapsMax{end+1} = blank;
% figure;subplot(1,2,1);imagesc(allCellMapsMax{1});subplot(1,2,2);imagesc(blank);
%%
%%%%%%%%%%%%%%%%%
%%%DO TURBOREG%%%
%%%%%%%%%%%%%%%%%

ioptions.RegisType = 2;
ioptions.meanSubtract = 0;
ioptions.complementMatrix = 0;
ioptions.normalizeType = 'divideByLowpass';
ioptions.RegisType = 2;
ioptions.SmoothX = 10;
ioptions.SmoothY = 10;
ioptions.minGain = 0;
ioptions.Levels = 6;
ioptions.Lastlevels = 1;
ioptions.Epsilon = 1.1921e-07;
ioptions.zapMean = 0;
ioptions.parallel = 1;
ioptions.cropCoords = [];
ioptions.closeMatlabPool = 0;
ioptions.removeEdges = 0;
ioptions.registrationFxn = 'transfturboreg';
ioptions.turboregRotation = 1;

MprFullTB = {length(allCellMapsMax)};
for i = 1:length(allCellMapsMax)
    for j = 1:length(allCellMapsMax)
        ioptions.altMovieRegister = allCellMaps{j};
        [~,coords] = turboregMovie(cat(3,allCellMapsMax{i},allCellMapsMax{j}),'options',ioptions); 
        MprFullTB{i,j} = turboregMovie(allCellMapsMax{j},'precomputedRegistrationCooords',coords);
    end
end

%%
%%%%%%%%%%%%%%%%%%
%%%DO NORMCORRE%%%
%%%%%%%%%%%%%%%%%%

d1 = size(allCellMapsMax{1},1);
d2 = size(allCellMapsMax{1},2);

options = NoRMCorreSetParms('d1',d1,'d2',d2,'grid_size',[d1/2,d2/2,1],'bin_width',1,'min_patch_size',[5,5,1],'upd_template',0,'correct_bidir',0,'use_parallel',0,'max_shift',200,'shifts_method','linear');

MprFullNC = {length(allCellMapsMax)};
for i = 1:length(allCellMapsMax)
    template = allCellMapsMax{i};
    for j = 1:length(allCellMapsMax)
        MprFullNC{i,j} = normcorre(allCellMapsMax{j},options,template);
    end
end 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%COMPUTE CELLMAP CORRELATIONS%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

corrMatNC = zeros(length(allCellMapsMax));
corrMatTB = zeros(length(allCellMapsMax));
corrMatRaw = zeros(length(allCellMapsMax));
for i = 1:length(allCellMapsMax)
    for j = 1:length(allCellMapsMax)
        corrMatNC(i,j) = corr2(allCellMapsMax{i},MprFullNC{i,j});
        corrMatTB(i,j) = corr2(allCellMapsMax{i},MprFullTB{i,j});
        corrMatRaw(i,j) = corr2(allCellMapsMax{i},allCellMapsMax{j});
    end
end
figure;subplot(2,3,1);imagesc((abs(corrMatNC)-abs(corrMatRaw))./abs(corrMatRaw),[0,80]);title('NONLINEAR');colorbar;axis('square');subplot(2,3,2);imagesc((abs(corrMatTB)-abs(corrMatRaw))./abs(corrMatRaw),[0,80]);title('TURBOREG');colorbar;axis('square');subplot(2,3,3);imagesc(abs(corrMatRaw),[0,1]);title('RAW');axis('square');colorbar;
subplot(2,3,4);imagesc(corrMatNC,[0,1]);title('NONLINEAR');colorbar;axis('square');subplot(2,3,5);imagesc(corrMatTB,[0,1]);title('TURBOREG');colorbar;axis('square');subplot(2,3,6);imagesc(corrMatRaw,[0,1]);title('RAW');axis('square');colorbar;

%%
%%%%%%%%%%%%%%%%%%%%%
%%%VIEW ALIGNMENTS%%%
%%%%%%%%%%%%%%%%%%%%%
a = figure;
b = figure;
k = 1;
for i = 1:length(allCellMapsMax)
    template = allCellMapsMax{i};
    for j = 1:length(allCellMapsMax)
%         figure(a);subplot(2,3,1);imagesc(allCellMapsMax{j});title([num2str(j),' original; correlation: ',num2str(corrMatRaw(i,j))]);ylabel('NONLINEAR');axis square;subplot(2,3,2);imagesc(MprFullNC{i,j});title([num2str(j),' registered; correlation: ',num2str(corrMatNC(i,j))]);axis square;subplot(2,3,3);imagesc(template);title(['aligned to ', num2str(i)]);axis square;
%         figure(b);subplot(1,2,1);imagescNCellMaps({allCellMapsMax{j},MprFullNC{i,j},template},{[num2str(j),' original; correlation: ',num2str(corrMatRaw(i,j))],[num2str(j),' registered; correlation: ',num2str(corrMatNC(i,j))],['aligned to ', num2str(i)]});ylabel('NONLINEAR');axis square;       
        figure(b);subplot(1,2,1);imagescNCellMaps({MprFullNC{i,j},template},{[num2str(j),' registered; correlation: ',num2str(corrMatNC(i,j))],['aligned to ', num2str(i)]});title('NONLINEAR');axis square;       
        k = k + 1;
%         figure(a);subplot(2,3,4);imagesc(allCellMapsMax{j});title([num2str(j),' original; correlation: ',num2str(corrMatRaw(i,j))]);ylabel('TURBOREG');axis square;subplot(2,3,5);imagesc(MprFullTB{i,j});title([num2str(j),' registered; correlation: ',num2str(corrMatTB(i,j))]);axis square;subplot(2,3,6);imagesc(template);title(['aligned to ', num2str(i)]);axis square;
%         figure(b);subplot(1,2,2);imagescNCellMaps({allCellMapsMax{j},MprFullTB{i,j},template},{[num2str(j),' original; correlation: ',num2str(corrMatRaw(i,j))],[num2str(j),' registered; correlation: ',num2str(corrMatTB(i,j))],['aligned to ', num2str(i)]});ylabel('TURBOREG');axis square;               
        
        figure(b);subplot(1,2,2);imagescNCellMaps({MprFullTB{i,j},template},{[num2str(j),' registered; correlation: ',num2str(corrMatTB(i,j))],['aligned to ', num2str(i)]});title('TURBOREG');axis square;               
        
        pause;
        k = k + 1;
        
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%FIND CELL CORRESPONDENCES%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Align cells separately
templateIdx = round(length(allCellMaps)/2);

d1 = size(allCellMapsMax{1},1);
d2 = size(allCellMapsMax{1},2);

options = NoRMCorreSetParms('d1',d1,'d2',d2,'grid_size',[d1/2,d2/2,1],'bin_width',1,'min_patch_size',[5,5,1],'upd_template',0,'correct_bidir',0,'use_parallel',0,'max_shift',200,'shifts_method','linear');
template = allCellMapsMax{templateIdx};
alignedImgNC = cell([1,length(allCellMaps)]);
for j = 1:length(allCellMapsMax)
    [~,shifts] = normcorre(allCellMapsMax{j},options,template);
    for cellIdx = 1:size(allCellMaps{j},3)
        alignedImgNC{j} = cat(3,alignedImgNC{j},apply_shifts(allCellMaps{j}(:,:,cellIdx),shifts,options));
    end
end
%%
%Find cell centroids
coords = {};
for j=1:length(allCellMapsMax)
    [xCoords yCoords] = findCentroid(alignedImgNC{j});
    coords{j} = [xCoords; yCoords]';
end
%%
%Run cell correspondence finder (Biafra's)
globalIDoptions.analysisType = 'pairwise';
globalIDoptions.trialToAlign = templateIdx;
globalIDoptions.maxDistance = 5;
globalIDoptions.nCorrections = 3;
% globalIDoptions.threshold = 0.5000;
% globalIDoptions.inputSignals = {};
% globalIDoptions.additionalAlignmentImages = [];
globalIDoptions.trialIDs = [];
% globalIDoptions.runMotionCorrection = 1;
% globalIDoptions.altInputImagesToRegister = [];
% globalIDoptions.alignWithCentroids = 1;
[OutStruct] = computeGlobalIdsPairwise([],coords,globalIDoptions,length(allCellMapsMax));
globalIDs = OutStruct.globalIDs;
%%
%Look at results
assessAlignmentAcrossDays(allCellMaps,globalIDs)

figure;
for i = 1:length(alignedImgNC)
    template = alignedImgNC{i};
    for j = 1:length(alignedImgNC)
        clf;
        imagescNCellMaps({alignedImgNC{j},alignedImgNC{i}});axis square; %,{[num2str(j),' registered; correlation: ',num2str(corrMatNC(i,j))],['aligned to ', num2str(i)]});ylabel('NONLINEAR');axis square;         
        hold on  
        commonCenters = findcommonCenters(globalIDs,alignedImgNC,i,j);
        plot(commonCenters(:,:,2),commonCenters(:,:,1),'.','markersize',10)
        for k = 1:size(commonCenters,1)
            plot(squeeze(commonCenters(k,:,2)),squeeze(commonCenters(k,:,1)),'w:');  
        end
        title(['Correspondences between ', num2str(i),' and ',num2str(j), ' aligned to ', num2str(templateIdx)])
        pause;
    end
end
  