%clear all
% dataDir = 'C:\Users\David\Desktop\idibaps\caImagging\Data\';
% animalNames = {'m2028'};%,'m2029'};
%[myData,allValidCellMax,cellsPerBatch] = getData(dataDir,animalNames);
% load('C:\Users\David\Desktop\temp_CellINfo\201502282-icx\cellInfo\cellImages.mat')
%load('C:\Users\David\Desktop\temp_CellINfo\201502282-icx\cellInfo\inputMovie.mat')
%load('C:\Users\David\Desktop\temp_CellINfo\201502282-icx\cellInfo\inputSignals.mat')
%load('C:\Users\David\Desktop\temp_CellINfo\201502282-icx\cellInfo\testpeaksArray.mat')

% cellImages = permute(emAnalysisOutput.cellImages,[3,1,2]);
% inputSignals = double(emAnalysisOutput.scaledProbability);
% hinf = hdf5info([myPath,'/2015_02_28_p000_mouse2028_NULL000_turboreg_crop_dfof_downsampleTime_1.h5']);
% inputMovie = hdf5read(hinf.GroupHierarchy.Datasets);
allValidCellMax = zeros([1,size(inputSignals,1)]);

customSort = 0;

allValidCellMax(allValidCellMax==3)=0;

pred = allValidCellMax';%zeros([1,length(allValidCellMax)])';%rulePrd(myData) & myData.imMovCorr > 0.5;
%pred = rulePrdParametric(myData,gaLimits.lLimits,gaLimits.uLimits);
%pred = rulePrd(myData);
%plotconfusion(pred',allValidCellMax)

if ~customSort
    predCrop = pred(1:size(inputSignals,1));
    groundTruthCrop = allValidCellMax(1:size(inputSignals,1));
    myDataCrop = {};
%     nameFields = fields(myData);
    dispIdx = 1:length(predCrop);%find(predCrop' == 1 & groundTruthCrop == 0);
%     for j = 1:length(nameFields)
%         myDataCrop.(nameFields{j}) = myData.(nameFields{j})(dispIdx);
%     end
    options.minValConstant = -400;
    options.inputMovie = inputMovie;
    options.valid = predCrop(dispIdx);
    options.showROITrace = 0;
    signalSorterCustom(cellImages(dispIdx,:,:),inputSignals(dispIdx,:),groundTruthCrop(dispIdx),myDataCrop,'',[],'options',options)
end


if customSort
    fig = figure;

    cellImages = permute(cellImages,[2,3,1]);

    for j = 100:size(inputSignals,1)

        if pred(j)% == allValidCellMax(j)
            correct = true;
        else
            correct = false;
        end

        if pred(j) == 0 && allValidCellMax(j) == 1
            disp('--FALSE NEGATIVE--')
        elseif pred(j) == 1 && allValidCellMax(j) == 0
            disp('--FALSE POSITIVE--')
        end   

        if true%~correct
            cellShape = cellImages(:,:,j);
            threshold = mean(cellShape(:))+3*std(cellShape(:));
            cellShapeT = cellShape > threshold;

            [x,y] = find(cellShapeT');
            [perX,perY] = find(bwperim(cellShapeT'));
            xHigh = max(x);
            yHigh = max(y);
            xLow = min(x);
            yLow = min(y);

            windowSize = 10;
            if ~isempty(testpeaksArray{j}) && ~isempty([xHigh,yHigh,xLow,yLow]);

                windowFrames = bsxfun(@plus,testpeaksArray{j}',-windowSize:windowSize);
                windowFrames(windowFrames<=0)=1;
                windowFrames(windowFrames>=size(inputMovie,3))=size(inputMovie,3);

                movieEvents = zeros([size(inputMovie,1),size(inputMovie,2),size(windowFrames,2),size(windowFrames,1)]);
                for i = 1:length(testpeaksArray{j})
                   movieEvents(:,:,:,i) = inputMovie(:,:,windowFrames(i,:));
                end
                movieEventsAvg = nanmean(movieEvents,4);

                movieCrop = movieEventsAvg(yLow:yHigh,xLow:xHigh,:);
                movieMin = min(movieCrop(:));
                movieMax = max(movieCrop(:));


            %Mean signal

            inputSignal = inputSignals(j,:);
            windowTraces = reshape(inputSignal(windowFrames),size(windowFrames));
            windowTracesAvg = mean(windowTraces,1);

            %Evaluating cell

                %expratio
                meanSignalPre =  windowTracesAvg(1:windowSize+1);
                meanSignalPost =  windowTracesAvg(windowSize+1:end);

                f = {NaN,NaN};
                ind = 0;
                for meanSignal = [flip(meanSignalPre)',meanSignalPost']
                    ind = ind + 1;
                    [~,locs] = findpeaks(diff(meanSignal));
                    if isempty(locs) %meanSignal has no minima
                        [~,locs] = max(diff(meanSignal));
                        locs = flip(locs);
                    end
                    trimMeanSignal = meanSignal(1:locs(1));
                    x = 1:size(trimMeanSignal);
                    try
                        f{ind} = fit(x',trimMeanSignal,'exp1');
                    catch
                        warning('Could not fit the trace')
                        f{ind} = struct('a' , NaN, 'b', NaN);
                    end
                end

                fittedPre = f{1}.a*exp(f{1}.b*(1:0.01:size(meanSignalPre,2)));
                fittedPost = f{2}.a*exp(f{2}.b*(1:0.01:size(meanSignalPost,2)));

                expRatio = f{2}.b/f{1}.b;
            end
                %dataFeatures(1,j) = expRatio;
                %disp(['Exp ratio: ', num2str(expRatio)])
                disp(['Exp ratio: ', num2str(myData.expRatio(j))])
                disp(['Correlation: ', num2str(myData.imMovCorr(j))])
                %dataFeatures(2,j) = c2;

            %Plotting

            subplot(2,2,1)
            imagesc(cellShape(yLow:yHigh,xLow:xHigh)); 
            set(gcf,'currentch','3');
            keyIn = get(gcf,'CurrentCharacter');
            frame = -1;
            while strcmp(keyIn,'3')
                frame = frame +1;
                k = mod(frame,(windowSize*2+1))+1;
                keyIn = get(gcf,'CurrentCharacter');
                subplot(2,2,2)
                if ~isempty(testpeaksArray{j}) && ~isempty([xHigh,yHigh,xLow,yLow])
                    normalizedImg = (movieCrop(:,:,k)-movieMin)/(movieMax-movieMin);
                    imagesc(normalizedImg);
                    caxis([0 1]);
                    hold on
                    plot(1+perX-min(perX),1+perY-min(perY),'.w','markerSize',30)  
                    hold off

                subplot(2,2,3:4)
                plot(linspace(1,ceil(length(windowTracesAvg)/2),size(fittedPre,2)),flip(fittedPre),'c','linewidth',2);
                hold on
                plot(linspace(ceil(length(windowTracesAvg)/2),length(windowTracesAvg),size(fittedPost,2)),fittedPost,'c','linewidth',2);
                plot(windowTracesAvg(1:k),'b')
                xlim([1,length(windowTracesAvg)])
                ylim([0,max(windowTraces(:))])
                hold off
                end
                if correct
                    whitebg([0,0.7,0.1]);
                else
                    whitebg([0.6,0.1,0.1]);
                end
                set(subplot(2,2,3:4),'Color','White')
                if pred(j) == 0 && allValidCellMax(j) == 1
                    text(8.5,1.8,'--FALSE NEGATIVE--','fontSize',20)
                elseif pred(j) == 1 && allValidCellMax(j) == 0
                    text(8.5,1.8,'--FALSE POSITIVE--','fontSize',20)
                end
                pause(0.01)
            end
            figure(fig)
            clc
        end
    end
end