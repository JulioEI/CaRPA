function [ cellInfo ] = simpleGetCellInfo(inputSignals,testpeaksArray,cellImages,inputMovie, varargin)

% fig = figure;

nCells = size(inputSignals,1);

cellInfo.expRatio = NaN([nCells,1]);
cellInfo.imMovCorr = NaN([nCells,1]);

for j = 1:nCells
    if ~isempty(testpeaksArray{j})
        %disp(['CELL ',num2str(j)])
        cellShape = cellImages(:,:,j);
        threshold = mean(cellShape(:))+3*std(cellShape(:));
        cellShapeT = cellShape > threshold;

        [x,y] = find(cellShapeT');   
        xHigh = max(x);
        yHigh = max(y);
        xLow = min(x);
        yLow = min(y);

        windowSize = 10;
        windowFrames = bsxfun(@plus,testpeaksArray{j}',-windowSize:windowSize);
        windowFrames(windowFrames<=0)=1;
        windowFrames(windowFrames>=size(inputMovie,3))=size(inputMovie,3);

        movieEvents = zeros([size(inputMovie,1),size(inputMovie,2),size(windowFrames,2),size(windowFrames,1)]);
        for i = 1:length(testpeaksArray{j})
           movieEvents(:,:,:,i) = inputMovie(:,:,windowFrames(i,:));
        end
        movieEventsAvg = mean(movieEvents,4);

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

    %         fittedPre = f{1}.a*exp(f{1}.b*(1:0.01:size(meanSignalPre,2)));
    %         fittedPost = f{2}.a*exp(f{2}.b*(1:0.01:size(meanSignalPost,2)));

            expRatio = f{2}.b/f{1}.b;
            
            cellInfo.expRatio(j) = expRatio;
            %disp(['Exp ratio: ', num2str(expRatio)])
            
            %corr
            A = cellShape(yLow:yHigh,xLow:xHigh);
            B = (movieCrop(:,:,windowSize+1)-movieMin)/(movieMax-movieMin);

            AMean = nanmean(A(:));
            BMean = nanmean(B(:));
            c2 = nansum((A(:)-AMean).*(B(:)-BMean))./sqrt(nansum((A(:)-AMean).^2).*nansum((B(:)-BMean).^2));     
            %disp(['Correlation: ', num2str(c2)])
            cellInfo.imMovCorr(j) = c2;
    else
        %disp('No events found')
    end
%         %MI
%         mi = mutualinfo(uint8(100*A(:)),uint8(100*B(:)));
%         disp(['Mutual information: ', num2str(mi)])

%     if expRatio < 1 && c2 > 0.8
%         correct = true;
%     else
%         correct = false;
%     end        
        
%     %Plotting
%     
%     subplot(2,2,1)
%     imagesc(cellShape(yLow:yHigh,xLow:xHigh));
%   
%     set(gcf,'currentch','3');
%     keyIn = get(gcf,'CurrentCharacter');
%     frame = -1;
%     while strcmp(keyIn,'3')
%         frame = frame +1;
%         k = mod(frame,(windowSize*2+1))+1;
%         keyIn = get(gcf,'CurrentCharacter');
%         subplot(2,2,2)
%         normalizedImg = (movieCrop(:,:,k)-movieMin)/(movieMax-movieMin);
%         imagesc(normalizedImg);
%         caxis([0 1]);
%         
%         subplot(2,2,3:4)
%         plot(linspace(1,ceil(length(windowTracesAvg)/2),size(fittedPre,2)),flip(fittedPre),'c','linewidth',2);
%         hold on
%         plot(linspace(ceil(length(windowTracesAvg)/2),length(windowTracesAvg),size(fittedPost,2)),fittedPost,'c','linewidth',2);
%         plot(windowTracesAvg(1:k),'b')
%         xlim([1,length(windowTracesAvg)])
%         ylim([0,max(windowTraces(:))])
%         hold off
%         
%         if correct
%             whitebg([0,0.7,0.1]);
%         else
%             whitebg([0.6,0.1,0.1]);
%         end
%         set(subplot(2,2,3:4),'Color','White')
%         
%         pause(0.01)
%     end
%     figure(fig)
%     clc
% end

end

