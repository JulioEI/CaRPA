clear all
folders = {'D:\Storage\Processing\Mouse-2023'...
    ,'D:\Storage\Processing\Mouse-2025'...
    ,'E:\Processing\Mouse-2026'...
    ,'E:\Processing\Mouse2028'...
    ,'E:\Processing\Mouse2029'...
    ,'E:\Processing\Mouse2010_linear_track2'...
    ,'E:\Processing\Mouse-2019'};
names = {'m2023','m2025','m2026','m2028','m2029','m2010','m2019'};
frequencies = {'5','2','2','2','5','2','5'};
sessions = cell([1,length(names)]);
for k = 6%1:length(folders)
    PCCorr.(names{k}).PC = PCAnalysis(folders{k});
    sessions{k} = PCCorr.(names{k}).PC.chooseSessions;
end
%%
for k = 6%:length(folders)
    [~,PCCorr.(names{k}).meanVal,PCCorr.(names{k}).steVal] = PCCorr.(names{k}).PC.PFxCorr('dT',0.05,'sessions',sessions{k},'rField','spikeDeconv');
    PCCorr.(names{k}).turnover = PCCorr.(names{k}).PC.quantifyTurnoverAlignment('sessions',sessions{k});
    PCCorr.(names{k}).frequency = frequencies{k};
end

%%
%Do this for your animals before runing this script
% PCCorr.m2025.PC = PCAnalysis;
% [~,PCCorr.m2025.meanVal,PCCorr.m2025.steVal] = PCCorr.m2025.PC.PFxCorr('dT',0.05);
% PCCorr.m2025.frequency = 2;
% PCStruct = PCCorr;
% OR

% PCTurnover.m2025.PC = PCAnalysis;
% PCTurnover.m2025.meanVal = PCTurnover.m2025.PC.quantifyTurnoverAlignment;
% PCTurnover.m2025.frequency = 2;
% PCStruct = PCTurnover;
%%


%Create cumulative correlation matrices for all the animals
metrics = {'meanVal','turnover'};
figure;
for m = 1:length(metrics)
    animals = fields(PCCorr);
    for k = 1:length(animals)
        meanCorr = PCCorr.(animals{k}).(metrics{m});
        cumCorr = nan(length(meanCorr));
        for j = 1:length(meanCorr)
            cumCorr(j,1:(length(meanCorr)-j+1)) = meanCorr(j,j:end);
        end
        PCCorr.(animals{k}).cumCorr = cumCorr;
    end

    %Join the cumcorr of animals with same freq (padding nans)
    comulativeCorr5d = [];
    comulativeCorr2d = [];
    for k = 1:length(animals)
        if PCCorr.(animals{k}).frequency == '5'
            cumCorr = PCCorr.(animals{k}).cumCorr;
            prevSize = [size(comulativeCorr5d,1),size(comulativeCorr5d,2)];
            padMat = nan([(size(cumCorr,1) + prevSize(1)),max(size(cumCorr,2),prevSize(2))]);
            padMat(1:prevSize(1),1:prevSize(2)) = comulativeCorr5d;
            padMat((prevSize(1)+1):end,1:size(cumCorr,2)) = cumCorr;
            comulativeCorr5d = padMat;
        elseif PCCorr.(animals{k}).frequency == '2'
            cumCorr = PCCorr.(animals{k}).cumCorr;
            prevSize = [size(comulativeCorr2d,1),size(comulativeCorr2d,2)];
            padMat = nan([(size(cumCorr,1) + prevSize(1)),max(size(cumCorr,2),prevSize(2))]);
            padMat(1:prevSize(1),1:prevSize(2)) = comulativeCorr2d;
            padMat((prevSize(1)+1):end,1:size(cumCorr,2)) = cumCorr;
            comulativeCorr2d = padMat;
        else
            error('freq not understood')
        end
    end

    %Plot
    confidence5d = 1.96*nanstd(comulativeCorr5d)./sqrt(sum(~isnan(comulativeCorr5d)));
    confidence2d = 1.96*nanstd(comulativeCorr2d)./sqrt(sum(~isnan(comulativeCorr2d)));
    hold on;
    subplot(2,2,m*2-1)
    h = shadedErrorBar(1:5:size(comulativeCorr5d,2)*5,nanmean(comulativeCorr5d),confidence5d,'lineprops',{'-o','linewidth',2},'patchSaturation',0.2);
    h.mainLine.DisplayName = '5day';
    h = shadedErrorBar(1:2:size(comulativeCorr2d,2)*2,nanmean(comulativeCorr2d),confidence2d,'lineprops',{'-o','linewidth',2},'patchSaturation',0.2);
    h.mainLine.DisplayName = '2day';
    legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'));
    xlabel('days passed')
    ylim([0,1])
    title([metrics{m}, 'across days'])
    % ylabel('mean correlation')

    subplot(2,2,m*2)
    h = shadedErrorBar(1:size(comulativeCorr5d,2),nanmean(comulativeCorr5d),confidence5d,'lineprops',{'-o','linewidth',2},'patchSaturation',0.2);
    h.mainLine.DisplayName = '5day';
    h = shadedErrorBar(1:size(comulativeCorr2d,2),nanmean(comulativeCorr2d),confidence2d,'lineprops',{'-o','linewidth',2},'patchSaturation',0.2);
    h.mainLine.DisplayName = '2day';
    legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'));
    xlabel('sessions passed')
    ylim([0,1])
    title([metrics{m}, 'across sessions'])
end

%%
%Ttest
pVal = zeros([1,size(comulativeCorr5d,2)]);
for distanceSessions = 1:size(comulativeCorr5d,2)
    distance5d = comulativeCorr5d(:,distanceSessions);
    distance2d = comulativeCorr2d(:,distanceSessions);
    
    [~,pVal(distanceSessions)] = ttest2(distance5d,distance2d);
end

%%
%Find alignment of one animal respect to confidence quartiles
animals = fields(PCCorr);
for k = 1:length(animals)
    scores = PCCorr.(names{k}).PC.getClassifierConfidence(sessions{k});
    globalIDs = PCCorr.(names{k}).PC.alignSessions(sessions{k});
    %Look at scores per same cell
    scoreId = zeros(size(globalIDs));
    for celli = 1:size(globalIDs,1)
        for sess = 1:size(globalIDs,2)
            if globalIDs(celli,sess) ~= 0
                scoreId(celli,sess) = scores{sess}(globalIDs(celli,sess));
            end
        end
    end
    %Compute the variance when found
    scoreId(scoreId == 0) = nan;
    PCCorr.(names{k}).scoreStdPerCell = nanstd(scoreId,0,2);
    PCCorr.(names{k}).scoreMeanPerCell = nanmean(scoreId,2);
    PCCorr.(names{k}).scoreMeanPerCell(isnan(PCCorr.(names{k}).scoreMeanPerCell)) = 0;
    quartiles = [0,quantile(PCCorr.(names{k}).scoreMeanPerCell,.25),quantile(PCCorr.(names{k}).scoreMeanPerCell,.50),quantile(PCCorr.(names{k}).scoreMeanPerCell,.75),quantile(PCCorr.(names{k}).scoreMeanPerCell,1)];

    for j = 2:length(quartiles)
        quantileScoreIdx = find(PCCorr.(names{k}).scoreMeanPerCell >= quartiles(j-1) & PCCorr.(names{k}).scoreMeanPerCell < quartiles(j));
        disp(length(quantileScoreIdx))
        quartileScores = repmat({quantileScoreIdx},[1,length(scores)]);
        PCCorr.(names{k}).(['q',num2str(j-1)]) = PCCorr.(names{k}).PC.quantifyTurnoverAlignment('sessions',sessions{k},'filterPF',quartileScores);
    end
end

%%
quartilScore2D = zeros([1,(length(quartiles)-1)]);
quartilScore5D = zeros([1,(length(quartiles)-1)]);
f1 = figure;
f2 = figure;
for m = 1:(length(quartiles)-1)
    animals = fields(PCCorr);
    for k = 1:length(animals)
        meanCorr = PCCorr.(animals{k}).(['q',num2str(m)]);
        cumCorr = nan(length(meanCorr));
        for j = 1:length(meanCorr)
            cumCorr(j,1:(length(meanCorr)-j+1)) = meanCorr(j,j:end);
        end
        PCCorr.(animals{k}).cumCorr = cumCorr;
    end

    %Join the cumcorr of animals with same freq (padding nans)
    comulativeCorr5d = [];
    comulativeCorr2d = [];
    for k = 1:length(animals)
        if PCCorr.(animals{k}).frequency == '5'
            cumCorr = PCCorr.(animals{k}).cumCorr;
            prevSize = [size(comulativeCorr5d,1),size(comulativeCorr5d,2)];
            padMat = nan([(size(cumCorr,1) + prevSize(1)),max(size(cumCorr,2),prevSize(2))]);
            padMat(1:prevSize(1),1:prevSize(2)) = comulativeCorr5d;
            padMat((prevSize(1)+1):end,1:size(cumCorr,2)) = cumCorr;
            comulativeCorr5d = padMat;
        elseif PCCorr.(animals{k}).frequency == '2'
            cumCorr = PCCorr.(animals{k}).cumCorr;
            prevSize = [size(comulativeCorr2d,1),size(comulativeCorr2d,2)];
            padMat = nan([(size(cumCorr,1) + prevSize(1)),max(size(cumCorr,2),prevSize(2))]);
            padMat(1:prevSize(1),1:prevSize(2)) = comulativeCorr2d;
            padMat((prevSize(1)+1):end,1:size(cumCorr,2)) = cumCorr;
            comulativeCorr2d = padMat;
        else
            error('freq not understood')
        end
    end

    quartilScore2D(m) = mean(nanmean(comulativeCorr5d));
    quartilScore5D(m) = mean(nanmean(comulativeCorr2d));
    
    %Plot
    figure(f1)
    confidence5d = 1.96*nanstd(comulativeCorr5d)./sqrt(sum(~isnan(comulativeCorr5d)));
    confidence2d = 1.96*nanstd(comulativeCorr2d)./sqrt(sum(~isnan(comulativeCorr2d)));
    hold on;
%     subplot(ceil(sqrt(length(quartiles)-1)),ceil(sqrt(length(quartiles)-1)),m)
%     h = shadedErrorBar(1:5:size(comulativeCorr5d,2)*5,nanmean(comulativeCorr5d),confidence5d,'lineprops',{'-o','linewidth',2},'patchSaturation',0.2);
%     h.mainLine.DisplayName = '5day';
    h = shadedErrorBar(1:2:size(comulativeCorr2d,2)*2,nanmean(comulativeCorr2d),confidence2d,'lineprops',{'-o','linewidth',2},'patchSaturation',0.2);
    h.mainLine.DisplayName = ['q',num2str(m)];%2day';
    legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'));
    xlabel('days passed')
    ylim([0,1])
    title('5day')%title([['q',num2str(m)], ' across days'])
    ylabel('recurrence probability')
% 
%     figure(f2)
%     subplot(ceil(sqrt(length(quartiles)-1)),ceil(sqrt(length(quartiles)-1)),m)
%     h = shadedErrorBar(1:size(comulativeCorr5d,2),nanmean(comulativeCorr5d),confidence5d,'lineprops',{'-o','linewidth',2},'patchSaturation',0.2);
%     h.mainLine.DisplayName = '5day';
%     h = shadedErrorBar(1:size(comulativeCorr2d,2),nanmean(comulativeCorr2d),confidence2d,'lineprops',{'-o','linewidth',2},'patchSaturation',0.2);
%     h.mainLine.DisplayName = '2day';
%     legend(findobj(gca, '-regexp', 'DisplayName', '[^'']'));
%     xlabel('sessions passed')
%     ylim([0,1])
%     ylabel('recurrence probability')
%     title([['q',num2str(m)], ' across sessions'])
end
%%

%%
figure;
bar([quartilScore2D;quartilScore5D]')
legend('2day','5day')
xlabel('quartiles')
ylabel('Mean recurrence probability')
%%


