function options=displayNeighborMontages(currentInd,cellParams,cellImages,cvxHulls,eventTimes,cellTraces,DFOF,valid,varargin)

options.maxDist=10;
options.eventMontages=[];
options.doppList=zeros(0,2);
options=getOptions(options,varargin);

if isempty(options.eventMontages)
    disp('Calculating montage...')
    [mainCellEventMontage,~] = getEventMontage(eventTimes{currentInd}, cellTraces(currentInd,:), DFOF, cellParams(currentInd,:), cellImages(:,:,currentInd));
else
    mainCellEventMontage=options.eventMontages{currentInd};
end

figure(470)
[~,maxEventInd]=max(cellTraces(currentInd,eventTimes{currentInd}));
for cInd=1:length(eventTimes)
    if cInd~=currentInd && norm(cellParams(cInd,:)-cellParams(currentInd,:))<=options.maxDist
        numSharedEvents=length(intersect(eventTimes{currentInd},eventTimes{cInd})) + ...
            length(intersect(eventTimes{currentInd}+1,eventTimes{cInd})) + ...
            length(intersect(eventTimes{currentInd}-1,eventTimes{cInd}));
        if numSharedEvents>0
        
            h=figure(470);
            subplot(131)
            imagesc(DFOF(:,:,eventTimes{currentInd}(maxEventInd))); hold on;
            plot(cvxHulls{currentInd}(:,1),cvxHulls{currentInd}(:,2),'m')
            plot(cvxHulls{cInd}(:,1),cvxHulls{cInd}(:,2),'w')
            set(gca, 'Fontsize', 14)
            title('pink: original, white: neighbor')
            
            subplot(132)
            imagesc(mainCellEventMontage)
            set(gca, 'Fontsize', 14)
            title(sprintf('Original montage: %d events', length(eventTimes{currentInd})))
            
            subplot(133)
            if isempty(options.eventMontages)
                disp('Calculating neighbor montage')
                [eventMontage,~] = getEventMontage(eventTimes{cInd}, cellTraces(cInd,:), DFOF, cellParams(cInd,:), cellImages(:,:,cInd));
            else
                eventMontage=options.eventMontages{cInd};
            end
            imagesc(eventMontage)
            set(gca, 'Fontsize', 14)
            switch valid(cInd)
                case 0
                    cellLabel='n';
                case 0.5
                    cellLabel='s';
                case 1
                    cellLabel='y';
                case 2
                    cellLabel='c';
                otherwise
                    cellLabel='unknown';
            end
            isDopp=0;
            for dInd=1:size(options.doppList,1)
                if sum(options.doppList(dInd,:)==[cInd currentInd])==2 ||...
                       sum(options.doppList(dInd,:)==[currentInd cInd])==2
                   isDopp=1;
                end
            end
            if isDopp
                title('DOPPELGANGER PAIR MARKED')
            else
                title(sprintf('Neighbor montage: %s, %d shared events, %d total', cellLabel, numSharedEvents, length(eventTimes{cInd})))
            end
            set(h, 'CurrentCharacter', 'k');
            k=waitforbuttonpress;
            
            while k==0
                k=waitforbuttonpress;
            end
            
            reply=get(h, 'CurrentCharacter');
            
            if strcmpi(reply, 'd')
                options.doppList(end+1,:)=[currentInd, cInd];
            end
            
        end
    end
end