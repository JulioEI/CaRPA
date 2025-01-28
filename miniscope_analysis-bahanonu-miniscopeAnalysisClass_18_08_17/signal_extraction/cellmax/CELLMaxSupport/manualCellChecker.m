function [valid, options] = manualCellChecker(DFOF, cellTraces, cellImages, eventTimes, varargin)

%Goes through each cell and shows the CellImage, the CellTrace, the average
%burst, and the eventImage.
%Then asks user for a reply: y,n,s,m
%If the user is unsure (m), DFOF around the burst is played.
%
%valid is a vector indicating the users choice.
%   1: y, a valid cell
%   0.5: s, a cell to use for across day alignment but not within day statistics
%   0: n, not a cell
%   -1: any other input, these can be rechecked later if the user mistyped
%
%Written by Lacey Kitch and Maggie Carr Larkin, 2014
%--------------------------------------------------------------------------

% Preallocate output
nImages=size(cellImages,3);
valid = -1*ones(nImages,1);

% Set options
options.areaThresh=300;
options.skipNoEventCells=1;
options.playFullMovieClips=0;
options.displayOptions.plotTraces=1;
options.displayOptions.markCentroids=1;
options.currentInd=1;
options.cellOrder=1:nImages;        %%%%% change here to randomize order
options.framerate=4;
options.displayOptions.framerate=options.framerate*2;
options.meanEventImageCorrs=-1*ones(nImages,1);
options.eventImageCorrs=cell(nImages,1);
options.cvxHulls=[];
options.areas=[];
options.cellParams=[];
options.eventMontages=[];
options.doppList=zeros(0,2);
options=getOptions(options,varargin);
if isfield(options, 'valid')
    valid=options.valid;
end

% calculate convex hulls and centroids from images
if isempty(options.cvxHulls) || isempty(options.areas)
    [cvxHulls,areas,~] = getConvexHull(cellImages);
    options.cvxHulls=cvxHulls;
    options.areas=areas;
else
    cvxHulls=options.cvxHulls;
    areas=options.areas;
end

% get cell parameters
if isempty(options.cellParams)
    cellParams=getImageCentroids(cellImages,size(DFOF(:,:,1)));
    options.cellParams=cellParams;
else
    cellParams=options.cellParams;
end

% calculate all event montages
if isempty(options.eventMontages)
    eventMontages=cell(length(eventTimes),1);
    for cInd=1:length(eventTimes)
        if ~isempty(eventTimes{cInd})
            [eventMontages{cInd}, eventImageCorrs] = getEventMontage(eventTimes{cInd}, cellTraces(cInd,:), DFOF, cellParams(cInd,:), cellImages(:,:,cInd));
            options.meanEventImageCorrs(cInd)=mean(eventImageCorrs);
            options.eventImageCorrs{cInd}=eventImageCorrs;
        else
            options.meanEventImageCorrs(cInd)=0;
        end
    end
    options.eventMontages=eventMontages;
else
    eventMontages=options.eventMontages;
end

% get movie height and width
Height=size(DFOF,1);
Width=size(DFOF,2);

% Go through each cell
finished=0;
h=figure(1);
badCells=zeros(nImages,1);
lastDir=1;
while ~finished

     %try
        i=options.cellOrder(options.currentInd);
        
    
        if ~(areas(i)>options.areaThresh) && ~(options.skipNoEventCells && isempty(eventTimes{i}))

            figure(h)
            subplot(2,2,1); %Show the candidate cell's image
            imagesc(cellImages(:,:,i));
            axis off

            subplot(2,2,2); %Show the average spike waveform
            plot(0,0)
            if ~isempty(eventTimes{i})
                spikematrix = nan(length(eventTimes{i}),60);
                for s = 1:length(eventTimes{i})
                    if eventTimes{i}(s)>20 && eventTimes{i}(s) < size(cellTraces,2)-40
                        spikematrix(s,:) = cellTraces(i,eventTimes{i}(s)-19:eventTimes{i}(s)+40);
                    end
                end
                plot(-19:1:40,spikematrix','k',-19:1:40,nanmean(spikematrix),'r')
                try
                    set(gca,'xlim',[-19 40],'ylim',[min(spikematrix(:)) max(spikematrix(:))],'xtick',[],'ytick',[])
                catch
                    disp('There was an error in the try block...')
                end
            end
            clear s spikematrix

            subplot(2,2,[3 4]); %Plot the trace
            plot(cellTraces(i,:),'k')
            if ~isempty(eventTimes{i})
                hold on
                plot(eventTimes{i},cellTraces(i,eventTimes{i}),'r*');
                hold off
            end
            if valid(i)==-1
                currState='unknown';
            elseif valid(i)==1
                currState='y';
            elseif valid(i)==0
                currState='n';
            elseif valid(i)==3
                currState='c';
            end
            suptitle(sprintf('%d of %d, Cell? (y/n/c), Current: %s, Fwd (f), Back (b), Neighbors(l), Movie (m), Quit (q)', options.currentInd, length(options.cellOrder), currState))
            set(gca,'ylim',[min(cellTraces(i,:)) max(cellTraces(i,:))])
            set(h, 'CurrentCharacter', 'k');

            figure(2);
            imagesc(eventMontages{i})
            title(sprintf('%d events', length(eventTimes{i})))
            set(gca, 'Fontsize', 14)

            figure(h);
            waitforbuttonpress();
            reply=get(h, 'CurrentCharacter');

            if strcmpi(reply,'m')
                [~,burstorder] = sort(cellTraces(i,eventTimes{i}),'descend');

                options.displayOptions.lims=[max(1,floor(cellParams(i,2)-20)) min(Height,ceil(cellParams(i,2)+20)) max(1,floor(cellParams(i,1)-20)) min(Width,ceil(cellParams(i,1)+20))];
                movieRound = 1;

                veryCloseNeighbors = (cellParams(:,2)>(cellParams(i,2)-10) & ...
                    cellParams(:,2)<(cellParams(i,2)+10) & ...
                    cellParams(:,1)>(cellParams(i,1)-10) & ...
                    cellParams(:,1)<(cellParams(i,1)+10));

                while strcmpi(reply,'m')

                    if isempty(eventTimes{i}) || options.playFullMovieClips
                        framePadding=round(options.framerate*10);
                        movieind=mod(framePadding*(movieRound-1)+1,size(DFOF,3)-framePadding);
                    else
                        movieind = eventTimes{i}(burstorder(movieRound));
                        framePadding=round(options.framerate*4);
                    end
                    if movieind > framePadding && movieind < (size(DFOF,3)-framePadding)
                        options.displayOptions.cvxHulls=cvxHulls(veryCloseNeighbors);
                        options.displayOptions.skipSim=1;
                        options.displayOptions.specificCvxHull=cvxHulls{i};
                        options.displayOptions.writeAVI=0;
                        options.displayOptions.framerate=2*options.framerate;
                        options.displayOptions.grayscale=1;
                        %options.displayOptions.aviName='/Users/ljkitch/Desktop/CheckMovies/CellEventMovie';
                        makeResultsMovie(DFOF(:,:,movieind-framePadding:movieind+framePadding-1), cellParams(i,:), cellTraces(i,movieind-framePadding:movieind+framePadding), options.displayOptions);
                        figure(h);
                        set(h, 'CurrentCharacter', 'k');
                        waitforbuttonpress();
                        reply=get(h, 'CurrentCharacter');
                    end
                    if movieRound < length(eventTimes{i}) || isempty(eventTimes{i})
                        movieRound = movieRound+1;
                    else
                        movieRound = 1;
                    end
                end
            end
            clear burstorder movieind round
            if strcmpi(reply, 'l')
                options=displayNeighborMontages(i,cellParams,cellImages,cvxHulls,eventTimes,cellTraces,DFOF,valid,'options',options);
            elseif strcmpi(reply,'y')
                valid(i) = 1; %Valid cells
                options.currentInd=options.currentInd+1;
                lastDir=1;
            elseif strcmpi(reply,'n')
                valid(i) = 0; %Invalid
                options.currentInd=options.currentInd+1;
                lastDir=1;
            elseif strcmpi(reply, 'c')
                valid(i)=3;
                options.currentInd=options.currentInd+1;
                lastDir=1;
            elseif strcmpi(reply, 'd')
                valid(i)=5;
                options.currentInd=options.currentInd+1;
                lastDir=1;
            elseif strcmpi(reply, 'f')
                options.currentInd=options.currentInd+1;
                lastDir=1;
            elseif strcmpi(reply, 'b')
                options.currentInd=options.currentInd-1;
                lastDir=-1;
            elseif strcmpi(reply,'q')
                finished=1;
            end

        else
            valid(i)=0;
            badCells(i)=1;
            if sum(badCells)==nImages
                finished=1;
                disp('All the images are either no big or have no events.... Finishing.')
            end
            options.currentInd=options.currentInd+lastDir;
        end
        if options.currentInd>nImages
            options.currentInd=1;
        elseif options.currentInd<1
            options.currentInd=nImages;
        end
%     catch
%         finished=1;
%         disp('Error in loop, quitting checker. Might be done - check outputs.')
%     end
end
options.valid=valid;