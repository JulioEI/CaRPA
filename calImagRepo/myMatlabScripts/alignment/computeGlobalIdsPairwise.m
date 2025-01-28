function [OutStruct] = computeGlobalIdsPairwise(OutStruct,coords,options,nTrials)
    % matches obj coordinates across trials, assigns them a global ID or creates a new one if no global ID is found

    % initialize the global centroid
    coordsGlobal = coords{options.trialToAlign};

    % initialize global ref matrix
    % OutStruct.globalIDs(:,1) = 1:size(coordsGlobal,1);
    OutStruct.globalIDs(:,options.trialToAlign) = 1:length(coords{options.trialToAlign});

    ignoreDistanceReplace = 1e4;

    % variables with 'local' in them are those for the current trial being compared to global
    for i=1:nTrials
        if i==options.trialToAlign
            OutStruct.globalIDs(1:length(coords{options.trialToAlign}),i) = 1:length(coords{options.trialToAlign});
            if isempty(options.trialIDs)
                OutStruct.trialIDs{i} = i;
            else
                OutStruct.trialIDs{i} = options.trialIDs{i};
            end
            continue
        end
        localCoords = coords{i};
        nGlobalCentroids = size(coordsGlobal,1);
        nLocalCentroids = size(localCoords,1);

        % obtain distance matrix of the centroid
        distanceMatrix = squareform(pdist([coordsGlobal; localCoords]));
        % avoid matching objects from the same set
        distanceMatrix(1:nGlobalCentroids, 1:nGlobalCentroids)=ignoreDistanceReplace;
        distanceMatrix(nGlobalCentroids+(1:nLocalCentroids), nGlobalCentroids+(1:nLocalCentroids))=ignoreDistanceReplace;
        % avoid matching objects to themselves
        distanceMatrix.*~diag(ones(1,nLocalCentroids+nGlobalCentroids));
        % figure(888+i);imagesc(distanceMatrix);

        % global are in rows, local in columns
        distanceMatrixCut = distanceMatrix(1:nGlobalCentroids,(nGlobalCentroids+1):end);
        % figure(999+i);imagesc(distanceMatrixCut);
        % find the minimum for the index
        [minDistances, minLocalIdx] = min(distanceMatrixCut,[],2);

        % remove duplicate minimum matching values
        minLocalIdxUnique = unique(minLocalIdx);
        for uniqueIdx = 1:length(minLocalIdxUnique)
            duplicateIdx = find(minLocalIdx==minLocalIdxUnique(uniqueIdx));
            distanceIdx = minDistances(duplicateIdx);
            [minDistDup, minLocalIdxDup] = min(distanceIdx);
            ignoreIdx = duplicateIdx;
            ignoreIdx(minLocalIdxDup) = [];
            % ignoreIdx = duplicateIdx(duplicateIdx~=minLocalIdxDup);
            display([num2str(duplicateIdx(:)') ' | ' num2str(ignoreIdx(:)') ' | ' num2str(distanceIdx(:)') '| ' num2str(minLocalIdxDup(:)')])
            minDistances(ignoreIdx) = ignoreDistanceReplace;
        end

        % find the global and local indexes for distances that meet criteria
        matchedGlobalIdx = find(minDistances<options.maxDistance);
        % only take the nearest cell
        % matchedGlobalIdx = matchedGlobalIdx(1);
        matchedLocalIdx = minLocalIdx(matchedGlobalIdx);
        % remove extra indicies that exceed the local index number
        matchedLocalIdx(matchedLocalIdx>nLocalCentroids) = [];

        % find objs from current trials that don't match
        newGlobalIdx = setdiff(1:nLocalCentroids,matchedLocalIdx);

        % remove extra indicies that aren't in matched that are beyond local index size
        % matchedLocalIdx>nLocalCentroids
        % matchedLocalIdx
        % nLocalCentroids
        % matchedLocalIdx(matchedLocalIdx>nLocalCentroids) = 0;

        % add the indicies from the current trial to the global idx matrix
        OutStruct.globalIDs(matchedGlobalIdx,i) = matchedLocalIdx;

        % extend
        OutStruct.globalIDs(end+1:end+length(newGlobalIdx),i) = newGlobalIdx;

        % add unmatched local coords to the global coordinates list
        coordsGlobal(end+1:end+length(newGlobalIdx),:) = localCoords(newGlobalIdx,:);

        % average the global and matched local coordinates to get a new global coordinate for that global obj
        coordsGlobal(matchedGlobalIdx,1) = mean([coordsGlobal(matchedGlobalIdx,1) localCoords(matchedLocalIdx,1)],2);
        coordsGlobal(matchedGlobalIdx,2) = mean([coordsGlobal(matchedGlobalIdx,2) localCoords(matchedLocalIdx,2)],2);

        if isempty(options.trialIDs)
            OutStruct.trialIDs{i} = i;
        else
            OutStruct.trialIDs{i} = options.trialIDs{i};
        end
    end

    OutStruct.coordsGlobal = coordsGlobal;
end