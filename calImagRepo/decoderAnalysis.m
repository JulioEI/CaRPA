classdef decoderAnalysis < singleAnimalAnalysis
 %`
    properties
        dimensions = 1;
        cellPermutations = 10;
        dir = 0;
        decoder = 'SVM';
        kFolds = 10;
    end
    
    methods
        function obj = decoderAnalysis(folder)
            if nargin == 0
                folder = [];
            end
            obj@singleAnimalAnalysis(folder)
        end
        
%         function dPrime = computeDPrime(obj,varargin)
%             %Get parameters
%             param = obj.getDefaultParams({'binN','minVel','pxPerCm','dtCamera','dT','filterPF'},varargin);
%             
%             decoder = dataAnalysis.parseInput({varargin,'decoder',obj.decoder});
%             kFolds = dataAnalysis.parseInput({varargin,'kFolds',obj.kFolds});
%             
%             %Get the sessions
%             sessions = obj.chooseSessions(dataAnalysis.parseInput({varargin,'sessions',[]}))   
%             
%             %Get traces and positions
%             [pos,r] = obj.getXR('sessions',sessions,varargin{:});
%             
%             %Remove cells with non significant place fields
%             if param.filterPF
%                 [r,pos] = PFAnalysis.isPF(r,pos,param,varargin);               
%             end
%             
%             for k = 1:length(sessions)
%                 disp(repmat('-',[2,100]))
%                 disp(['Session ',num2str(k),'/',num2str(length(sessions))])
%                 [xCur,rCur,cmPerBin] = dataAnalysis.curateXandR(pos{k},r{k},param);
%                 
%                 dimensions = obj.dimensions;
%                 shuffleTrials = 0;
%                 
%                 xCur = xCur(:,1);
%                 if shuffleTrials
%                     rCur = decoderAnalysis.shuffleKeepingTrials1D(rCur,xCur);
%                 end
%                 N = length(xCur);
%                 foldSize = floor(N/kFolds);
%                 Nfold = foldSize*kFolds;
%                 errorPerFold = zeros([1,kFolds]);
% 
%                 for fold = 1:kFolds
%                     testIdx = (1+(fold-1)*foldSize):(fold*foldSize);
%                     testX = xCur(testIdx);
%                     testR = rCur(testIdx,:);
% 
%                     trainIdx = setdiff(1:Nfold,testIdx);
%                     trainX = xCur(trainIdx);
%                     trainR = rCur(trainIdx,:);
%                     
%                     %mdl = fitcecoc(trainR,trainX);
%                     %predXTest = mdl.predict(testR);
%                     for i = 1:obj.binN(1)
%                         for j = 1:obj.binN(1)
%                             CVSVMModel = fitcsvm(trainX,trainY,'Standardize',true);
%                         end
%                     end
%                     [ScoreCVSVMModel,transform] = fitSVMPosterior(mdl);
%                     [predXTest,predProb] = ScoreCVSVMModel.predict(testR);
%                     predProb = predProb(:,2);
%                     
%                     
%                 end
%                 
%             end
%             
%         end

        function cmPerBin = getCmPerBin(obj,varargin)
            %Get the sessions
            sessions = obj.chooseSessions(dataAnalysis.parseInput({varargin,'sessions',[]}));   
            %Get parameters
            param = obj.getDefaultParams({'binN','minVel','pxPerCm','dtCamera','dT','filterPF'},varargin);
            %Get traces and positions
            [pos,r] = obj.getXR('sessions',sessions,varargin{:}); 
            cmPerBin = cell([1,length(sessions)]);
             for k = 1:length(sessions)
                [~,~,~,cmPerBin{k}] = dataAnalysis.curateXandR(pos{k},r{k},param);
             end
        end
        
        function errors = decoderError(obj,varargin)
            %Decodes
            
            %Get the sessions
            sessions = obj.chooseSessions(dataAnalysis.parseInput({varargin,'sessions',[]}));   
            
            %Get parameters
            param = obj.getDefaultParams({'binN','minVel','pxPerCm','dtCamera','dT','filterPF'},varargin);
            dir = dataAnalysis.parseInput({varargin,'dir',obj.dir});
            decoder = dataAnalysis.parseInput({varargin,'decoder',obj.decoder});
            kFolds = dataAnalysis.parseInput({varargin,'kFolds',obj.kFolds});
            
            %Get traces and positions
            [pos,r] = obj.getXR('sessions',sessions,varargin{:});
            
            %Remove cells with non significant place fields
            if param.filterPF
                [r,pos] = PFAnalysis.isPF(r,pos,param,varargin);               
            end
            
            for k = 1:length(sessions)
                disp(repmat('-',[2,100]))
                disp(['Session ',num2str(k),'/',num2str(length(sessions))])
                errors(k).sessInf = sessions(k);
                cellN = dataAnalysis.parseInput({varargin,'cellN',size(r{k},2)});
                r{k} = r{k}(:,1:cellN);
                %Curate position and r
                if dir
                    [xR,xL,rR,rL] = dataAnalysis.separateLR1D(r{k},pos{k});
                    [xRCur,rRCur,~,cmPerBin] = dataAnalysis.curateXandR(xR,rR,param);
                    [xLCur,rLCur] = dataAnalysis.curateXandR(xL,rL,param);
                else
                    [xCur,rCur,~,cmPerBin] = dataAnalysis.curateXandR(pos{k},r{k},param);
                end
                decCells(k).cmPerBin = cmPerBin;
                dimensions = obj.dimensions;
                if dir
                    shuffleTrials = 1;
                    shuffledDataL = decoderAnalysis.decodePosKFoldParallel(xLCur,rLCur,kFolds,decoder,dimensions,shuffleTrials);
                    shuffledDataR = decoderAnalysis.decodePosKFoldParallel(xRCur,rRCur,kFolds,decoder,dimensions,shuffleTrials);
                    shuffleTrials = 0;
                    rawDataL = decoderAnalysis.decodePosKFoldParallel(xLCur,rLCur,kFolds,decoder,dimensions,shuffleTrials);
                    rawDataR = decoderAnalysis.decodePosKFoldParallel(xRCur,rRCur,kFolds,decoder,dimensions,shuffleTrials);
                else
                    permuteTrials = 0;
                    shuffleTrials = 1;
                    shuffledData = decoderAnalysis.decodePosKFoldParallel(xCur,rCur,kFolds,decoder,dimensions,shuffleTrials,permuteTrials);
                    shuffleTrials = 0;
                    rawData = decoderAnalysis.decodePosKFoldParallel(xCur,rCur,kFolds,decoder,dimensions,shuffleTrials,permuteTrials);
                end
                if dir
                    errors(k).rawData.L = rawDataL;
                    errors(k).rawData.R = rawDataR;
                    errors(k).shuffledData.L = shuffledDataL;
                    errors(k).shuffledData.R = shuffledDataR;
                else
                    errors(k).shuffledData = shuffledData;
                    errors(k).rawData = rawData;
                end
            end
        end
        
        function decCells = sanityCheck(obj,varargin)
            %Look at the varible when randomly permuting
            
            %Get parameters
            param = obj.getDefaultParams({'binN','minVel','pxPerCm','dtCamera','dT','filterPF'},varargin);
            
            dir = dataAnalysis.parseInput({varargin,'dir',obj.dir});
            decoder = dataAnalysis.parseInput({varargin,'decoder',obj.decoder});
            kFolds = dataAnalysis.parseInput({varargin,'kFolds',obj.kFolds});
            numEval = dataAnalysis.parseInput({varargin,'numEval',10});
            
            %Get the sessions
            sessions = obj.chooseSessions(dataAnalysis.parseInput({varargin,'sessions',[]})); 
            
            %Get traces and positions
            [pos,r] = obj.getXR('sessions',sessions,varargin{:});
            
            %Remove cells with non significant place fields
            if param.filterPF
                [r,pos] = PFAnalysis.isPF(r,pos,param,varargin);               
            end
            
            decCells = struct('sessInf', cell(1, length(sessions)),...
                'cellVect', cell(1, length(sessions)),...
                'cmPerBin',cell(1, length(sessions)),...
                'cellPerm',cell(1, length(sessions)),...
                'shuffledData',cell(1, length(sessions)),...
                'rawData',cell(1, length(sessions)));
            
            for k = 1:length(sessions)
                disp(repmat('-',[2,100]))
                disp(['Session ',num2str(k),'/',num2str(length(sessions))])
                decCells(k).sessInf = sessions(k);
                totalCells = size(r{k},2);
                
                %Create x axis
                cellVect = round(linspace(1,totalCells,numEval));
                decCells(k).cellVect = cellVect;
            
                %Curate position and r
                if dir
                    [xR,xL,rR,rL] = dataAnalysis.separateLR1D(r{k},pos{k});
                    [xRCur,rRCur,~,cmPerBin] = dataAnalysis.curateXandR(xR,rR,param);
                    [xLCur,rLCur] = dataAnalysis.curateXandR(xL,rL,param);
                else
                    [xCur,rCur,~,cmPerBin] = dataAnalysis.curateXandR(pos{k},r{k},param);
                end    
                decCells(k).cmPerBin = cmPerBin;
                
                %Compute the cell idx of each x
                permIdx = cell([length(cellVect),1]);
                selectedCells = randperm(totalCells);
                for cellN = 1:length(cellVect)
                    permIdx{cellN} = [1,selectedCells(1:cellVect(cellN))];
                end
                decCells(k).cellPerm = permIdx; 
                
                maxPerm = max(cellfun(@(x) size(x,1),permIdx));
                if dir
                    rawDataL = nan([maxPerm,length(cellVect),kFolds]);
                    rawDataR = nan([maxPerm,length(cellVect),kFolds]);
                    shuffledDataL = nan([maxPerm,length(cellVect),kFolds]);
                    shuffledDataR = nan([maxPerm,length(cellVect),kFolds]);
                else
                    rawData = nan([maxPerm,length(cellVect),kFolds]);
                    shuffledData = nan([maxPerm,length(cellVect),kFolds]);
                end
                dimensions = obj.dimensions;
                
                for cellN = 1:length(cellVect)
                    disp(repmat('.',[1,100]))
                    disp(['Looking at ',num2str(cellN),'/',num2str(length(cellVect)), ' (', num2str(cellVect(cellN)),') cells'])
                    xTimer = tic;
                    
                    thisCellNPerms = permIdx{cellN};
                    %textprogressbar(['Starting permutations... ']);
                    shuffleTrials = 1;
                    for perm = 1:size(thisCellNPerms,1)
                        %textprogressbar(100*perm./size(thisCellNPerms,1));
                        if dir
                            rLimR = rRCur(:,thisCellNPerms(perm,:));
                            rLimL = rLCur(:,thisCellNPerms(perm,:));
                            permuteTrials = 1;
                            shuffledDataL(perm,cellN,:) = decoderAnalysis.decodePosKFoldParallel(xLCur,rLimL,kFolds,decoder,dimensions,shuffleTrials,permuteTrials);
                            shuffledDataR(perm,cellN,:) = decoderAnalysis.decodePosKFoldParallel(xRCur,rLimR,kFolds,decoder,dimensions,shuffleTrials,permuteTrials);
                            permuteTrials = 0;
                            rawDataL(perm,cellN,:) = decoderAnalysis.decodePosKFoldParallel(xLCur,rLimL,kFolds,decoder,dimensions,shuffleTrials,permuteTrials);
                            rawDataR(perm,cellN,:) = decoderAnalysis.decodePosKFoldParallel(xRCur,rLimR,kFolds,decoder,dimensions,shuffleTrials,permuteTrials);
                        else
                            rLim = rCur(:,thisCellNPerms(perm,:));
                            permuteTrials = 1;
                            shuffledData(perm,cellN,:) = decoderAnalysis.decodePosKFoldParallel(xCur,rLim,kFolds,decoder,dimensions,shuffleTrials,permuteTrials);
                            permuteTrials = 0;
                            rawData(perm,cellN,:) = decoderAnalysis.decodePosKFoldParallel(xCur,rLim,kFolds,decoder,dimensions,shuffleTrials,permuteTrials);
                        end                     
                    end
                    %textprogressbar(['done!']);
                    disp(['...this run took ',datestr(toc(xTimer)/86400, 'HH:MM:SS')]);clear xTimer;
                end
                
                if dir
                    decCells(k).rawData.L = rawDataL;
                    decCells(k).rawData.R = rawDataR;
                    decCells(k).shuffledData.L = shuffledDataL;
                    decCells(k).shuffledData.R = shuffledDataR;
                else
                    decCells(k).shuffledData = shuffledData;
                    decCells(k).rawData = rawData;
                end  
            end
        end
        
        function decCells = decoderAcrossCellsAdaptative(obj,varargin)
            %Decoder error when increasing cell number with adaptative
            %permutation numbers (min 'minPerm' permutations)
            
            %Get parameters
            param = obj.getDefaultParams({'binN','minVel','pxPerCm','dtCamera','dT','filterPF'},varargin);
            
            dir = dataAnalysis.parseInput({varargin,'dir',obj.dir});
            decoder = dataAnalysis.parseInput({varargin,'decoder',obj.decoder});
            kFolds = dataAnalysis.parseInput({varargin,'kFolds',obj.kFolds});
            numEval = dataAnalysis.parseInput({varargin,'numEval',10});%30});
            minPerm = dataAnalysis.parseInput({varargin,'minPerm',1});%5});
            
            %Get the sessions
            sessions = obj.chooseSessions(dataAnalysis.parseInput({varargin,'sessions',[]})); 
            
            %Get traces and positions
            [pos,r] = obj.getXR('sessions',sessions,varargin{:});
            
            %Remove cells with non significant place fields
            if param.filterPF
                [r,pos] = PFAnalysis.isPF(r,pos,param,varargin);               
            end
            
            decCells = struct('sessInf', cell(1, length(sessions)),...
                'cellVect', cell(1, length(sessions)),...
                'cmPerBin',cell(1, length(sessions)),...
                'cellPerm',cell(1, length(sessions)),...
                'shuffledData',cell(1, length(sessions)),...
                'rawData',cell(1, length(sessions)));
            
            for k = 1:length(sessions)
                disp(repmat('-',[2,100]))
                disp(['Session ',num2str(k),'/',num2str(length(sessions))])
                decCells(k).sessInf = sessions(k);
                totalCells = size(r{k},2);
                
                %Create x axis
                cellVect = round(linspace(1,totalCells,numEval));
                decCells(k).cellVect = cellVect;
            
                %Curate position and r
                if dir
                    [xR,xL,rR,rL] = dataAnalysis.separateLR1D(r{k},pos{k});
                    [xRCur,rRCur,~,cmPerBin] = dataAnalysis.curateXandR(xR,rR,param);
                    [xLCur,rLCur] = dataAnalysis.curateXandR(xL,rL,param);
                else
                    [xCur,rCur,~,cmPerBin] = dataAnalysis.curateXandR(pos{k},r{k},param);
                end    
                decCells(k).cmPerBin = cmPerBin;
                
                %For each x compute the number of permutations
                permVect = floor(totalCells ./ cellVect);
                
                %Compute the cell idx of each x
                cellPerm = randperm(totalCells);
                permIdx = cell([length(cellVect),1]);
                for cellN = 1:length(cellVect)
                    if 0%permVect(cellN) >= minPerm
                        %Take independent subsets of cells
                        groupIdx = 0:cellVect(cellN):totalCells;
                        permIdxMat = zeros([permVect(cellN),cellVect(cellN)]);
                        for j = 2:(length(groupIdx))
                           permIdxMat(j-1,:) =  cellPerm((1+groupIdx(j-1)):groupIdx(j));
                        end
                        permIdx{cellN} = permIdxMat;
                    else
                        %Take overlaping subsets of cells
                        permIdx{cellN} = reshape(randsample(cellPerm,cellVect(cellN)*minPerm,true)',[minPerm,cellVect(cellN)]);
                    end
                end
                decCells(k).cellPerm = permIdx; 
                
                maxPerm = max(cellfun(@(x) size(x,1),permIdx));
                if dir
                    rawDataL = nan([maxPerm,length(cellVect),kFolds]);
                    rawDataR = nan([maxPerm,length(cellVect),kFolds]);
                    shuffledDataL = nan([maxPerm,length(cellVect),kFolds]);
                    shuffledDataR = nan([maxPerm,length(cellVect),kFolds]);
                else
                    rawData = nan([maxPerm,length(cellVect),kFolds]);
                    shuffledData = nan([maxPerm,length(cellVect),kFolds]);
                    %permData = nan([maxPerm,length(cellVect),kFolds]);
                end
                dimensions = obj.dimensions;
                
                for cellN = 1:length(cellVect)
                    disp(repmat('.',[1,100]))
                    disp(['Looking at ',num2str(cellN),'/',num2str(length(cellVect)), ' (', num2str(cellVect(cellN)),') cells'])
                    xTimer = tic;
                    
                    thisCellNPerms = permIdx{cellN};
                    textprogressbar(['Starting permutations... ']);
                    for perm = 1:size(thisCellNPerms,1)
                        textprogressbar(100*perm./size(thisCellNPerms,1));
                        if dir
                            error('oudated')
%                             rLimR = rRCur(:,thisCellNPerms(perm,:));
%                             rLimL = rLCur(:,thisCellNPerms(perm,:));
%                             
%                             shuffleTrials = 1;
%                             shuffledDataL(perm,cellN,:) = decoderAnalysis.decodePosKFoldParallel(xLCur,rLimL,kFolds,decoder,dimensions,shuffleTrials);
%                             shuffledDataR(perm,cellN,:) = decoderAnalysis.decodePosKFoldParallel(xRCur,rLimR,kFolds,decoder,dimensions,shuffleTrials);
%                             shuffleTrials = 0;
%                             rawDataL(perm,cellN,:) = decoderAnalysis.decodePosKFoldParallel(xLCur,rLimL,kFolds,decoder,dimensions,shuffleTrials);
%                             rawDataR(perm,cellN,:) = decoderAnalysis.decodePosKFoldParallel(xRCur,rLimR,kFolds,decoder,dimensions,shuffleTrials);
                        else
                            rLim = rCur(:,thisCellNPerms(perm,:));
                            
                            permuteTrials = 0;
                            shuffleTrials = 1;
                            shuffledData(perm,cellN,:) = decoderAnalysis.decodePosKFoldParallel(xCur,rLim,kFolds,decoder,dimensions,shuffleTrials,permuteTrials);
                            shuffleTrials = 0;
                            rawData(perm,cellN,:) = decoderAnalysis.decodePosKFoldParallel(xCur,rLim,kFolds,decoder,dimensions,shuffleTrials,permuteTrials);
                            %permuteTrials = 1;
                            %permData(perm,cellN,:) = decoderAnalysis.decodePosKFoldParallel(xCur,rLim,kFolds,decoder,dimensions,shuffleTrials,permuteTrials);
                        end                     
                    end
                    textprogressbar(['done!']);
                    disp(['...this run took ',datestr(toc(xTimer)/86400, 'HH:MM:SS')]);clear xTimer;
                end
                
                if dir
                    decCells(k).rawData.L = rawDataL;
                    decCells(k).rawData.R = rawDataR;
                    decCells(k).shuffledData.L = shuffledDataL;
                    decCells(k).shuffledData.R = shuffledDataR;
                else
                    decCells(k).shuffledData = shuffledData;
                    decCells(k).rawData = rawData;
                    %decCells(k).permData = permData;
                end
                clear xCur permVect perm k rCur rLim shuffledData rawData cmPerBin totalCells cellVect permIdx selectedCells maxPerm dimensions thisCellNPerms shuffleTrials permuteTrials
                disp(['Memory in use: ',num2str(monitor_memory_whos)])
            end
        end
        
        function decCells = decoderAcrossCells(obj,varargin)
            error('outdated')
            %Get parameters
            param = obj.getDefaultParams({'binN','minVel','pxPerCm','dtCamera','dT','filterPF'},varargin);
            
            cellPermutations = dataAnalysis.parseInput({varargin,'cellPermutations',obj.cellPermutations});
            dir = dataAnalysis.parseInput({varargin,'dir',obj.dir});
            decoder = dataAnalysis.parseInput({varargin,'decoder',obj.decoder});
            kFolds = dataAnalysis.parseInput({varargin,'kFolds',obj.kFolds});
            numEval = dataAnalysis.parseInput({varargin,'numEval',30});

            %Get the sessions
            sessions = obj.chooseSessions(dataAnalysis.parseInput({varargin,'sessions',[]}))   
            
            %Get traces and positions
            [pos,r] = obj.getXR('sessions',sessions,varargin{:});
            
            %Remove cells with non significant place fields
            if param.filterPF
                [r,pos] = PFAnalysis.isPF(r,pos,param,varargin);               
            end
            
            for k = 1:length(sessions)
                disp(repmat('-',[2,100]))
                disp(['Session ',num2str(k),'/',num2str(length(sessions))])
                decCells(k).sessInf = sessions(k);
                decCells(k).cellPerm = [];
                
                
                %Search for the optimal way to select num of cells in x axis
%                 cellVect = (1:size(r{k},2)).^log2(length(1:size(r{k},2)))./max((1:size(r{k},2)).^log2(length(1:size(r{k},2)))).*size(r{k},2);
%                 cellVect = round(cellVect(round(linspace(1,length(cellVect),floor(length(cellVect)*0.1)))));
%                 cellVect = cellVect + (1:length(cellVect));
%                 cellVect(cellVect>size(r{k},2)) = size(r{k},2);
                %Add another point 10 cells before to see if it plateaus
                %cellVect = [cellVect(1:end-1),cellVect(end)-10,cellVect(end)];
                
                evalVectLin = round(linspace(1,size(r{k},2),numEval));
                cellVect = evalVectLin;
                
                decCells(k).cellVect = cellVect;
                %figure;plot(cellVect,'o-')
                
%                     textprogressbar(['Iterating over cells... ']);

                if dir
                    [xR,xL,rR,rL] = dataAnalysis.separateLR1D(r{k},pos{k});
                    [xRCur,rRCur,~,cmPerBin] = dataAnalysis.curateXandR(xR,rR,param);
                    [xLCur,rLCur] = dataAnalysis.curateXandR(xL,rL,param);
                else
                    [xCur,rCur,~,cmPerBin] = dataAnalysis.curateXandR(pos{k},r{k},param);
                end    
                decCells(k).cmPerBin = cmPerBin;
                                    
                cellPermMat = [];
                sizeRK2 = size(r{k},2);
                if dir
                    rLimL = cell([cellPermutations,length(cellVect)]);
                    rLimR = cell([cellPermutations,length(cellVect)]);
                else
                    rLim = cell([cellPermutations,length(cellVect)]);
                end
                for perm = 1:cellPermutations
                    cellPerm = randperm(sizeRK2);
                    cellPermMat = [cellPermMat ; cellPerm];
                    if dir                       
                       for cellN = 1:length(cellVect)
                            rLimL{perm,cellN} = rRCur(:,cellPerm(1:cellVect(cellN)));
                            rLimR{perm,cellN} = rLCur(:,cellPerm(1:cellVect(cellN)));
                       end
                    else
                        for cellN = 1:length(cellVect)
                            rLim{perm,cellN} = rCur(:,cellPerm(1:cellVect(cellN)));
                        end
                    end
                end
                decCells(k).cellPerm = cellPermMat;   
                
                if dir
                    rawDataL = zeros([cellPermutations,length(cellVect),kFolds]);
                    rawDataR = zeros([cellPermutations,length(cellVect),kFolds]);
                    shuffledDataL = zeros([cellPermutations,length(cellVect),kFolds]);
                    shuffledDataR = zeros([cellPermutations,length(cellVect),kFolds]);
                else
                    rawData = zeros([cellPermutations,length(cellVect),kFolds]);
                    shuffledData = zeros([cellPermutations,length(cellVect),kFolds]);
                end
                dimensions = obj.dimensions;
                for perm = 1:cellPermutations
                    disp(repmat('.',[1,100]))
                    disp(['Permutation ',num2str(perm),'/',num2str(cellPermutations)])
                    permTimer = tic;
                    %textprogressbar(100*cellN/length(cellVect));
                    if dir
                        parfor cellN = 1:length(cellVect)
                            shuffleTrials = 1;
                            shuffledDataL(perm,cellN,:) = decoderAnalysis.decodePosKFold(xLCur,rLimL{perm,cellN},kFolds,decoder,dimensions,shuffleTrials);
                            shuffledDataR(perm,cellN,:) = decoderAnalysis.decodePosKFold(xRCur,rLimR{perm,cellN},kFolds,decoder,dimensions,shuffleTrials);
                            shuffleTrials = 0;
                            rawDataL(perm,cellN,:) = decoderAnalysis.decodePosKFold(xLCur,rLimL{perm,cellN},kFolds,decoder,dimensions,shuffleTrials);
                            rawDataR(perm,cellN,:) = decoderAnalysis.decodePosKFold(xRCur,rLimR{perm,cellN},kFolds,decoder,dimensions,shuffleTrials);
                        end
                    else
                        parfor cellN = 1:length(cellVect)
                            shuffleTrials = 1;
                            shuffledData(perm,cellN,:) = decoderAnalysis.decodePosKFold(xCur,rLim{perm,cellN},kFolds,decoder,dimensions,shuffleTrials);
                            shuffleTrials = 0;
                            rawData(perm,cellN,:) = decoderAnalysis.decodePosKFold(xCur,rLim{perm,cellN},kFolds,decoder,dimensions,shuffleTrials);
                        end
                    end
                    %textprogressbar(['... done']);
                    disp(['...this permutation took ',datestr(toc(permTimer)/86400, 'HH:MM:SS')]);clear permTimer;
                end                
                
                if dir
                    decCells(k).rawData.L = rawDataL;
                    decCells(k).rawData.R = rawDataR;
                    decCells(k).shuffledData.L = shuffledDataL;
                    decCells(k).shuffledData.R = shuffledDataR;
                else
                    decCells(k).shuffledData = shuffledData;
                    decCells(k).rawData = rawData;
                end  
            end
        end
    end
    
    methods (Static)
        
        function [newR,newX] = shuffleSegments(r,x,shuffleTrials)
            minSegLen = 3;
            
            segments = {1};
            dX = [1;diff(x)];
            segmentIdx = 1;
            dir = find(dX,1);
            for k = 2:length(dX)
                if dir == -dX(k)&& length(segments{segmentIdx}) > minSegLen%%sign(dX(k-1)) ~= sign(dX(k)) && dX(k) ~= 0 && dX(k-1) ~= 0 
                    dir = -dir;
                    segmentIdx = segmentIdx + 1;
                    segments{segmentIdx} = [];
                    %Put half the last bin points in this segment
                    thisIdx = segments{segmentIdx-1}(x(segments{segmentIdx-1})==x(segments{segmentIdx-1}(end)));
                    pivotPoint = round(length(thisIdx)/2);
                    segments{segmentIdx-1} = segments{segmentIdx-1}(1:end-pivotPoint);
                    segments{segmentIdx} = thisIdx((1+pivotPoint):end);
                end
                segments{segmentIdx} = [segments{segmentIdx},k];
            end
%             figure;hold on;
%             plot(x,'k--')
%             for seg = 1:length(segments)
%                 plot(segments{seg},x(segments{seg}),'-o');
%             end
%             yticks(1:max(x));grid on
            
            %Take only one ponint at each bin, do the same for raw, and shufle by chaning whole segments
            segmentsSingle = cell(size(segments));
            binN = max(x);
            for seg = 1:length(segments)
                xSeg = x(segments{seg});
                segmentsSingle{seg} = zeros([binN,1]);
                for k = 1:binN %Bins go from 1 to N! This also flips all the segments in the same dir!
                    [singleIdx] = find(xSeg==k);
                    if ~isempty(singleIdx)
                        segmentsSingle{seg}(k) = segments{seg}(singleIdx(round(length(singleIdx)/2)));
                    else
                        segmentsSingle{seg}(k) = nan;
                    end
                end
            end
            badSegments = isnan(cellfun(@prod,segmentsSingle));
            segmentsSingle(badSegments) = [];
            
%             figure;hold on;
%             plot(x,'k--')
%             for seg = 1:length(segmentsSingle)
%                 plot(segmentsSingle{seg},x(segmentsSingle{seg}),'-o');
%             end
%             yticks(1:max(x));grid on

            %Swap segments
            finalIdx = cat(1,segmentsSingle{:});
            if shuffleTrials
                newX = x(finalIdx);%This wont change when shuffling
                disp([num2str(100*length(finalIdx)/length(x)),'% of data keept in order to shuffle'])                
                rSeg = cell(size(segmentsSingle));
                for celli = 1:size(r,2)
                    randSeg = randperm(length(segmentsSingle));
                    for seg = 1:length(segmentsSingle)
                        rSeg{seg} = [rSeg{seg},r(segmentsSingle{randSeg(seg)},celli)];
                    end            
                end
                newR = cat(1,rSeg{:});
            else
                disp([num2str(100*length(finalIdx)/length(x)),'% of data keept to make a fair comparison with shuffle'])
                newR = r(finalIdx,:);
                newX = x(finalIdx);
            end

        end
        
        function shuffledR = shuffleKeepingTrials1D(r,x)
            %x should be binned
            %Find timeframes in the spatial bins
            binN = max(x);
            timeFramesInBin = cell([1,binN]);
            for k = 1:binN
                timeFramesInBin{k} = find(x==k);
            end
            %Shufle the timeframes for each neuron
            shuffledR = r;
            for j = 1:size(r,2)
                for k = 1:binN
                    originalX = r(timeFramesInBin{k},j);
                    shuffledR(timeFramesInBin{k},j) = originalX(randperm(length(originalX)));
                end
            end
        end
        
%         function errorPerFold = decodePosKFold(x,r,kFolds,decoder,dimensions,shuffleTrials)
%             error('deprecated')
%             switch dimensions
%             case 1
%                 x = x(:,1);
%                 
%                 %Randomly permute
%                 randPermIdx = randperm(length(x));
%                 x = x(randPermIdx);
%                 r = r(randPermIdx,:);
%                 
%                 if shuffleTrials
%                     r = decoderAnalysis.shuffleKeepingTrials1D(r,x);
%                 end
%                 N = length(x);
%                 foldSize = floor(N/kFolds);
%                 Nfold = foldSize*kFolds;
%                 errorPerFold = zeros([1,kFolds]);
% 
%                 for k = 1:kFolds
%                     testIdx = (1+(k-1)*foldSize):(k*foldSize);
%                     testX = x(testIdx);
%                     testR = r(testIdx,:);
% 
%                     trainIdx = setdiff(1:Nfold,testIdx);
%                     trainX = x(trainIdx);
%                     trainR = r(trainIdx,:);
% 
%                     switch decoder
%                     
%                     case 'poissonMixtureDecoder'
%                         predXTest = poissonMixtureDecoder(trainR,trainX,testR);
%                         errorPerFold(k) = mean(sqrt((predXTest-testX).^2));
%                         
%                     case 'NBnormal'
%                         try
%                             mdl = fitcnb(trainR,trainX);
%                             predXTest = mdl.predict(testR);
%                         catch
%                             warning('Could not fit NB')
%                             predXTest = randi([min(x) max(x)],size(testX)); %Cannot fit model if data has 0 variance, yield random predictions
%                         end
%                         errorPerFold(k) = mean(sqrt((predXTest-testX).^2));
%                         
%                     case 'NBmvmn' 
%                         mdl = fitcnb(trainR,trainX,'DistributionNames','mvmn');
%                         predXTest = mdl.predict(testR);
%                         errorPerFold(k) = mean(sqrt((predXTest-testX).^2));
%                         
%                     case 'NN'
%                         mdl = fitcknn(trainR,trainX);
%                         predXTest = mdl.predict(testR);
%                         errorPerFold(k) = mean(sqrt((predXTest-testX).^2));
%                         
%                     case 'SVM'
% %                         predProb = zeros([length(testR),max(x)]);
% %                         for bin = 1:max(x)
% %                             frameTrainInBin = trainX == bin;                           
% %                             trainBinary = zeros(size(trainX));
% %                             trainBinary(frameTrainInBin) = 1;
% %                             
% %                             CVSVMModel = fitcsvm(trainR,trainBinary,'Standardize',true);
% %                             ScoreCVSVMModel = fitSVMPosterior(CVSVMModel);
% %                             [~,predProbBin] = ScoreCVSVMModel.predict(testR);
% %                             predProb(:,bin) = predProbBin(:,2);
% % 
% %                         end
% %                         predXTest = zeros([size(predProb,1),1]);
% %                         predMat = (predProb == max(predProb,[],2));
% %                         for row = 1:size(predProb,1)
% %                             predXTest(row) = find(predMat(row,:));
% %                         end
%                         
%                         mdl = fitcecoc(trainR,trainX); %Other way?
%                         predXTest = mdl.predict(testR);  
%                         
%                         errorPerFold(k) = mean(sqrt((predXTest-testX).^2));
%                     end
%                 end
%             case 2
%                 error('2 dimensional decoder not implemented yet')
%             end 
%         end
        
        function errorPerFold = decodePosKFoldParallel(x,r,kFolds,decoder,dimensions,shuffleTrials,permuteTrials)
            switch dimensions
            case 1
                x = x(:,1);
%                 
%                 if shuffleTrials
%                     r = decoderAnalysis.shuffleKeepingTrials1D(r,x);
%                 end
%                 %[r,x] = decoderAnalysis.shuffleSegments(r,x,shuffleTrials);
%                 
%                 if permuteTrials
%                     randPermIdx = randperm(length(x));
%                     %Randomly permute
%                     x = x(randPermIdx);
%                     r = r(randPermIdx,:);
%                 end
                
                N = length(x);
                foldSize = floor(N/kFolds);
                Nfold = foldSize*kFolds;
                errorPerFold = zeros([1,kFolds]);
                
                for k = 1:kFolds
                    testIdx = (1+(k-1)*foldSize):(k*foldSize);
                    testX = x(testIdx);
                    testR = r(testIdx,:);

                    trainIdx = setdiff(1:Nfold,testIdx);
                    trainX = x(trainIdx);
                    trainR = r(trainIdx,:);
                    
                    if shuffleTrials
                        testR = decoderAnalysis.shuffleKeepingTrials1D(testR,testX);
                        trainR = decoderAnalysis.shuffleKeepingTrials1D(trainR,trainX);
                    end
                    
                    if permuteTrials
                        randPermIdxTest = randperm(length(testX));
                        randPermIdxTrain = randperm(length(trainX));
                        %Randomly permute
                        testX = testX(randPermIdxTest);
                        testR = testR(randPermIdxTest,:);
                        trainX = trainX(randPermIdxTrain);
                        trainR = trainR(randPermIdxTrain,:);                        
                    end

                    switch decoder
                    
                    case 'poissonMixtureDecoder'
                        predXTest = poissonMixtureDecoder(trainR,trainX,testR);
                        errorPerFold(k) = mean(sqrt((predXTest-testX).^2));
                        
%                     case 'NBnormal'
%                         try
%                             mdl = fitcnb(trainR,trainX);
%                             predXTest = mdl.predict(testR);
%                         catch
%                             warning('Could not fit NB')
%                             predXTest = randi([min(x) max(x)],size(testX)); %Cannot fit model if data has 0 variance, yield random predictions
%                         end
%                         errorPerFold(k) = mean(sqrt((predXTest-testX).^2));
                        
                    case 'NBmvmn' 
                        mdl = fitcnb(trainR,trainX,'DistributionNames','mvmn');
                        predXTest = mdl.predict(testR);
                        errorPerFold(k) = mean(sqrt((predXTest-testX).^2));
                        
                    case 'NN'
                        mdl = fitcknn(trainR,trainX);
                        errorPerFold(k) = mean(sqrt((mdl.predict(testR)-testX).^2));
                        
                    case 'SVM'                        
                        mdl = fitcecoc(trainR,trainX);
                        errorPerFold(k) = mean(sqrt((mdl.predict(testR)-testX).^2));
                    end
                end
            case 2
                error('2 dimensional decoder not implemented yet')
            end 
        end
        
    end
    
end

