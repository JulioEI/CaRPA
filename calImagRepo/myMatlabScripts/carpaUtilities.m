classdef carpaUtilities
        
    methods(Static)
              
        function [movieDS,movieDFOF,movieRegFilt,movieReg] = preprocessMovie(movie,varargin)
            movie = carpaUtilities.readHdf5(movie);
            
            freqLow = carpaUtilities.parseInput({varargin,'freqLow',1});
            freqHigh = carpaUtilities.parseInput({varargin,'freqHigh',4});
            filterReg = carpaUtilities.parseInput({varargin,'filterReg','divideByLowpass'});
%             filterDFOF = carpaUtilities.parseInput({varargin,'filterDFOF','kalmanTime'});
            dsTimeFactor = carpaUtilities.parseInput({varargin,'dsTimeFactor',4});
            
            movieReg = double(carpaUtilities.registerMovie(movie));
            cropRegMovie = double(carpaUtilities.cropRegisteredMovie(movieReg));
            
            movieRegFilt = double(carpaUtilities.filterMovie(cropRegMovie,filterReg,freqLow,freqHigh));
            movieDFOF = carpaUtilities.dfofMovie(movieRegFilt);
%             movieDFOFFilt = carpaUtilities.filterMovie(movieDFOF,filterDFOF);
            movieDS = carpaUtilities.downsampleTime(movieDFOF,dsTimeFactor);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%PROCESSING FUNCTIONS%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function cropRegMovie = cropRegisteredMovie(movie)
           %For each pixel, look at wheather in time it has been 0 at least once
            cropMask = ~all(movie,3);
            
            %Ensure we are not removing regions unconected to the borders
            rProps = regionprops(cropMask,'PixelList');
            for k = 1:length(rProps)
                pixelList = rProps(k).PixelList;
                if isempty(intersect(pixelList,[1,size(cropMask,1),size(cropMask,2)]))
                    cropMask(pixelList(:,1),pixelList(:,2)) = 0;
                end
            end
            
            %If there is still something to crop, crop it
            if any(cropMask(:))
                %Find largest rectangle
                bounds = FindLargestRectangles(~cropMask);
                cropRegMovie = movie(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),:);
            else
                cropRegMovie = movie;
            end
        end
        
        function movieREG = registerMovie(movie)
            disp('Registring...')
            movieREG = register1p(movie); %Eliminate first 10 samples?
        end
            
        function movieF = filterMovie(movie, filter, freqLow, freqHigh)
            disp('Filtering...')
            switch filter
                case 'divideByLowpass'
                    movieF = normalizeMovie(single(movie),'normalizationType','lowpassFFTDivisive','freqLow',freqLow,'freqHigh',freqHigh,'waitbarOn',0,'bandpassMask','gaussian');
                case 'bandpass'
                    movieF = normalizeMovie(single(movie),'normalizationType','fft','freqLow',freqLow,'freqHigh',freqHigh,'bandpassType','bandpass','showImages',0,'bandpassMask','gaussian');
                case 'myHighpass'
                    hLarge = fspecial('average', 40);
                    hSmall = fspecial('average', 2);
                    movieF = zeros(size(movie));
                    parfor t = 1:size(movie,3)
                        movieF(:,:,t) = filter2(hSmall,movie(:,:,t)) - filter2(hLarge, movie(:,:,t));
                    end
                case 'kalmanTime'
                    movieF = Kalman_Stack_Filter(movie,.8);
                case 'wiener'
                    movieF = zeros(size(movie));
%                     H = fspecial('log');
                    parfor k = 1:size(movie,3)
                        movieF(:,:,k) = wiener2(movie(:,:,k),[5 5]);%imfilter(movie(:,:,k),H,'replicate');
                    end
                otherwise
                    warning(['Filter ',filter, 'not implemented'])
            end
        end

        function movieDF = dfofMovie(movie)
            disp('DFOFing...')
            movieDF = movie./mean(movie,3) - 1;
        end
        
        function movieDS = downsampleTime(movie,factor)
            disp('Downsampling in time...')
            downX = size(movie,1);
            downY = size(movie,2);
            downZ = floor(size(movie,3)/factor);
            for frame=1:downY
               downsampledFrame = imresize(squeeze(movie(:,frame,:)),[downX downZ],'bilinear');		  
               movie(1:downX,frame,1:downZ) = downsampledFrame;
            end
            movie(:,:,(downZ+1):end) = 0;
            thisMovieTmp = movie(:,:,1:downZ);
            movieDS = thisMovieTmp;        
        end
        
        function movieDS = downsampleSpace(movie,factor)
            disp('Downsampling in space...')
            downX = floor(size(movie,1)/factor);
            downY = floor(size(movie,2)/factor);
            downZ = size(movie,3);
            for frame=1:downZ
               downsampledFrame = imresize(squeeze(movie(:,:,frame)),[downX downY],'bilinear');
               movie(1:downX,1:downY,frame) = downsampledFrame;
            end
            movieDS = movie(1:downX,1:downY,:);
        end
        

        
                %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%UTILITY FUNCTIONS%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        function allFiles = getFilesInFolder(path)
            allItems = dir(path);
            allItems = allItems(3:end); %remove WINDOWS metadata folders
            allFiles = allItems(~[allItems.isdir]);
        end


        
        function fullMovie = readHdf5(movie,varargin)
            keyboard;
            % search_for_dataset
            if isstr(movie)
                data_path = carpaUtilities.getDatasetPath(movie);
                fullMovie = h5read(movie, data_path);
            end

            if iscell(movie)
                if isstr(movie{1}) 
                    %Check out if dimensions are compatible
                    movieSizes = zeros([length(movie),3]);
                    dataPathList = cell(size(movie));
                    for k = 1:length(movie)
                        dataPathList{k} = carpaUtilities.getDatasetPath(movie{k});
                        datasetInfo = h5info(movie{k}, dataPathList{k});  % Get the info of the dataset
                        movieSizes(k,:) = datasetInfo.Dataspace.Size;  % Get the size of the dataset
                    end

                    if length(unique(movieSizes(:,1))) > 1 || length(unique(movieSizes(:,2))) > 1  %same size movies
                        warning('MOVIE DIMENSIONS ARE NOT COMPATIBLE')
                        %Quick and dirty fix
                        dimX = min(movieSizes(:,1));
                        dimY = min(movieSizes(:,2));
                        offsetX = [floor((movieSizes(:,1)-dimX)/2),floor((movieSizes(:,1)-dimX)/2)+mod((movieSizes(:,1)-dimX),2)];
                        offsetY = [floor((movieSizes(:,2)-dimY)/2),floor((movieSizes(:,2)-dimY)/2)+mod((movieSizes(:,2)-dimY),2)];

                        fullMovie = zeros(dimX,dimY, sum(movieSizes(:,3)), 'uint16');
                        for k = 1:length(movie)
                            hinf = hdf5info(movie{k});
                            
                            % sliceMovie = hdf5read(hinf.GroupHierarchy.Datasets);
                            sliceMovie = h5read(movie{k}, dataset_path);

                            sliceMovie = sliceMovie((1+offsetX(k,1)):(end-offsetX(k,2)),(1+offsetY(k,1)):(end-offsetY(k,2)),:);
                            
                            % Calculate the starting and ending index for the 3rd dimension of fullMovie
                            start_frame = sum(movieSizes(1:k-1, 3)) + 1;  % Start frame index based on previous movie sizes
                            end_frame = sum(movieSizes(1:k, 3));  % End frame index for the current movie slice
                        
                            % Assign the sliceMovie to the corresponding part of fullMovie
                            fullMovie(:, :, start_frame:end_frame) = sliceMovie;  % Place sliceMovie in the right spot
                            % fullMovie(:,:,idx_start:idx_end) = sliceMovie;

                        end
                    else  %same size movies
                        fullMovie = zeros(movieSizes(1,1), movieSizes(1,2), sum(movieSizes(:,3)), 'uint16');
                        for  k = 1:length(movie)
                            % original
                            sliceMovie = h5read(movie{k}, dataPathList{k});
                            start_frame = sum(movieSizes(1:k-1, 3)) + 1;  % Start frame index
                            end_frame = sum(movieSizes(1:k, 3));  % End frame index
                            fullMovie(:, :, start_frame:end_frame) = sliceMovie;

                            %  % Process in chunks (e.g., 100 frames at a time)
                            % chunkSize = 100;  % added Set chunk size based on available memory
                            % for chunkStart = 1:chunkSize:movieSizes(k, 3)
                            %     sliceMovie = h5read(movie{k}, dataPathList{k});
                            %     chunkEnd = min(chunkStart + chunkSize - 1, movieSizes(k, 3));  % Handle last chunk
                            % 
                            %     % Extract the chunk of frames
                            %     frameChunk = sliceMovie(:, :, chunkStart:chunkEnd);
                            % 
                            %     % Calculate start and end frame indices for the fullMovie
                            %     start_frame = sum(movieSizes(1:k-1, 3)) + chunkStart;
                            %     end_frame = sum(movieSizes(1:k-1, 3)) + chunkEnd;
                            % 
                            %     % Place the chunk of frames into fullMovie
                            %     fullMovie(:, :, start_frame:end_frame) = frameChunk;
                            % end
                            % 
                            % clear sliceMovie;  % Free memory after processing the current movie part

                       end
                    end
                else
                    fullMovie = [];
                    for moviePart = movie
                        fullMovie = cat(3,fullMovie,moviePart{1});
                    end
                end
            end
        end
       
        function simplifiedName = simplifyName(name)
            words = regexp(name,'([^\W_]*)','tokens');
            simplifiedName = '';
            for word = words
                simplifiedName = [simplifiedName,word{1}{1}];
            end
        end
        
        function match = checkRegexp(sourceText,regexpList)
            match = zeros(size(sourceText));
            k = 1;
            for text = sourceText
                regexpMatch = 0;
                for regexpSeq = regexpList
                    if regexp(text{1},regexpSeq{1})
                        regexpMatch = regexpMatch + 1;
                    end
                end
                if regexpMatch > 0
                    match(k) = 1;
                end
                k = k + 1;
            end
        end
        
        function playMovie(movie,varargin)
            
            deltaT = carpaUtilities.parseInput({varargin,'deltaT',1/20});
            movie = carpaUtilities.readHdf5(movie);
            
            upLim = quantile(movie(:),0.995);
            dwLim = quantile(movie(:),0.005);
            
            figure;
            for k = 1:size(movie,3)
                imagesc(movie(:,:,k),[dwLim,upLim]);
                pause(deltaT)
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function getTracesEventsTraj(extractionFile,decisions,calcium,missingFrames,behavior,experiment,backgroundTreshold)
            if isstr(extractionFile)
                extractionFileF = load(extractionFile);
                dummyField = fields(extractionFileF);
                extractionFile = extractionFileF.(dummyField{1});
            end
            if isstr(decisions)
                decisionsF = load(decisions);
                dummyField = fields(decisionsF);
                decisions = decisionsF.(dummyField{1});
            end
            options = {'backgroundTreshold',backgroundTreshold};
            tracesEvents = getTracesAndEvents(extractionFile, decisions, calcium, missingFrames, behavior, options);
            if ~isempty(tracesEvents)
                [pathstr,name,~] = fileparts(calcium);
                saveNamePart = regexp(name,'(.*\d+)','tokens');
                saveName = [pathstr,filesep,char(saveNamePart{1}),'-',experiment,num2str(backgroundTreshold),'bg','-TracesAndEvents.mat'];
                save(saveName,'tracesEvents','-v7.3')
            end
        end
        
        function pcaicaExtraction(file,experiment,varargin)
            
            [pathstr,name,~] = fileparts(file);
            
            nPC = carpaUtilities.parseInput({varargin,'nPC',750});
            nIC = carpaUtilities.parseInput({varargin,'nIC',500});
            inputDatasetName = carpaUtilities.parseInput({varargin,'inputDatasetName','/1'});
            pcaicaOutputUnits = carpaUtilities.parseInput({varargin,'pcaicaOutputUnits','fl'});
            
            [PcaOutputSpatial, PcaOutputTemporal, PcaOutputSingularValues, PcaInfo] = run_pca({file}, nPC, 'movie_dataset_name',inputDatasetName);
            
            calcInf = hdf5info(file);
            movieDims.x = calcInf.GroupHierarchy.Datasets.Dims(1);
            movieDims.y = calcInf.GroupHierarchy.Datasets.Dims(2);
            
            [IcaFilters, IcaTraces, IcaInfo] = run_ica(PcaOutputSpatial, PcaOutputTemporal, PcaOutputSingularValues, movieDims.x, movieDims.y, nIC, 'output_units',pcaicaOutputUnits);
            IcaTraces = permute(IcaTraces,[2 1]);
            
            pcaicaAnalysisOutput.IcaInfo = IcaInfo;
            pcaicaAnalysisOutput.PcaInfo = PcaInfo;
            pcaicaAnalysisOutput.IcaFilters = IcaFilters;
            pcaicaAnalysisOutput.IcaTraces = IcaTraces;
            pcaicaAnalysisOutput.nPCs = nPC;
            pcaicaAnalysisOutput.nICs = nIC;
            pcaicaAnalysisOutput.movieFilename = file;
            
            saveNamePart = regexp(name,'(.*\d+)','tokens');
            saveName = [pathstr,filesep,char(saveNamePart{1}),'-',experiment,'-pcaicaAnalysis.mat'];
            save(saveName,'pcaicaAnalysisOutput','-v7.3');
            
        end
        
        function cellMaxExtraction(file,experiment,varargin)
            
            [pathstr,name,~] = fileparts(file);
            
            emOptions.useParallel = carpaUtilities.parseInput({varargin,'useParallel',1});
            emOptions.movieDatasetName = carpaUtilities.parseInput({varargin,'movieDatasetName','/1'});
            
            CELLMaxoptions.initMethod = 'grid';
            CELLMaxoptions.gridSpacing = 9;
            CELLMaxoptions.gridWidth = 4;
            CELLMaxoptions.inputSizeManual = 0;
            CELLMaxoptions.subsampleMethod = 'resampleRemaining';
            CELLMaxoptions.percentFramesPerIteration = 0.5000;
            CELLMaxoptions.percentRemainingSubsample = 0.7500;
            CELLMaxoptions.maxSqSize = 101;
            CELLMaxoptions.threshForElim = 0.0050;
            CELLMaxoptions.localICimgs = [];
            CELLMaxoptions.localICtraces = [];
            CELLMaxoptions.minIters = 200;
            CELLMaxoptions.maxIters = 460;
            CELLMaxoptions.numSigmasThresh = 0.5000;
            CELLMaxoptions.nParallelWorkers = inf;
            CELLMaxoptions.generateNovelSeed = 1;
            CELLMaxoptions.numFramesRandom = 2000;
            CELLMaxoptions.readMovieChunks = 1;
            CELLMaxoptions.movieFilename = file;
            emOptions.CELLMaxoptions = carpaUtilities.parseInput({varargin,'CELLMaxoptions',CELLMaxoptions});
            
            startTime = tic;
            [emAnalysisOutput, ~] = CELLMax_Wrapper(file,'options',emOptions);
            
            emOptions.CELLMaxoptions.sqSizeX = [];
            emOptions.CELLMaxoptions.sqSizeY = [];

            emAnalysisOutput.dsCellTraces = emAnalysisOutput.cellTraces;
            emOptions.CELLMaxoptions.numSignalsDetected = size(emAnalysisOutput.dsCellTraces,1);
            emOptions.versionCellmax = emAnalysisOutput.versionCellmax;
            emOptions.time.startTime = startTime;
            emOptions.time.endTime = toc(startTime);
            
            saveNamePart = regexp(name,'(.*\d+)','tokens');
            saveName = [pathstr,filesep,char(saveNamePart{1}),'-',experiment,'-emAnalysis.mat'];
            save(saveName,'emAnalysisOutput','-v7.3','emOptions');
            rmdir(fullfile(pathstr,'tmpImages'),'s')
        end
        
        function decideCells(moviePath,analysisOutput,type)
            cInf = cellInfo(moviePath,permute(analysisOutput.cellImages,[3,1,2]),analysisOutput.scaledProbability,[]);
            cInf.setTresholds('showProgress',0,'tresholds',{'getOverlap','>0.5','getglobalSNR','>1','getCellShapePropieties.Eccentricity','<=0.9','getCellShapePropieties.Area','>100','getsScore','<0.02'});
            cInf.saveDecisions(type,'skipConfirmation',1)
        end
        
        function createHDF5(filePath,movie,varargin)
            datasetName = carpaUtilities.parseInput({varargin,'datasetName','/1'});
            try
                createHdf5File(filePath,datasetName,movie);
            catch
                disp('filename may be too large or directory may not exist, fix manually')
                pause;
            end
        end
        
        function idx = compareNamesTillExperiment(name,names,experiment)
            if iscell(name)
                error('first argument should be a str')
            end
            
            if ~iscell(names)
                names = {names};
            end
            
            try
                nameCut = [char(extractBefore(name,experiment)),experiment];
            catch
                error(['REVISE FILENAME OF ',name, ' AND MAKE SURE IT MATCHES THE EXPERIMENT ', experiment])
            end
            namesCut = cell([1,length(names)]);
            for k = 1:length(names)
                try
                    namesCut{k} = [char(extractBefore(names{k},experiment)),experiment];
                end
            end
            
            idx = find(strcmp(nameCut,namesCut));
            
        end

        function createLogFile(obj,dateTime)
            % This function creates a new log file in the specified folder with the specified XML format.
            
            % Extract information from obj
            logFolder = obj.folderStruct.path;  % Assuming 'folderStruct' is a field in obj
            dateStr = datestr(now, 'yyyyMMdd');  % Current date
            mouseCode = obj.mouse;
            mouseDate = obj.dateParser;
            mouseTime = obj.timeParser;
            mouseExperiment = obj.folderStruct.experiment;
        
            % Define the log file name using the extracted information
            logFileName = fullfile(logFolder, sprintf('Mouse-%s-%s-%s-log.xml', mouseCode,dateTime,mouseExperiment));
        
            % Open a new file to write the log
            fid = fopen(logFileName, 'w');
        
            if fid == -1
                error('Unable to create log file.');
            end
             
            filePath = strcat(obj.folderStruct.path, '\', obj.folderStruct.concat.fileName);
            info = h5info(filePath);  % Get HDF5 file info
            datasetSize = info.Datasets.Dataspace.Size;
            numFrames = datasetSize(3);



        
            % Write the XML structure to the file
            fprintf(fid, ['<recording>\n' ...
                          '<attrs>\n' ...
                          '<attr name="version">2.0.32-20150915-065415</attr>\n' ...
                          '<attr name="width">1080</attr>\n' ...
                          '<attr name="height">1080</attr>\n' ...
                          '<attr name="left">360</attr>\n' ...
                          '<attr name="top">0</attr>\n' ...
                          '<attr name="fps">20.00</attr>\n' ...
                          '<attr name="exposure">49.664</attr>\n' ...
                          '<attr name="time">[]</attr>\n' ...
                          '<attr name="frames">%d</attr>\n' ...
                          '<attr name="dropped_count">0</attr>\n' ...
                          '<attr name="dropped">[]</attr>\n' ...
                          '<attr name="downsample">1X</attr>\n' ...
                          '<attr name="gain">[]</attr>\n' ...
                          '<attr name="led_power">[]</attr>\n' ...
                          '<attr name="led_delay_value">[]</attr>\n' ...
                          '<attr name="led_session">[]</attr>\n' ...
                          '<attr name="led_project">[]</attr>\n' ...
                          '<attr name="record_start">[]</attr>\n' ...
                          '<attr name="record_end">[]</attr>\n' ...
                          '<attr name="record_triggered">False</attr>\n' ...
                          '<attr name="record_sched_batch">[]</attr>\n' ...
                          '<attr name="record_sched_name">[]</attr>\n' ...
                          '<attr name="record_sched_step">[]</attr>\n' ...
                          '<attr name="record_sched_cycle">[]</attr>\n' ...
                          '<attr name="rois">[]</attr>\n' ...
                          '<attr name="meta_data">[]</attr>\n' ...
                          '<attr name="camera_chip_version">0x1402L</attr>\n' ...
                          '<attr name="sensor_board_serial_number">A-201020632</attr>\n' ...
                          '<attr name="hardware_serial_number">AA-1211150428</attr>\n' ...
                          '</attrs>\n' ...
                          '<decompressed>\n' ...
                          '<file frames="%d">%s</file>\n' ...
                          '</decompressed>\n' ...
                          '</recording>\n'],numFrames,numFrames, obj.folderStruct.concat.fileName);
        
            % Close the file
            fclose(fid);
        
            % Display a message that the file was created
            disp(['Log file created: ', logFileName]);
        end


        
        function maskedFiles = maskSelectedFiles(selectedFiles,maskFiles)
            allSelectedFiles = cat(2,selectedFiles{:});
            allMaskFiles = cat(2,maskFiles{:});
            if ~isempty(allMaskFiles)
                allMasks = cat(1,allMaskFiles.date);
            else
                allMasks = {};
            end
            %group by days
            maskedFiles = {};
            folderList = {};
            for k = 1:length(allSelectedFiles)
                if isempty(allMasks) || ~any(ismember(allSelectedFiles(k).date,allMasks,'rows'))
                    if ~any(strcmp(folderList,allSelectedFiles(k).folder))
                        folderList = [folderList, allSelectedFiles(k).folder];
                        maskedFiles = [maskedFiles,{allSelectedFiles(k)}];
                    else
                        maskedFiles{end} = [maskedFiles{end},allSelectedFiles(k)];
                    end
                end
            end 
        end
        
        function new_files = delete_list(files,folder,this_regexp)
            file_mask = {};
            for k_reg = 1:length(this_regexp)
                newstr = strrep(this_regexp{k_reg},'$','');
                splitstr = split(newstr,'*');
                file_mask = [file_mask,['*',splitstr{end}]];
            end
            
            to_delete = uigetfile(fullfile(folder,char(join(file_mask,';'))),'Select the files that will be deleted.','MultiSelect','on');
            if ~iscell(to_delete)
                to_delete = {to_delete};
            end
            new_files = setdiff(files,cellfun(@(x) fullfile(folder,x), to_delete,'UniformOutput',0));
        end
        
        function new_files = merge_list(files,folder,this_regexp)
            file_mask = {};
            for k_reg = 1:length(this_regexp)
                newstr = strrep(this_regexp{k_reg},'$','');
                splitstr = split(newstr,'*');
                file_mask = [file_mask,['*',splitstr{end}]];
            end

            not_done_merging = true;
            allMerged = {};
            while not_done_merging
                merged = uigetfile(fullfile(folder,char(join(file_mask,';'))),'Select the files that will be merged. Press cancel when done','MultiSelect','on');
                if ~iscell(merged)
                    if merged == 0
                        not_done_merging = false;
                    else
                        error('Need to merge at least 2 files')
                    end
                else
                    allMerged =[allMerged,{merged}];
                end
            end

            all_idx = {};
            for k_merg = 1:length(allMerged)
                merged = allMerged{k_merg};
                idx = cellfun(@(x) find(~cellfun(@isempty,(strfind(files,fullfile(folder,x))))), merged);
                if any(diff(idx) > 1)
                    error('Cannot merge non continuous files')
                end
                all_idx = [all_idx,idx];
            end
            if length(unique([all_idx{:}])) ~= sum(cellfun(@length,all_idx))
                error('Overlaping merge groups')
            end

            new_files = {};
            idx = all_idx{1};
            if idx(1)~=1
                idx_before_merg = 1:(idx(1)-1);
                new_files = [new_files,files(idx_before_merg)];
            end
            for k_merg = 1:length(all_idx)
                idx = all_idx{k_merg};
                new_files = [new_files,{files(idx)}];
                if k_merg ~= length(all_idx)
                    new_files = [new_files,files((idx(end)+1):(all_idx{k_merg+1}-1))];
                end
            end
            idx = all_idx{end};
            if idx(end)~=length(files)
                idx_before_merg = (idx(end)+1):length(files);
                new_files = [new_files,files(idx_before_merg)];
            end
        end
    
        function [position,velocity] = compute_trajectory(file,varargin)
            [position,velocity] = getMouseTrajectory(file,varargin{1}{:});
%             try
%                 [position,velocity,score] = getMouseTrajectory(file);
%             catch
%                 score = 1;
%             end
%             if score > .1
%                 disp('First tracker with grave errors, trying python...')
%                 %disp(behavior{block})
%                 [positionPy,velocityPy,scorePy] = pyTrackerInterface(file);
%                 if scorePy > .1
%                     disp('First and second trackers with grave errors, taking best one...')
%                 end
%                 if scorePy <= score
%                     position = positionPy;
%                     velocity = velocityPy;
%                 end
%             end
        end
        function datasetPath = getDatasetPath(filePath, varargin)
            verbose = carpaUtilities.parseInput({varargin,'verbose',false});
            % This function reads the information of an HDF5 file and finds a dataset with 3 dimensions
            info = h5info(filePath);  % Get HDF5 file info

            % Initialize variable to store the path of the dataset with 3 dimensions
            currentPath = '';
            % Explore each group in the HDF5 file
            for i = 1:length(info.Groups)
                % Try to find the dataset with exactly 3 dimensions
                datasetPath = carpaUtilities.exploreGroup(info.Groups(i), filePath, currentPath, verbose);
                if ~isempty(datasetPath)
                    if verbose; disp(['Found dataset path: ', datasetPath]); end
                    return;  % If dataset is found, exit early
                end
            end
            
            % Check for datasets at the root level
            if isfield(info, 'Datasets') && ~isempty(info.Datasets)
                for j = 1:length(info.Datasets)
                    currentPath = '/';
                    datasetPath = carpaUtilities.checkDatasetDimensions(info.Datasets(j), currentPath, verbose);
                    if ~isempty(datasetPath)
                        if verbose; disp(['Found dataset path: ', datasetPath]); end
                        return;  % Exit early if found
                    end
                end
            end

            % After exploring all groups and root-level datasets, return the path of the dataset
            if isempty(datasetPath)
                if verbose; disp('No dataset with 3 dimensions found.'); end
            end
        end
                % exploreGroup
        function datasetPath = exploreGroup(group, filePath, currentPath, verbose)
            datasetPath  = '';
            % Display the current group name
            if verbose; disp(['Group: ', group.Name]); end
            % Check if this group contains datasets
            if isfield(group, 'Datasets') && ~isempty(group.Datasets)
                % Iterate over datasets in this group
                for j = 1:length(group.Datasets)
                    datasetPath = carpaUtilities.checkDatasetDimensions(group.Datasets(j), [currentPath, group.Name], verbose);
                    if ~isempty(datasetPath)
                        return;  % Exit early if found
                    end
                end
            else
                % If the group does not contain any datasets
                if verbose; disp('  No datasets found in this group.'); end
            end
            
            % Recursively explore subgroups
            if isfield(group, 'Groups')
                for k = 1:length(group.Groups)
                    datasetPath = carpaUtilities.exploreGroup(group.Groups(k), filePath, currentPath, verbose);
                    if ~isempty(datasetPath)
                        return;  % Exit early if dataset is found
                    end
                end
            end
        end

        % checkDatasetDimensions
        function datasetPath = checkDatasetDimensions(dataset, currentPath, verbose)

            if currentPath(end) ~= '/'
                currentPath = [currentPath, '/'];
            end
            % Display the size of the dataset
            datasetSize = carpaUtilities.getDatasetSize(dataset);
            if verbose; fprintf('  Dataset: %s', dataset.Name); end
            if numel(datasetSize) == 3
                if verbose; fprintf(' | size %s \n',num2str(datasetSize)); end
                % Check if the dataset has exactly 3 dimensions
                % Check the third dimension (frames)
                numFrames = datasetSize(3);
                % If the third dimension is more than 100 frames, proceed automatically
                if numFrames > 100
                    datasetPath = [currentPath, dataset.Name];  % Append "/"
                else
                    % Ask for user confirmation
                    disp(['    Third dimension has ', num2str(numFrames), ' frames (<= 100).']);
                    user_response = input('Would you like to confirm this dataset? (y/n): ', 's');
                    if lower(user_response) == 'y'
                        % If the group is the root ("/"), don't add an extra "/"
                        datasetPath = [groupName, dataset.Name];  % Append "/"
                    else
                        datasetPath = '';
                    end
                end
            elseif isempty(datasetSize)
                if verbose; fprintf(' | size not available\n'); end
                datasetPath = '';
            else
                if verbose; fprintf(' | non 3D size found\n'); end
                datasetPath = '';
            end
        end

        function datasetSize = getDatasetSize(dataset)
            % Try to get the dataset size, return empty if size is not available
            if isfield(dataset, 'Dataspace') && isfield(dataset.Dataspace, 'Size')
                datasetSize = dataset.Dataspace.Size;  
            elseif isfield(dataset, 'Dims')
                datasetSize = dataset.Dims;  
            else
                datasetSize = '';  
            end
        end


    end
    
    methods(Static, Access = private)

        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%UTILITY FUNCTIONS%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        function [varargout] = parseInput(varargin)
            %{{inputs},optional(arg1 default1, arg2 default2, ...)}
            tmpVar = varargin{1};
            if iscell(tmpVar{1})
                argumentList = tmpVar(2:2:end);
                defaultVal = tmpVar(3:2:end);
                input = tmpVar{1};
                inputArg = input(1:2:end);
                inputVal = input(2:2:end);
                for k = 1:length(argumentList)
                    idx = strcmp(argumentList{k},inputArg);
                    if sum(idx) == 0
                        varargout{k} = defaultVal{k};
                    else
                        varargout{k} = inputVal{idx};
                    end
                end                
            else
                for k = 2:2:length(tmpVar)
                    varargout{k} = tmpVar{k};
                end
            end
        end
    end
        
end
