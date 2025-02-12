function obj = modelVerifyDataIntegrity(obj)
    % get information for each folder
    % biafra ahanonu
    % branched from controllerAnalysis: 2014.08.01 [16:09:16]
    % inputs
        %
    % outputs
        %

    % changelog
        %
    % TODO
        %

    scnsize = get(0,'ScreenSize');
    signalExtractionMethodStr = {'hdf5_datasetnames','movies_signals','duplicates','movieStatistics','stimulusIndex'};
    [fileIdxArray, ok] = listdlg('ListString',signalExtractionMethodStr,'ListSize',[scnsize(3)*0.2 scnsize(4)*0.25],'Name','which signal extraction method?');
    analysisType = signalExtractionMethodStr{fileIdxArray};
    [fileIdxArray idNumIdxArray nFilesToAnalyze nFiles] = obj.getAnalysisSubsetsToAnalyze();
    currentDateTimeStr = datestr(now,'yyyymmdd_HHMM','local');

    switch analysisType
        case 'stimulusIndex'
            for idNumIdx = 1:length(idNumIdxArray)
                idNum = idNumIdxArray(idNumIdx);
                obj.stimNum = idNum;
                try
                    % =====================
                    display(repmat('=',1,7))
                    display([num2str(idNum) '/' num2str(nIDs) ': analyzing ' nameArray{idNum}])

                    % ===============================================================
                    % obtain stimulus information
                    idArray(idNum)
                    stimVector = obj.modelGetStim(idArray(idNum));
                    if isempty(stimVector); continue; end;
                    framesToAlign = find(stimVector);
                    % maxTrials = 100;
                    if length(framesToAlign)>=maxTrials
                        % framesToAlign = framesToAlign(1:20);
                        framesToAlign = framesToAlign(randperm(length(framesToAlign),maxTrials));
                    end
                    nPoints = size(IcaTraces,2);
                    timeVector = [-preOffset:postOffset];
                    framesToAlign(find((framesToAlign<preOffset))) = [];
                    framesToAlign(find((framesToAlign>(nPoints-postOffset)))) = [];
                    [~, ~] = openFigure(776, '');
                        options.videoTrialRegExp = [obj.date{obj.fileNum} '_' obj.protocol{obj.fileNum} '_' obj.fileIDArray{obj.fileNum}];
                        vidList = getFileList(obj.videoDir,options.videoTrialRegExp);
                        [xPlot yPlot] = getSubplotDimensions(length(framesToAlign));
                        downsampleFactor = 4;
                        length(framesToAlign)
                        for trialNo = 1:length(framesToAlign)
                            % subplot(xPlot,yPlot,trialNo)
                            % behaviorMovie2 = loadMovieList(vidList,'convertToDouble',0,'frameList',bsxfun(@plus,framesToAlign(trialNo),0:2)*downsampleFactor,'treatMoviesAsContinuous',1);
                            thisMovie1 = convertInputMovieToCell(loadMovieList(vidList,'convertToDouble',0,'frameList',bsxfun(@plus,framesToAlign(trialNo),-2:2)*downsampleFactor,'treatMoviesAsContinuous',1));
                            thisMovie1 = cat(2,thisMovie1{:});
                            imagesc(thisMovie1);
                            ginput(1);
                        end
                        continue
                catch err
                    display(repmat('@',1,7))
                    disp(getReport(err,'extended','hyperlinks','on'));
                    display(repmat('@',1,7))
                end
            end
        case 'cell_extraction'
            % file = dir(filesToLoad{fileNo});
            % file.date
        case 'movies_signals'
            movieSettings = inputdlg({...
                    'processed imaging movie regexp:',...
                    'raw imaging movie regexp:',...
                    'recording file regexp:',...
                    'recording file regexp ext:'...
                },...
                'Settings for HDF5 data checking',[1 100],...
                {...
                    obj.fileFilterRegexp,...
                    'concat',...
                    'recording.*.(txt|xml)',...
                    '(.txt|.xml)'...
                }...
            );
            fileFilterRegexp = movieSettings{1};
            rawFileFilterRegexp = movieSettings{2};
            recordingFileFilterRegexp = movieSettings{3};
            recordingFileFilterRegexpExt = movieSettings{4};
            originalSignalExtractionMethod = obj.signalExtractionMethod;
            nFolders = length(obj.inputFolders);
            for folderNo = 1:nFolders
                display(repmat('=',1,21))
                display([num2str(folderNo) '/' num2str(nFolders) ': ' obj.inputFolders{folderNo}])
                signalSets = {{obj.rawICfiltersSaveStr,obj.rawICtracesSaveStr,obj.rawPCAICAStructSaveStr,obj.sortedICdecisionsSaveStr},{obj.rawEMStructSaveStr,obj.sortedEMStructSaveStr},{obj.rawROItracesSaveStr},{obj.rawEXTRACTStructSaveStr,obj.sortedEXTRACTStructSaveStr},{obj.rawCNMFStructSaveStr,obj.sortedCNMFStructSaveStr}};
                signalSetStr = {'PCA-ICA','CELLMax','ROI','EXTRACT','CNMF'};
                for signalSetNo = 1:length(signalSets)
                    filesToLoad = getFileList(obj.inputFolders{folderNo},signalSets{signalSetNo});
                    if isempty(filesToLoad)
                        display(['no ' signalSetStr{signalSetNo} ' files: ' obj.inputFolders{folderNo}])
                    else
                        for fileNo = 1:length(filesToLoad)
                            file = dir(filesToLoad{fileNo});
                            display([file.date ' | ' filesToLoad{fileNo}]);
                        end
                    end
                end
                % datenum

                filesToLoad = getFileList(obj.inputFolders{folderNo},fileFilterRegexp);
                if isempty(filesToLoad)
                    display(['no ' fileFilterRegexp ' files: ' obj.inputFolders{folderNo}])
                else
                    for fileNo = 1:length(filesToLoad)
                        file = dir(filesToLoad{fileNo});
                        display([file.date ' | ' filesToLoad{fileNo}]);
                    end
                end
                % if isempty(filesToLoad)
                %     display(['no dfof files: ' obj.inputFolders{folderNo}])
                %     filesToLoad = getFileList(obj.inputFolders{folderNo},rawFileFilterRegexp);
                %     [movieDims] = loadMovieList(filesToLoad,'getMovieDims',1,'inputDatasetName',obj.inputDatasetName);
                %     display(['raw movie length: ' num2str(sum(movieDims.z))]);
                % else
                %     filesToLoad = getFileList(obj.inputFolders{folderNo},'crop');
                %     [movieDims] = loadMovieList(filesToLoad,'getMovieDims',1,'inputDatasetName',obj.inputDatasetName);
                %     display(['processed movie length: ' num2str(movieDims.z)]);
                % end
            end
        case 'hdf5_datasetnames'
            nFolders = length(obj.inputFolders);
            % get user input
                if ~isempty(obj.videoDir)
                    if iscell(obj.videoDir)
                        videoDir = strjoin(obj.videoDir,',');
                    else
                        videoDir = obj.videoDir;
                    end
                else
                    videoDir = '';
                end
            movieSettings = inputdlg({...
                    'processed imaging movie regexp:',...
                    'raw imaging movie regexp:',...
                    'recording file regexp:',...
                    'recording file regexp ext:',...
                    'video file path:'...
                },...
                'Settings for HDF5 data checking',1,...
                {...
                    obj.fileFilterRegexp,...
                    'concat',...
                    'recording.*.(txt|xml)',...
                    '(.txt|.xml)',...
                    strjoin(obj.videoDir,',')...
                }...
            );
            fileFilterRegexp = movieSettings{1};obj.fileFilterRegexp = fileFilterRegexp;
            rawFileFilterRegexp = movieSettings{2};
            recordingFileFilterRegexp = movieSettings{3};
            recordingFileFilterRegexpExt = movieSettings{4};
            videFilePath = movieSettings{5};obj.videoDir = videFilePath;

            videoTrialRegExp = '';
            videoTrialRegExpIdx = 1;
            if ~isempty(videFilePath)
                videoTrialRegExpList = {'yyyy_mm_dd_pNNN_mNNN_assayNN','yymmdd-mNNN-assayNN','yymmdd_mNNN_assayNN','subject_assay','yymmdd_mNNN'};
                scnsize = get(0,'ScreenSize');
                [videoTrialRegExpIdx, ok] = listdlg('ListString',videoTrialRegExpList,'ListSize',[scnsize(3)*0.2 scnsize(4)*0.25],'Name','video string type (N = number)');
                % videoTrialRegExpList = {'yyyy_mm_dd_pNNN_mNNN_assayNN','yymmdd-mNNN-assayNN','subject_assay'};
                % scnsize = get(0,'ScreenSize');
                % [videoTrialRegExpIdx, ok] = listdlg('ListString',videoTrialRegExpList,'ListSize',[scnsize(3)*0.2 scnsize(4)*0.25],'Name','video string type (N = number)');
            else

            end

            movieIntegrityTable = table(...
                {'0000.00.00'},...
                {'tmp'},...
                {'tmp'},...
                0,...
                0,...
                0,...
                0,...
                0,...
                0,...
                'VariableNames',{...
                'date',...
                'folder',...
                'folderName',...
                'movieProcessedFramesCount',...
                'movieRawFramesCount',...
                'movieRawFramesCountDownsampled',...
                'textRawFramesCount',...
                'textRawFramesCountDownsampled',...
                'textDroppedFramesCount'})

            movieBehaviorIntegrityTable = table(...
                {'0000.00.00'},...
                {'tmp'},...
                {'tmp'},...
                {'tmp'},...
                0,...
                0,...
                0,...
                0,...
                0,...
                0,...
                0,...
                'VariableNames',{...
                'date',...
                'folder',...
                'folderName',...
                'trialName',...
                'movieRawFramesCount',...
                'textRawFramesCount',...
                'textDroppedFramesCount',...
                'camera01FrameCount',...
                'camera02FrameCount',...
                'camera01FrameCountDiff',...
                'camera02FrameCountDiff'});

            for folderNo = fileIdxArray
                addRowTwo = size(movieBehaviorIntegrityTable,1)+1;
                addRow = size(movieIntegrityTable,1)+1;

                % display('==========')
                display(repmat('=',1,21))
                display([num2str(folderNo) '/' num2str(nFolders) ': ' obj.inputFolders{folderNo}])
                obj.fileNum = folderNo;
                local_getVideoRegexp();
                % ==============================
                processsedFilesToLoad = getFileList(obj.inputFolders{folderNo},fileFilterRegexp);
                rawFilesToLoad = getFileList(obj.inputFolders{folderNo},rawFileFilterRegexp);
                filesToLoad = [processsedFilesToLoad(:);rawFilesToLoad(:)];
                totalProcessedFramesCount = 0;
                totalRawFramesCount = 0;
                addNum = 1;
                for fileNo = 1:length(filesToLoad)
                    if exist(filesToLoad{fileNo},'dir')~=0
                        continue
                    end
                    [pathstr,name,ext] = fileparts(filesToLoad{fileNo});
                    if strcmp(ext,'.h5')|strcmp(ext,'.hdf5')
                        fileinfo = hdf5info(filesToLoad{fileNo});
                        [movieDims] = loadMovieList(filesToLoad{fileNo},'getMovieDims',1,'inputDatasetName',fileinfo.GroupHierarchy.Datasets.Name);
                    else
                        [movieDims] = loadMovieList(filesToLoad{fileNo},'getMovieDims',1);
                    end
                    % num2str(sum(movieDims.z))
                    file = dir(filesToLoad{fileNo});

                    display([file.date ' | ' fileinfo.GroupHierarchy.Datasets.Name ' | ' num2str([movieDims.x movieDims.y movieDims.z]) ' | ' filesToLoad{fileNo}]);

                    if sum(strcmp(rawFilesToLoad,filesToLoad{fileNo}))>0
                        totalRawFramesCount = totalRawFramesCount + movieDims.z;

                        % add to behavior only movie
                        addRowTwoTmp = addRowTwo-1+addNum;
                        warning off;
                        movieBehaviorIntegrityTable.date{addRowTwoTmp,1} = datestr(now,'yyyy.mm.dd','local');
                        movieBehaviorIntegrityTable.folder{addRowTwoTmp,1} = obj.inputFolders{folderNo};
                        movieBehaviorIntegrityTable.folderName{addRowTwoTmp,1} = obj.folderBaseSaveStr{folderNo};
                        movieBehaviorIntegrityTable.movieRawFramesCount(addRowTwoTmp,1) = movieDims.z;
                        warning on;
                        addNum = addNum + 1;
                    end
                    if sum(strcmp(processsedFilesToLoad,filesToLoad{fileNo}))>0
                        totalProcessedFramesCount = totalProcessedFramesCount + movieDims.z;
                    end

                end
                display(['total processed frames: ' num2str(totalProcessedFramesCount)])
                display(['total raw frames: ' num2str(totalRawFramesCount)  '(' num2str(floor(totalRawFramesCount/obj.DOWNSAMPLE_FACTOR)) ')'])
                % add results to table
                movieIntegrityTable.date{addRow,1} = datestr(now,'yyyy.mm.dd','local');
                movieIntegrityTable.folder{addRow,1} = obj.inputFolders{folderNo};
                movieIntegrityTable.folderName{addRow,1} = obj.folderBaseSaveStr{folderNo};
                movieIntegrityTable.movieProcessedFramesCount(addRow,1) = totalProcessedFramesCount;
                movieIntegrityTable.movieRawFramesCount(addRow,1) = totalRawFramesCount;
                movieIntegrityTable.movieRawFramesCountDownsampled(addRow,1) = floor(totalRawFramesCount/obj.DOWNSAMPLE_FACTOR);
                % ==============================
                % get information about number of frames and dropped frames from recording files
                filesToLoad = getFileList(obj.inputFolders{folderNo},recordingFileFilterRegexp);
                totalRawFramesCount = 0;
                totalDropCount = 0;
                pastBaseTxtFileName = '';
                for fileNo = 1:length(filesToLoad)
                    % check whether have dubplicate xml and txt, skip if so
                    currentBaseTxtFileName = regexprep(filesToLoad{fileNo},recordingFileFilterRegexpExt,'');
                    if strcmp(pastBaseTxtFileName,currentBaseTxtFileName)==1
                        continue
                    end
                    pastBaseTxtFileName = currentBaseTxtFileName;

                    try
                        logInfo = getLogInfo(filesToLoad{fileNo});
                        if isfield(logInfo,'null')
                            continue;
                        end
                    catch err
                        continue
                        display(repmat('@',1,7))
                        disp(getReport(err,'extended','hyperlinks','on'));
                        display(repmat('@',1,7))
                    end

                    fileType = logInfo.fileType;
                    switch fileType
                        case 'inscopix'
                            movieFrames = logInfo.FRAMES;
                            droppedCount = logInfo.DROPPED_COUNT;
                        case 'inscopixXML'
                            movieFrames = logInfo.frames;
                            droppedCount = logInfo.dropped_count;
                        otherwise
                            % do nothing
                    end
                    if ischar(droppedCount)
                        droppedCount = str2num(droppedCount);
                    end
                    if ischar(movieFrames)
                        movieFrames = str2num(movieFrames);
                    end
                    file = dir(filesToLoad{fileNo});

                    display([file.date ' | ' num2str([movieFrames droppedCount]) ' | ' filesToLoad{fileNo}]);

                    totalRawFramesCount = totalRawFramesCount + movieFrames + droppedCount;
                    totalDropCount = totalDropCount+droppedCount;

                    addRowTwoTmp = addRowTwo-1+fileNo;
                    movieBehaviorIntegrityTable.textRawFramesCount(addRowTwoTmp,1) = movieFrames + droppedCount;
                    movieBehaviorIntegrityTable.textDroppedFramesCount(addRowTwoTmp,1) = droppedCount;
                end
                display(['total raw frames: ' num2str(totalRawFramesCount) '(' num2str(floor(totalRawFramesCount/obj.DOWNSAMPLE_FACTOR)) ')'])
                % ==============================
                if ~isempty(videFilePath)
                    behaviorVideosFilesToLoad = getFileList(videFilePath,videoTrialRegExp);
                    numCameras = round(length(behaviorVideosFilesToLoad)/length(rawFilesToLoad));
                    addNum = 1;
                    for fileNo = 1:numCameras:length(behaviorVideosFilesToLoad)
                        if exist(behaviorVideosFilesToLoad{fileNo},'dir')~=0
                            continue
                        end
                        if numCameras==1
                            [movieDims] = loadMovieList(behaviorVideosFilesToLoad{fileNo},'getMovieDims',1,'inputDatasetName',fileinfo.GroupHierarchy.Datasets.Name);
                            camera01FrameCount = movieDims.z;
                            camera02FrameCount = NaN;
                            file = dir(behaviorVideosFilesToLoad{fileNo});
                            display([file.date ' | ' num2str([movieDims.x movieDims.y movieDims.z]) ' | ' behaviorVideosFilesToLoad{fileNo}]);
                        elseif numCameras>1
                            [movieDims] = loadMovieList(behaviorVideosFilesToLoad{fileNo},'getMovieDims',1,'inputDatasetName',fileinfo.GroupHierarchy.Datasets.Name);
                            [movieDims2] = loadMovieList(behaviorVideosFilesToLoad{fileNo+1},'getMovieDims',1,'inputDatasetName',fileinfo.GroupHierarchy.Datasets.Name);
                            camera01FrameCount = movieDims.z;
                            camera02FrameCount = movieDims2.z;

                            file = dir(behaviorVideosFilesToLoad{fileNo});
                            display([file.date ' | ' num2str([movieDims.x movieDims.y movieDims.z]) ' | ' behaviorVideosFilesToLoad{fileNo}]);
                            file = dir(behaviorVideosFilesToLoad{fileNo+1});
                            display([file.date ' | ' num2str([movieDims.x movieDims.y movieDims.z]) ' | ' behaviorVideosFilesToLoad{fileNo+1}]);
                        end

                        thisFileInfo = getFileInfo(behaviorVideosFilesToLoad{fileNo});
                        thisTrial = thisFileInfo.trial;

                        addRowTwoTmp = addRowTwo-1+addNum;
                        movieBehaviorIntegrityTable.trialName{addRowTwoTmp,1} = thisTrial;
                        movieBehaviorIntegrityTable.camera01FrameCount(addRowTwoTmp,1) = camera01FrameCount;
                        movieBehaviorIntegrityTable.camera02FrameCount(addRowTwoTmp,1) = camera02FrameCount;
                        movieBehaviorIntegrityTable.camera01FrameCountDiff(addRowTwoTmp,1) = movieBehaviorIntegrityTable.textRawFramesCount(addRowTwoTmp,1) - camera01FrameCount;
                        movieBehaviorIntegrityTable.camera02FrameCountDiff(addRowTwoTmp,1) = movieBehaviorIntegrityTable.textRawFramesCount(addRowTwoTmp,1) - camera02FrameCount;
                        addNum = addNum+1;
                    end

                end
                % ==============================
                % add results to table
                movieIntegrityTable.textRawFramesCount(addRow,1) = totalRawFramesCount;
                movieIntegrityTable.textRawFramesCountDownsampled(addRow,1) = floor(totalRawFramesCount/obj.DOWNSAMPLE_FACTOR);
                movieIntegrityTable.textDroppedFramesCount(addRow,1) = totalDropCount;
                % ==============================

            end

            movieIntegrityTable = movieIntegrityTable(2:end,:);
            movieIntegrityTable.processedTextEqualFrames = movieIntegrityTable.movieProcessedFramesCount==movieIntegrityTable.textRawFramesCountDownsampled;
            movieIntegrityTable.processedTextDiffFrames = movieIntegrityTable.movieProcessedFramesCount - movieIntegrityTable.textRawFramesCountDownsampled;
            movieIntegrityTable.rawTextMovieDiffFrames = (movieIntegrityTable.textRawFramesCount - movieIntegrityTable.textDroppedFramesCount) - movieIntegrityTable.movieRawFramesCount;

            movieIntegrityTablePath = [obj.dataSavePath filesep obj.protocol{1} '_database_movieIntegrityCheck_' currentDateTimeStr '.csv'];
            display(['writing to table: ' movieIntegrityTablePath])
            writetable(movieIntegrityTable,movieIntegrityTablePath,'FileType','text','Delimiter',',');

            movieBehaviorIntegrityTable = movieBehaviorIntegrityTable(2:end,:);
            movieBehaviorIntegrityTablePath = [obj.dataSavePath filesep obj.protocol{1} '_database_movieBehaviorIntegrityCheck_' currentDateTimeStr '.csv'];
            display(['writing to table: ' movieBehaviorIntegrityTablePath])
            writetable(movieBehaviorIntegrityTable,movieBehaviorIntegrityTablePath,'FileType','text','Delimiter',',');

        case 'duplicates'
            for thisSubjectStr=subjectList
                display(repmat('=',1,21))
                thisSubjectStr = thisSubjectStr{1};
                display(thisSubjectStr);

                validFoldersIdx = find(strcmp(thisSubjectStr,obj.subjectStr));
                % filter for folders chosen by the user
                validFoldersIdx = intersect(validFoldersIdx,fileIdxArray);
                if isempty(validFoldersIdx)
                    continue;
                end
                subjAssays = obj.assay(validFoldersIdx);
                size(subjAssays)
                size(unique(subjAssays))
                [~,idx] = unique(subjAssays);
                subjAssays(setdiff(1:length(subjAssays),idx))
                validFoldersIdx(setdiff(1:length(subjAssays),idx))
            end
            return
            % body
        case 'movieStatistics'
            movieSettings = inputdlg({...
                    'imaging movie regexp (separate multiple by comma):',...
                    'frame rate (Hz)',...
                    'frame list (empty = all)',...
                    'highpass cutoff frequency (Hz)',...
                    'frequency filter order',...
                    'picture path',...
                    'HDF5 dataset name'
                },...
                'view movie settings',1,...
                {...
                    obj.fileFilterRegexp,...
                    '5',...
                    '',...
                    '1',...
                    '3',...
                    obj.picsSavePath,...
                    obj.inputDatasetName,...
                }...
            );
            fileFilterRegexpList = strsplit(movieSettings{1},',');;
            framesPerSecond = str2num(movieSettings{2});
            frameList = str2num(movieSettings{3});
            hzCutoff = str2num(movieSettings{4});
            filterOrder = str2num(movieSettings{5});
            obj.picsSavePath = movieSettings{6};
            obj.inputDatasetName = movieSettings{7};

            obj.detailStats.linePlot = [];
            obj.detailStats.linePlot.xvalue = [];
            obj.detailStats.linePlot.yvalue = [];
            obj.detailStats.linePlot.fileFilterRegexp = {};
            obj.detailStats.linePlot.statType = {};
            obj.detailStats.linePlot.subject = {};
            obj.detailStats.linePlot.assay = {};
            obj.detailStats.linePlot.assayType = {};
            obj.detailStats.linePlot.assayNum = {};

            nRegExp = length(fileFilterRegexpList)
            for fileFilterRegexpNum = 1:nRegExp
                fileFilterRegexp = fileFilterRegexpList{fileFilterRegexpNum};
                for thisFileNumIdx = 1:nFilesToAnalyze
                    thisFileNum = fileIdxArray(thisFileNumIdx);
                    obj.fileNum = thisFileNum;
                % for thisFileNum = 1:nFiles
                    % obj.fileNum = thisFileNum;
                    display(repmat('=',1,21))
                    display([num2str(thisFileNumIdx) '/' num2str(nFilesToAnalyze) ' (' num2str(thisFileNum) '/' num2str(nFiles) '): ' obj.fileIDNameArray{obj.fileNum}]);

                    % if ~isempty(options.videoDir)
                    %   vidList = getFileList(options.videoDir,trialRegExp);
                    %   if isempty(vidList)
                    %       display(['cannot find movie: ' options.videoDir '\' trialRegExp]);
                    %       return
                 %        else
                 %            display(['going to load: ' options.videoDir '\' trialRegExp]);
                 %        end
                 %    end
                    movieList = getFileList(obj.inputFolders{obj.fileNum},fileFilterRegexp);

                    % get the movie
                    primaryMovie = loadMovieList(movieList,'convertToDouble',0,'frameList',frameList,'inputDatasetName',obj.inputDatasetName);

                    primaryMovie = single(primaryMovie);
                    testFrame = primaryMovie(:,:,1);
                    % if nanmean(testFrame(:))>0.5&&~strcmp(class(primaryMovie),'uint16')
                    %     display('zero centering...')
                    %     primaryMovie = primaryMovie-1;
                    % end

                    display(['std: ' num2str(nanstd(primaryMovie(:)))])
                    [figHandle, ~] = openFigure(177+thisFileNumIdx*fileFilterRegexpNum, '');
                        imagesc(nanstd(primaryMovie,[],3));
                        colorbar
                    obj.modelSaveImgToFile([],'movieStatisticsStd_','current',[]);

                    % [primaryMovie options] = getCurrentMovie(movieList,options);
                    lengthMovie = size(primaryMovie,3)/framesPerSecond/60;
                    primaryMovie = reshape(primaryMovie,[],size(primaryMovie,3),1);
                    % D = nanvar(C,[],1);
                    movieAvg = squeeze(nanmean(primaryMovie,1));
                    movieMax = squeeze(nanmax(primaryMovie,[],1));
                    movieMin = squeeze(nanmin(primaryMovie,[],1));
                    secondList = (1:size(primaryMovie,2))/framesPerSecond/60;
                    frameList2 = 1:size(primaryMovie,2);
                    primaryMovie = primaryMovie(round(1*end/3):round(2*end/3),:);
                    try
                        movieVar = squeeze(nanvar(primaryMovie,[],1));
                    catch
                        movieVar = 0;
                    end
                    clear primaryMovie;

                    [figHandle, ~] = openFigure(17799+thisFileNumIdx*fileFilterRegexpNum, '');
                    plot(frameList2,movieAvg,'k'); hold on;
                    plot(frameList2,movieMax);
                    plot(frameList2,movieMin);
                    % title(['mean']);
                    ylabel('pixel value');xlabel('frame'); box off;
                    axis tight;
                    % set(gca,'xlim',[0 lengthMovie]);
                    legend({'mean','max','min'},'best')
                    % save images
                    set(figHandle,'PaperUnits','inches','PaperPosition',[0 0 16 9])
                    % figure(figHandle)
                    obj.modelSaveImgToFile([],'movieStatisticsMeanMinMax_','current',[]);

                    [figHandle, ~] = openFigure(1776+thisFileNumIdx*fileFilterRegexpNum, '');
                    subplot(3,1,1)
                    plot(secondList,movieAvg,'k'); hold on;
                    plot(secondList,movieMax);
                    plot(secondList,movieMin);
                    % title(['mean']);
                    ylabel('pixel value');xlabel('minutes'); box off;
                    axis tight;
                    set(gca,'xlim',[0 lengthMovie]);
                    legend({'mean','max','min'})

                    % movieAvg = movieAvg/nanmean(movieAvg);
                    movieAvg = movieAvg - nanmean(movieAvg);

                    fs = framesPerSecond;
                    Wn = (hzCutoff/(fs/2));
                    N = filterOrder;
                    [B,A] = butter(N,Wn,'high');
                    movieAvg = filter(B,A,movieAvg);

                    subplot(3,1,2)
                    plot(secondList,movieVar,'b')
                    % title('variance');
                    ylabel('variance');xlabel('minutes'); box off;
                    axis tight;
                    set(gca,'xlim',[0 lengthMovie]);
                    % suptitle(thisDirSaveStr);
                    supHandle = suptitle([strrep(obj.folderBaseSaveStr{obj.fileNum},'_',' ') ' | ' fileFilterRegexp]);
                    set(supHandle, 'Interpreter', 'none');

                    % plot the PSD of the whole movie average
                    subplot(3,1,3)
                    % [pxx,f] = periodogram(movieAvg,[],framesPerSecond*10,framesPerSecond);
                    [pxx,f] = periodogram(movieAvg,[],[],framesPerSecond);
                    size(pxx)
                    size(f)
                    warning off;
                    plot(f,10*log10(pxx./f),'r');
                    box off; axis tight;
                    xlabel('frequency (Hz)');ylabel('Power spectral density (dB/Hz)')

                    % save images
                    set(figHandle,'PaperUnits','inches','PaperPosition',[0 0 16 9])
                    % figure(figHandle)
                    obj.modelSaveImgToFile([],'movieStatistics_','current',[]);
                    warning on;
                    % tmpDirPath = strcat(obj.picsSavePath,filesep,'movieStatistics',filesep);
                    % mkdir(tmpDirPath);
                    % saveFile = char(strrep(strcat(tmpDirPath,trialRegExp,'.png'),'/',''));
                    % saveas(gcf,saveFile);
                    numPtsToAdd = length(f(:));
                    obj.detailStats.linePlot.xvalue(end+1:end+numPtsToAdd,1) = f;
                    obj.detailStats.linePlot.yvalue(end+1:end+numPtsToAdd,1) = pxx./f;
                    obj.detailStats.linePlot.fileFilterRegexp(end+1:end+numPtsToAdd,1) = {fileFilterRegexp};
                    obj.detailStats.linePlot.statType(end+1:end+numPtsToAdd,1) = {'powerSpectralDensity'};
                    obj.detailStats.linePlot.subject(end+1:end+numPtsToAdd,1) = {obj.subjectStr{obj.fileNum}};
                    obj.detailStats.linePlot.assay(end+1:end+numPtsToAdd,1) = {obj.assay{obj.fileNum}};
                    obj.detailStats.linePlot.assayType(end+1:end+numPtsToAdd,1) = {obj.assayType{obj.fileNum}};
                    obj.detailStats.linePlot.assayNum(end+1:end+numPtsToAdd,1) = {obj.assayNum{obj.fileNum}};
                end
            end

            savePath = [obj.dataSavePath obj.protocol{obj.fileNum} '_movieStatistics.tab'];
            display(['saving data to: ' savePath])
            writetable(struct2table(obj.detailStats.linePlot),savePath,'FileType','text','Delimiter','\t');

            return
            % body
        otherwise
            body
    end
    return
    function local_getVideoRegexp()
        switch videoTrialRegExpIdx
            case 1
                videoTrialRegExp = [obj.date{obj.fileNum} '_' obj.protocol{obj.fileNum} '_' obj.fileIDArray{obj.fileNum}];
            case 2
                dateTmp = strsplit(obj.date{obj.fileNum},'_');
                videoTrialRegExp = strcat(dateTmp{1}(end-1:end),dateTmp{2},dateTmp{3},'-',obj.subjectStr{obj.fileNum},'-',obj.assay{obj.fileNum});
            case 3
                dateTmp = strsplit(obj.date{obj.fileNum},'_');
                videoTrialRegExp = strcat(dateTmp{1}(end-1:end),dateTmp{2},dateTmp{3},'_',obj.subjectStr{obj.fileNum},'_',obj.assay{obj.fileNum});
            case 4
                videoTrialRegExp = [obj.subjectStr{obj.fileNum} '_' obj.assay{obj.fileNum}]
            case 5
                dateTmp = strsplit(obj.date{obj.fileNum},'_');
                videoTrialRegExp = strcat(dateTmp{1}(end-1:end),dateTmp{2},dateTmp{3},'_',obj.subjectStr{obj.fileNum});
            otherwise
                videoTrialRegExp = fileFilterRegexp
        end
    end
end