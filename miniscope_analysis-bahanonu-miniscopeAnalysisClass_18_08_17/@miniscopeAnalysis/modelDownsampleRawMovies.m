function obj = modelDownsampleRawMovies(obj)
    % downsamples inscopix files and moves them to appropriate folders.
    % biafra ahanonu
    % inputs
        %
    % outputs
        %

    % changelog
        % 2015.11.25 - Added ability to downsample TIFFs, at the moment assumes that they will be done the inscopix way.

	try
		downsampleSettings = inputdlg(...
			{...
				'Folder(s) where raw HDF5s are located. Use comma to separate multiple source folders:',...
				'Folder to save downsampled HDF5s to:',...
				'Decompression source root folder(s). Use comma to separate multiple source folders:',...
				'Downsample factor (in x-y):',...
				'Regexp for HDF5/TIFF files:',...
				'HDF5 hierarchy name where movie is stored:',...
				'Max chunk size (MB)',...
				'Regexp for folders:',...
				'Regexp for base filename (use txt log name):',...
				'Extension for base filename:',...
				'Location to store raw HDF5 (leave empty to skip):',...
				'Raw downsample factor (in x-y):',...
                'Output HDF5 dataset name:'
			},...
			'downsample settings',[1 100],...
			{...
                obj.downsampleRawOptions.folderListInfo,...
                obj.downsampleRawOptions.downsampleSaveFolder,...
                obj.downsampleRawOptions.downsampleSrcFolder,...
                obj.downsampleRawOptions.downsampleFactor,...
                obj.downsampleRawOptions.fileFilterRegexp,...
                obj.downsampleRawOptions.datasetName,...
                obj.downsampleRawOptions.maxChunkSize,...
                obj.downsampleRawOptions.srcFolderFilterRegexp,...
                obj.downsampleRawOptions.srcSubfolderFileFilterRegexp,...
                obj.downsampleRawOptions.srcSubfolderFileFilterRegexpExt,...
                obj.downsampleRawOptions.downsampleSaveFolderTwo,...
                obj.downsampleRawOptions.downsampleFactorTwo,...
                obj.downsampleRawOptions.outputDatasetName,...
            });
				% 'A:\data\processing\',...
				% 'B:\data\processing\',...
				% 'E:\data\raw\',...
				% '4',...
				% 'recording.*.hdf5',...
				% '/images',...
				% '25000',...
				% '201\d',...
				% 'recording.*.(txt|xml)',...
				% '(.txt|.xml)',...
				% 'E:\data\raw\',...
				% '2',...
    %             '/1',...
        obj.downsampleRawOptions.folderListInfo = downsampleSettings{1};
        obj.downsampleRawOptions.downsampleSaveFolder = downsampleSettings{2};
        obj.downsampleRawOptions.downsampleSrcFolder = downsampleSettings{3};
        obj.downsampleRawOptions.downsampleFactor = downsampleSettings{4};
        obj.downsampleRawOptions.fileFilterRegexp = downsampleSettings{5};
        obj.downsampleRawOptions.datasetName = downsampleSettings{6};
        obj.downsampleRawOptions.maxChunkSize = downsampleSettings{7};
        obj.downsampleRawOptions.srcFolderFilterRegexp = downsampleSettings{8};
        obj.downsampleRawOptions.srcSubfolderFileFilterRegexp = downsampleSettings{9};
        obj.downsampleRawOptions.srcSubfolderFileFilterRegexpExt = downsampleSettings{10};
        obj.downsampleRawOptions.downsampleSaveFolderTwo = downsampleSettings{11};
        obj.downsampleRawOptions.downsampleFactorTwo = downsampleSettings{12};
        obj.downsampleRawOptions.outputDatasetName = downsampleSettings{13};

		folderListInfo = strsplit(downsampleSettings{1},',');
        cellfun(@display,folderListInfo);
		nFolders = length(folderListInfo);
        % ~isempty(regexp(downsampleOptions.fileFilterRegexp,'(.hdf5|.h5)'))
        % ~isempty(regexp(downsampleOptions.fileFilterRegexp,'(.tiff|.tif)'))
        % ============================
        for folderNo = 1:nFolders
            display(repmat('=',1,21))
            display([num2str(folderNo) '\' num2str(nFolders) ': ' folderListInfo{folderNo}])
            % downsample if all decompressed files are in the same folder
            downsampleOptions.downsampleSaveFolder = downsampleSettings{2};
            downsampleOptions.datasetName = downsampleSettings{6};
            downsampleOptions.maxChunkSize = str2num(downsampleSettings{7});
            downsampleOptions.downsampleFactor = str2num(downsampleSettings{4});
            downsampleOptions.fileFilterRegexp = downsampleSettings{5};
            downsampleOptions.outputDatasetName = downsampleSettings{13};

            downsampleOptions.downsampleSaveFolderTwo = downsampleSettings{11};
            downsampleOptions.downsampleFactorTwo = str2num(downsampleSettings{12});

            movieList = getFileList([folderListInfo{folderNo} filesep], downsampleOptions.fileFilterRegexp);

            if ~isempty(regexp(downsampleOptions.fileFilterRegexp,'(.hdf5|.h5)'))
                display('Downsampling HDF5s...')
                downsampleHDFMovieFxnObj(movieList,downsampleOptions);
            elseif ~isempty(regexp(downsampleOptions.fileFilterRegexp,'(.tiff|.tif)'))
                display('Downsampling TIFFs...')
                folderPath = folderListInfo{folderNo};
                downsampleOptions.srcSubfolderFileFilterRegexp = downsampleSettings{9};
                downsampleOptions.srcSubfolderFileFilterRegexpExt = downsampleSettings{10};
                downsampleTiffMovieFxnObj(folderPath,downsampleOptions);
            end
            % ioptions.folderListInfo = [folderListInfo{folderNo} filesep];
            % ioptions.downsampleSaveFolder = [downsampleSettings{2} filesep];
            % ioptions.downsampleFactor = str2num(downsampleSettings{4});
            % ioptions.fileFilterRegexp = downsampleSettings{5};
            % ioptions.datasetName = downsampleSettings{6};
            % ioptions.maxChunkSize = str2num(downsampleSettings{7});
            % ioptions.runArg = 'downsampleMovie';
            % ostruct = controllerAnalysis('options',ioptions);
        end
		% re-create folder structure
		clear ioptions
        % run HDF5 or TIFF
		%
        downsampleSrcFolder = strsplit(downsampleSettings{3},',');
		downsampleSaveFolder = [downsampleSettings{2} filesep];
		% used to determine which folders to copy from src to dest
		ioptions.srcFolderFilterRegexp = downsampleSettings{8};
		% this regexp is used to search the destination directory
		ioptions.srcSubfolderFileFilterRegexp = downsampleSettings{9};
		%
		ioptions.srcSubfolderFileFilterRegexpExt = downsampleSettings{10};
		[success destFolders] = moveFilesToFolders(downsampleSrcFolder,char(downsampleSaveFolder(:))','options',ioptions);
		destFolders = unique(destFolders);
        % ============================
        if ~isempty(downsampleSettings{11})
            % movie second set of movie files
            downsampleSaveFolder = [downsampleSettings{11} filesep];
            % used to determine which folders to copy from src to dest
            ioptions.srcFolderFilterRegexp = downsampleSettings{8};
            % this regexp is used to search the destination directory
            ioptions.srcSubfolderFileFilterRegexp = downsampleSettings{9};
            %
            ioptions.srcSubfolderFileFilterRegexpExt = downsampleSettings{10};
            [success destFolders2] = moveFilesToFolders(downsampleSrcFolder,char(downsampleSaveFolder(:))','options',ioptions);
        end
        % ============================
        % if ~isempty(downsampleSettings{11})
        %     % now same the semi-downsampled raw files
        %     for folderNo = 1:nFolders
        %         % downsample if all decompressed files are in the same folder
        %         downsampleOptions.downsampleSaveFolder = [downsampleSettings{11} filesep];
        %         downsampleOptions.datasetName = downsampleSettings{6};
        %         downsampleOptions.maxChunkSize = str2num(downsampleSettings{7});
        %         downsampleOptions.downsampleFactor = str2num(downsampleSettings{12});
        %         downsampleOptions.fileFilterRegexp = downsampleSettings{5};

        %         movieList = getFileList([folderListInfo{folderNo} filesep], downsampleOptions.fileFilterRegexp);

        %         downsampleHDFMovieFxnObj(movieList,downsampleOptions);
        %     end
        %     % re-create folder structure
        %     clear ioptions
        %     %
        %     downsampleSaveFolder = [downsampleSettings{11} filesep];
        %     % used to determine which folders to copy from src to dest
        %     ioptions.srcFolderFilterRegexp = downsampleSettings{8};
        %     % this regexp is used to search the destination directory
        %     ioptions.srcSubfolderFileFilterRegexp = downsampleSettings{9};
        %     %
        %     ioptions.srcSubfolderFileFilterRegexpExt = downsampleSettings{10};
        %     [success destFolders2] = moveFilesToFolders(strsplit(downsampleSettings{3},','),char(downsampleSaveFolder(:))','options',ioptions);
        %     % destFolders = unique(destFolders);
        % end
        % destFolders
        % class(destFolders)

        % add folder info to inputFolders class property
		if isempty(obj.inputFolders)
            if iscell(destFolders)
                nAddFolders = length(destFolders);
                for addFolderNo = 1:nAddFolders
                    obj.inputFolders{addFolderNo,1} = destFolders{addFolderNo};
                end
            else
                nAddFolders = 1;
                obj.inputFolders{1,1} = destFolders;
            end
			obj.dataPath = obj.inputFolders;
		else
            if iscell(destFolders)
                nAddFolders = length(destFolders);
                for addFolderNo = 1:nAddFolders
                    obj.inputFolders{end+addFolderNo,1} = destFolders{addFolderNo};
                end
            else
                nAddFolders = 1;
                obj.inputFolders{end+1,1} = destFolders;
            end
            % obj.inputFolders = [obj.inputFolders(:); destFolders(:)];
			% obj.inputFolders = cat(obj.inputFolders,destFolders);
            obj.dataPath = obj.inputFolders;
			% obj.dataPath = cat(obj.dataPath,destFolders);
		end
	catch err
		obj.foldersToAnalyze = [];
		display(repmat('@',1,7))
		disp(getReport(err,'extended','hyperlinks','on'));
		display(repmat('@',1,7))
	end
end
function downsampleTiffMovieFxnObj(folderPath,options)

    % get list of regular expressions
    fileRegexpListOriginal = getFileList([folderPath filesep],options.srcSubfolderFileFilterRegexp);
    fileRegexpList = regexp(fileRegexpListOriginal, options.srcSubfolderFileFilterRegexp,'match');
    cellfun(@display,fileRegexpList);
    % fileRegexpList = getFileList([folderPath filesep], options.srcFolderFilterRegexp);
    nFileRegexp = length(fileRegexpList);
    % go over each set of regexp
    for fileRegexpNo=1:nFileRegexp
        thisFileRegexp = fileRegexpList{fileRegexpNo};
        thisFileRegexp = regexprep(thisFileRegexp,options.srcSubfolderFileFilterRegexpExt,'');
        thisFileRegexp = thisFileRegexp{1};
        display(thisFileRegexp)
        % get a list of tiff associated with regexp
        [folderPath filesep]
        [thisFileRegexp '.*.tif']
        movieList = getFileList([folderPath filesep], [thisFileRegexp '.*.tif']);
        % due to naming convention, shift last to 1st
        movieList = circshift(movieList,[0 1]);
        cellfun(@display,movieList);
        % get path string
        [pathstr,name,ext] = fileparts(fileRegexpListOriginal{fileRegexpNo});
        newFilename = [options.downsampleSaveFolder filesep 'concat_' name '.h5'];
        display(['save HDF5: ' newFilename])
        % load and downsample each tiff
        nMovies = length(movieList);
        for movieNo=1:nMovies
            display(repmat('+',1,21))
            display(['downsampling ' num2str(movieNo) '/' num2str(nMovies)])
            inputFilePath = movieList{movieNo};
            display(['input: ' inputFilePath]);
            % [pathstr,name,ext] = fileparts(inputFilePath);

            % downsampleFilename = [pathstr '\concat_' name '.h5']
            doptions.downsampleDimension = 'space';
            doptions.downsampleType = 'bilinear';
            doptions.downsampleFactor = options.downsampleFactor;
            [inputMovie] = downsampleMovie(inputFilePath, 'options',doptions);
            % save the movie
            % save or append each tiff to corresponding HDF5 file
            if movieNo==1
                createHdf5File(newFilename, options.outputDatasetName, inputMovie);
            else
                appendDataToHdf5(newFilename, options.outputDatasetName, inputMovie);
            end
        end
        srcFilenameTxt = [pathstr filesep name '.txt'];
        srcFilenameXml = [pathstr filesep name '.xml'];
        destFilenameTxt = [options.downsampleSaveFolder filesep name '.txt']
        destFilenameXml = [options.downsampleSaveFolder filesep name '.xml']
        if exist(srcFilenameTxt,'file')
            display([srcFilenameTxt ' > ' destFilenameTxt])
            copyfile(srcFilenameTxt,destFilenameTxt)
        elseif exist(srcFilenameXml,'file')
            display([srcFilenameXml ' > ' destFilenameXml])
            copyfile(srcFilenameXml,destFilenameXml)
        end
    end
end
function downsampleHDFMovieFxnObj(movieList,options)
    % downsamples an HDF5 movie, normally the raw recording files

	% display(movieList)
    cellfun(@display,movieList);
    nMovies = length(movieList);
	for movieNo=1:nMovies
        display(repmat('+',1,21))
        display(['downsampling ' num2str(movieNo) '/' num2str(nMovies)])
		inputFilePath = movieList{movieNo};
		display(['input: ' inputFilePath]);
        [pathstr,name,ext] = fileparts(inputFilePath);
        downsampleFilename = [pathstr '\concat_' name '.h5']
        srcFilenameTxt = [pathstr filesep name '.txt']
        srcFilenameXml = [pathstr filesep name '.xml']

        % currentDateTimeStr = datestr(now,'yyyymmdd_HHMM','local');
        % mkdir([thisDir filesep 'processing_info'])
        % diarySaveStr = [thisDir filesep 'processing_info' filesep currentDateTimeStr '_preprocess.log'];
        % display(['saving diary: ' diarySaveStr])
        % diary(diarySaveStr);

        try
            % see if file already exists at the destination
	        if ~exist(downsampleFilename,'file')
                % check whether to downsample in different location
	        	if isempty(options.downsampleSaveFolder)
	        		downsampleHdf5Movie(inputFilePath, 'inputDatasetName', options.datasetName, 'maxChunkSize', options.maxChunkSize,'downsampleFactor',options.downsampleFactor,'outputDatasetName',options.outputDatasetName,'saveFolderTwo',options.downsampleSaveFolderTwo,'downsampleFactorTwo',options.downsampleFactorTwo);
	        	else
	        		downsampleHdf5Movie(inputFilePath, 'inputDatasetName', options.datasetName, 'maxChunkSize', options.maxChunkSize,'saveFolder',options.downsampleSaveFolder,'downsampleFactor',options.downsampleFactor,'outputDatasetName',options.outputDatasetName,'saveFolderTwo',options.downsampleSaveFolderTwo,'downsampleFactorTwo',options.downsampleFactorTwo);

                    % copy information files to new downsample folder location
	        		destFilenameTxt = [options.downsampleSaveFolder filesep name '.txt']
	        		destFilenameXml = [options.downsampleSaveFolder filesep name '.xml']
		        	if exist(srcFilenameTxt,'file')
		        		copyfile(srcFilenameTxt,destFilenameTxt)
		        	elseif exist(srcFilenameXml,'file')
		        		copyfile(srcFilenameXml,destFilenameXml)
	        		end
                    if ~isempty(options.downsampleSaveFolderTwo)
                        destFilenameTxt = [options.downsampleSaveFolderTwo filesep name '.txt']
                        destFilenameXml = [options.downsampleSaveFolderTwo filesep name '.xml']
                        if exist(srcFilenameTxt,'file')
                            copyfile(srcFilenameTxt,destFilenameTxt)
                        elseif exist(srcFilenameXml,'file')
                            copyfile(srcFilenameXml,destFilenameXml)
                        end
                    end
	        	end
	        elseif ~isempty(options.downsampleSaveFolder)&~exist([options.downsampleSaveFolder '\concat_' name '.h5'],'file')
	        	downsampleHdf5Movie(inputFilePath, 'inputDatasetName', options.datasetName, 'maxChunkSize', options.maxChunkSize,'saveFolder',options.downsampleSaveFolder,'downsampleFactor',options.downsampleFactor,'outputDatasetName',options.outputDatasetName,'saveFolderTwo',options.downsampleSaveFolderTwo,'downsampleFactorTwo',options.downsampleFactorTwo);
	        	destFilenameTxt = [options.downsampleSaveFolder filesep name '.txt']
	        	destFilenameXml = [options.downsampleSaveFolder filesep name '.xml']
	        	if exist(srcFilenameTxt,'file')
	        		copyfile(srcFilenameTxt,destFilenameTxt)
	        	elseif exist(srcFilenameXml,'file')
	        		copyfile(srcFilenameXml,destFilenameXml)
        		end

                % copy information files to new downsample folder location
                if ~isempty(options.downsampleSaveFolderTwo)
                    destFilenameTxt = [options.downsampleSaveFolderTwo filesep name '.txt']
                    destFilenameXml = [options.downsampleSaveFolderTwo filesep name '.xml']
                    if exist(srcFilenameTxt,'file')
                        copyfile(srcFilenameTxt,destFilenameTxt)
                    elseif exist(srcFilenameXml,'file')
                        copyfile(srcFilenameXml,destFilenameXml)
                    end
                end
	        else
	            display(['skipping: ' inputFilePath])
	        end
        catch err
        	display(repmat('@',1,7))
        	disp(getReport(err,'extended','hyperlinks','on'));
        	display(repmat('@',1,7))
        end
	end
end
function vestigialFxn()
    % ============================
    % move files to their correct folders
    % moveSettings = inputdlg({...
    %   'regexp for folders:',...
    %   'regexp for base filename (use txt log name):',...
    %   'extension for base filename:',...
    %   'source root folder(s), use comma to separate multiple source folders:',...
    %   'destination folder (where HDF5s were downsampled):'...
    %   },...
    %   'downsample settings',1,{...
    %   '2014',...
    %   'recording.*.txt',...
    %   '.txt',...
    %   'E:\',...
    %   char(ioptions.downsampleSaveFolder(:))'});
    % re-create folder structure
    % used to determine which folders to copy from src to dest
    % ioptions.srcFolderFilterRegexp = moveSettings{1};
    % % this regexp is used to search the destination directory
    % ioptions.srcSubfolderFileFilterRegexp = moveSettings{2};
    % %
    % ioptions.srcSubfolderFileFilterRegexpExt = moveSettings{3};
    % [success destFolders] = moveFilesToFolders(strsplit(moveSettings{4},','),moveSettings{5},'options',ioptions);
end