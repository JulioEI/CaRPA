function obj = modelModifyMovies(obj)
	% either modify an existing movie (e.g. to remove edges) or store a temporary file to be loaded later for other files
	% biafra ahanonu
	% started: 2016.02.19
	% inputs
		%
	% outputs
		%

	% changelog
		%
	% TODO
		%

	try
		% get user input
		movieSettings = inputdlg({...
				'start:end frames for preview (blank = all frames)',...
				'start:end frames for saving (blank = all frames)',...
				'file regexp:',...
				'replacement file regexp:',...
				'analyze specific folder (leave blank to copy to same folder)',...
				'input HDF5 dataset name',...
				'output HDF5 dataset name'...
			},...
			'copy files to /archive/ folder',1,...
			{...
				'1:50',...
				'',...
				obj.fileFilterRegexp,...
				'manualCut',...
				'',...
				obj.inputDatasetName,...
				obj.inputDatasetName...
			}...
		);
		setNo = 1;
		frameList = str2num(movieSettings{setNo});setNo = setNo+1;
		frameListSave = str2num(movieSettings{setNo});setNo = setNo+1;
		fileFilterRegexp = movieSettings{setNo};setNo = setNo+1;
		replaceFileFilterRegexp = movieSettings{setNo};setNo = setNo+1;
		saveSpecificFolder = movieSettings{setNo};setNo = setNo+1;
		inputDatasetName =  movieSettings{setNo};setNo = setNo+1;
		outputDatasetName =  movieSettings{setNo};setNo = setNo+1;
		% analyzeMovieFiles =  str2num(movieSettings{6});
		obj.fileFilterRegexp = replaceFileFilterRegexp;

        % get files to analyze
		[fileIdxArray idNumIdxArray nFilesToAnalyze nFiles] = obj.getAnalysisSubsetsToAnalyze();

		% start Miji
		Miji;

		% loop over all directories, get masks from user first then batch
		% cropping whole movies
        movieMaskArray = {};
		for thisFileNumIdx = 1:nFilesToAnalyze
			thisFileNum = fileIdxArray(thisFileNumIdx);
			obj.fileNum = thisFileNum;
			display(repmat('=',1,21))
			display([num2str(thisFileNumIdx) '/' num2str(nFilesToAnalyze) ' (' num2str(thisFileNum) '/' num2str(nFiles) '): ' obj.fileIDNameArray{obj.fileNum}]);
			% =====================
			try
				movieList = getFileList(obj.inputFolders{obj.fileNum}, fileFilterRegexp);
				moviePath = movieList{1};
				[frameListTmp] = subfxnVerifyFrameList(frameList,moviePath,inputDatasetName);
				[primaryMovie] = loadMovieList(moviePath,'convertToDouble',0,'frameList',frameListTmp(:),'inputDatasetName',inputDatasetName);

				% ranMovie = rand([1000 1000 20]);
				MIJ.createImage(obj.folderBaseSaveStr{obj.fileNum}, primaryMovie, true);
	            for foobar=1:2; MIJ.run('In [+]'); end
	        	for foobar=1:2; MIJ.run('Enhance Contrast','saturated=0.35'); end
				% uiwait(msgbox('select region of movie to keep','Success','modal'));
	            movieDecision = questdlg('Should movie be cropped? If YES, draw ROI of area to keep.', ...
						'Movie decision', ...
						'yes','no','yes');
	            if strcmp(movieDecision,'yes')
	                try
	                    MIJ.run('Set Slice...', 'slice=1');
	                    MIJ.run('Set...', 'value=1');
	                    MIJ.run('Make Inverse');
	                    MIJ.run('Set...', 'value=0');
	                    MIJ.run('Select None');
	                    MIJ.run('Make Substack...', 'delete slices=1');
	                    movieMaskArray{thisFileNumIdx} = MIJ.getCurrentImage;
	                    MIJ.run('Close All Without Saving');
	                catch err
	                    movieMaskArray{thisFileNumIdx} = [];
	                    try
	                        MIJ.run('Close All Without Saving');
	                    catch
	                    end
	                    display(repmat('@',1,7))
	                    disp(getReport(err,'extended','hyperlinks','on'));
	                    display(repmat('@',1,7))
	                end
	            else
	                MIJ.run('Close All Without Saving');
	                movieMaskArray{thisFileNumIdx} = [];
	            end
	       catch err
	           movieMaskArray{thisFileNumIdx} = [];
	           display(repmat('@',1,7))
	           disp(getReport(err,'extended','hyperlinks','on'));
	           display(repmat('@',1,7))
	       end
        end

		MIJ.exit;

        % go through and crop each movie then save
        display(repmat('*',1,21))
        display(repmat('*',1,21))
        display('CROPPING THE MOVIES!');
		for thisFileNumIdx = 1:nFilesToAnalyze
			thisFileNum = fileIdxArray(thisFileNumIdx);
			obj.fileNum = thisFileNum;
			display(repmat('=',1,21))
			display([num2str(thisFileNumIdx) '/' num2str(nFilesToAnalyze) ' (' num2str(thisFileNum) '/' num2str(nFiles) '): ' obj.fileIDNameArray{obj.fileNum}]);

            % get previously stored mask
            movieMask = movieMaskArray{thisFileNumIdx};
            if isempty(movieMask)
               continue;
            end
            movieMask(isnan(movieMask)) = 1;
            movieMask = logical(movieMask);

            % load entire movie and crop
            movieList = getFileList(obj.inputFolders{obj.fileNum}, fileFilterRegexp);
			moviePath = movieList{1};
			[frameListSaveTmp] = subfxnVerifyFrameList(frameListSave,moviePath,inputDatasetName);
			[primaryMovie] = loadMovieList(moviePath,'convertToDouble',0,'frameList',frameListSaveTmp,'inputDatasetName',inputDatasetName);
			primaryMovie = bsxfun(@times,primaryMovie,movieMask);

			% save the file in the new location
			[PATHSTR,NAME,EXT] = fileparts(moviePath);
			% newPathFile = [PATHSTR filesep NAME '_manualCrop' EXT];
			newPathFile = [PATHSTR filesep strrep(NAME,fileFilterRegexp,replaceFileFilterRegexp) EXT];
			[output] = writeHDF5Data(primaryMovie,newPathFile,'datasetname',				outputDatasetName);
        end

	catch err
		display(repmat('@',1,7))
		disp(getReport(err,'extended','hyperlinks','on'));
		display(repmat('@',1,7))
	end

end
function [frameListTmp] = subfxnVerifyFrameList(frameList,moviePath,inputDatasetName)
	if isempty(frameList)
		frameListTmp = frameList;
	else
		movieDims = loadMovieList(moviePath,'convertToDouble',0,'frameList',[],'inputDatasetName',inputDatasetName,'treatMoviesAsContinuous',1,'loadSpecificImgClass','single','getMovieDims',1);
		frameListTmp = frameList;
		% remove frames that are too large
		frameListTmp(frameListTmp>=movieDims.z) = [];
	end
end