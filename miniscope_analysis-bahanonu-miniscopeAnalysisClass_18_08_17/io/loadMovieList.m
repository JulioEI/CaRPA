function [outputMovie movieDims nPixels nFrames] = loadMovieList(movieList, varargin)
	% load movies, automatically detects type (avi, tif, or hdf5) and concatenates if multiple movies in a list
	% biafra ahanonu
	% started: 2013.11.01
	% inputs
	% 	movieList = either a char string containing a path name or a cell array containing char strings, e.g. 'pathToFile' or {'path1','path2'}
	% outputs
	% 	outputMovie
	% 	movieDims
	% 	nPixels
	% 	nFrames
	% options
	% options.supportedTypes = {'.h5','.hdf5','.tif','.tiff','.avi'};
	% % movie type
	% options.movieType = 'tiff';
	% % hierarchy name in hdf5 where movie is
	% options.inputDatasetName = '/1';
	% % convert file movie to double?
	% options.convertToDouble = 0;
	% % 'single','double'
	% options.loadSpecificImgClass = [];
	% % list of specific frames to load
	% options.frameList = [];
	% % should the waitbar be shown?
	% options.waitbarOn=1;
	% % just return the movie dimensions
	% options.getMovieDims = 0;
	% % treat movies in list as continuous with regards to frame
	% options.treatMoviesAsContinuous = 0;
	% NOTE: assume 3D movies with [x y frames] as dimensions, if movies are different sizes, use largest dimensions and align all movies to top-left corner

	% changelog
		% 2014.02.14 [14:14:39] now can load non-monotonic lists for avi and hdf5 files.
		% 2014.03.27 - several updates to speed up function, fixed several assumption issues (all movies same file type, etc.) and brought name scheme in line with other fxns
		% 2015.02.25 [15:32:10] fixed bug pertaining to treatMoviesAsContinuous not working properly if frameList was blank and only a single movie was input.
		% 2015.05.28 [02:39:51] bug fix with treatMoviesAsContinuous, the global frames weren't quite correct
		% 2016.02.23 [19:40:31] added in proper-ish ability to read in RGB grayscale AVI files
		% 2017.07.06 [20:13:53] improved getHdf5Info() subfxn to handle Jessica's HDF5 format in a more automatic fashion, e.g. if there is behavior or other group level data in the HDF5 file
	% TODO
		% Allow fallbacks for HDF5 dataset name, e.g. if can't find /1, look for /images
		% MAKE tiff loading recognize frameList input
		% add preallocation by pre-reading each movie's dimensions - DONE
		% determine file type by properties of file instead of extension (don't trust input...)
		% remove need to use tmpMovie....
		% verify movies are of supported load types, remove from list if not and alert user, should be an option (e.g. return instead) - DONE
		% allow user to input frames that are global across several files, e.g. [1:500 1:200 1:300] are the lengths of each movie, so if input [650:670] in frameList, should grab 150:170 from movie 2
		% add ability to degrade gracefully with HDF5 dataset names, so try several backup datasetnames if one doesn't work

	% ========================
	options.supportedTypes = {'.h5','.hdf5','.tif','.tiff','.avi'};
	% movie type
	options.movieType = 'tiff';
	% hierarchy name in hdf5 where movie is
	options.inputDatasetName = '/1';
	% fallback hierarchy name, e.g. '/images'
	options.inputDatasetNameBackup = [];
	% convert file movie to double?
	options.convertToDouble = 0;
	% 'single','double'
	options.loadSpecificImgClass = [];
	% list of specific frames to load
	options.frameList = [];
	% should the waitbar be shown?
	options.waitbarOn = 1;
	% just return the movie dimensions
	options.getMovieDims = 0;
	% treat movies in list as continuous with regards to frame
	options.treatMoviesAsContinuous = 0;
	% whether to display info
	options.displayInfo = 1;
	% pre-specify the size, if need to get around memory re-allocation issues
	options.presetMovieSize = [];
	% get options
	options = getOptions(options,varargin);
	% unpack options into current workspace
	% fn=fieldnames(options);
	% for i=1:length(fn)
	%     eval([fn{i} '=options.' fn{i} ';']);
	% end

	% ========================
	% allow usr to input just a string if a single movie
    if strcmp(class(movieList),'char')
        movieList = {movieList};
    end

	% ========================
	% remove unsupported files
	for iMovie=1:length(movieList)
		thisMoviePath = movieList{iMovie};
		[options.movieType supported] = getMovieFileType(thisMoviePath);
		if supported==0
			subfxnDisplay(['removing unsupported file from list: ' thisMoviePath],options);
		else
			tmpMovieList{iMovie} = movieList{iMovie};
		end
	end
	% if tmp doesn't exist, means no input files are valid, return
	if exist('tmpMovieList','var')
		movieList = tmpMovieList;
	else
		outputMovie = NaN;
		movieDims = NaN;
		nPixels = NaN;
		nFrames = NaN;
		return;
	end
	numMovies = length(movieList);

    % ========================
	% pre-read each file to allow pre-allocation of output file
	reverseStr = '';
	for iMovie=1:numMovies
		thisMoviePath = movieList{iMovie};
		[options.movieType supported] = getMovieFileType(thisMoviePath);
		if supported==0

		end
		switch options.movieType
			case 'tiff'
				tiffHandle = Tiff(thisMoviePath, 'r');
				tmpFrame = tiffHandle.read();
				tiffHandle.close(); clear tiffHandle
				xyDims=size(tmpFrame);

				dims.x(iMovie) = xyDims(1);
				dims.y(iMovie) = xyDims(2);
				dims.z(iMovie) = size(imfinfo(thisMoviePath),1);

				if dims.z(iMovie)==1
					fileInfo = imfinfo(thisMoviePath);
					try
						numFramesStr = regexp(fileInfo.ImageDescription, 'images=(\d*)', 'tokens');
					    nFrames = str2double(numFramesStr{1}{1});
					catch
						nFrames = 1;
					end
				    dims.z(iMovie) = nFrames;
				end
			case 'hdf5'
				%
				hinfo = hdf5info(thisMoviePath);
                [hReadInfo thisDatasetName] = getHdf5Info();
                hReadInfo
                hReadInfo.Name
				dims.x(iMovie) = hReadInfo.Dims(1);
				dims.y(iMovie) = hReadInfo.Dims(2);
				dims.z(iMovie) = hReadInfo.Dims(3);
				dims.one(iMovie) = hReadInfo.Dims(1);
				dims.two(iMovie) = hReadInfo.Dims(2);
				dims.three(iMovie) = hReadInfo.Dims(3);
                options.inputDatasetName = thisDatasetName;%CUSTOM_CODE!
				if ischar(options.inputDatasetName)
					tmpFrame = readHDF5Subset(thisMoviePath,[0 0 1],[dims.x(iMovie) dims.y(iMovie) 1],'datasetName',options.inputDatasetName,'displayInfo',options.displayInfo);
				else
					tmpFrame = readHDF5Subset(thisMoviePath,[0 0 1],[dims.x(iMovie) dims.y(iMovie) 1],'datasetName',thisDatasetName,'displayInfo',options.displayInfo);
				end
			case 'avi'
				xyloObj = VideoReader(thisMoviePath);
				dims.x(iMovie) = xyloObj.Height;
				dims.y(iMovie) = xyloObj.Width;
				dims.z(iMovie) = xyloObj.NumberOfFrames;
				tmpFrame = read(xyloObj, 1);
		end
		if isempty(options.loadSpecificImgClass)
			imgClass = class(tmpFrame);
		else
			imgClass = options.loadSpecificImgClass;
		end
		% change dims.z if user specifies a list of frames
		if (~isempty(options.frameList)|options.frameList>dims.z(iMovie))&options.treatMoviesAsContinuous==0
			dims.z(iMovie) = length(options.frameList);
		end
		if options.displayInfo==1
			reverseStr = cmdWaitbar(iMovie,numMovies,reverseStr,'inputStr','checking movies','waitbarOn',options.waitbarOn,'displayEvery',10);
		end
	end
	if options.getMovieDims==1
		outputMovie = dims;
		return;
	end
	% dims
	xDimMax = max(dims.x);
	yDimMax = max(dims.y);
	switch options.treatMoviesAsContinuous
		case 0
			zDimLength = sum(dims.z);
		case 1
			if isempty(options.frameList)
				zDimLength = sum(dims.z);
			else
				zDimLength = length(options.frameList);
			end
		otherwise
			% body
	end
	% pre-allocated output structure, convert to input movie datatype
	% if strcmp(imgClass,'single')|strcmp(imgClass,'double')
	% 	if isempty(options.loadSpecificImgClass)
	% 		outputMovie = nan([xDimMax yDimMax zDimLength],imgClass);
	% 	else
	% 		display('pre-allocating single matrix...')
	% 		outputMovie = ones([xDimMax yDimMax zDimLength],imgClass);
	% 		% j = whos('outputMovie');j.bytes=j.bytes*9.53674e-7;display(['movie size: ' num2str(j.bytes) 'Mb | ' num2str(j.size) ' | ' j.class]);
	% 		% return;
	% 		outputMovie(:,:,:) = 0;
	% 	end
	% else
	% 	outputMovie = zeros([xDimMax yDimMax zDimLength],imgClass);
	% end

	% size(outputMovie)

	if options.treatMoviesAsContinuous==1&~isempty(options.frameList)
		% totalZ = sum(dims.z);
		zdims = dims.z;
		frameList = options.frameList;
		zdimsCumsum = cumsum([0 zdims]);
		zdims = [1 zdims];
		for i=1:(length(zdims)-1)
		    g{i} = frameList>zdimsCumsum(i)&frameList<=zdimsCumsum(i+1);
		    globalFrame{i} = frameList(g{i}) - zdimsCumsum(i);
		    dims.z(i) = length(globalFrame{i});
		end
		cellfun(@max,globalFrame,'UniformOutput',false)
		cellfun(@min,globalFrame,'UniformOutput',false)
		% pause
	else
		globalFrame = [];
	end

	% reshape movie if not including larger movies
	removeIdx = [];
	numList = 1:numMovies;
	for iMovie=numList
		if isempty(globalFrame)
			thisFrameList = options.frameList;
		else
			thisFrameList = globalFrame{iMovie};
			if isempty(thisFrameList)
				subfxnDisplay(['no global frames:' num2str(iMovie) '/' num2str(numMovies) ': ' thisMoviePath],options);
				removeIdx(end+1) = iMovie;
			end
		end
	end

	subfxnDisplay('-------',options);
	keepIdx = setdiff(numList,removeIdx);
	xDimMax = max(dims.x(keepIdx));
	yDimMax = max(dims.y(keepIdx));

	% pre-allocated output structure, convert to input movie datatype
	if strcmp(imgClass,'single')|strcmp(imgClass,'double')
		if isempty(options.loadSpecificImgClass)
			outputMovie = nan([xDimMax yDimMax zDimLength],imgClass);
		else
			subfxnDisplay('pre-allocating single matrix...',options);
			outputMovie = ones([xDimMax yDimMax zDimLength],imgClass);
			% j = whos('outputMovie');j.bytes=j.bytes*9.53674e-7;display(['movie size: ' num2str(j.bytes) 'Mb | ' num2str(j.size) ' | ' j.class]);
			% return;
			outputMovie(:,:,:) = 0;
		end
	else
		outputMovie = zeros([xDimMax yDimMax zDimLength],imgClass);
	end
	subfxnDisplay('-------',options);
	% ========================
	for iMovie=1:numMovies
	    thisMoviePath = movieList{iMovie};

	    [options.movieType] = getMovieFileType(thisMoviePath);

	    if isempty(globalFrame)
	    	thisFrameList = options.frameList;
	    else
	    	thisFrameList = globalFrame{iMovie};
	    	if isempty(thisFrameList)
				subfxnDisplay(['no global frames:' num2str(iMovie) '/' num2str(numMovies) ': ' thisMoviePath],options);
	    		continue
	    	end
	    end

	    if options.displayInfo==1
	    	subfxnDisplay(['loading ' num2str(iMovie) '/' num2str(numMovies) ': ' thisMoviePath],options);
	    end
	    % depending on movie type, load differently
		switch options.movieType
			case 'tiff'
				if isempty(thisFrameList)
		    		tmpMovie = load_tif_movie(thisMoviePath,1);
		    		tmpMovie = tmpMovie.Movie;
		    	else
		    		tmpMovie = load_tif_movie(thisMoviePath,1,'Numberframe',thisFrameList);
		    		tmpMovie = tmpMovie.Movie;
		    	end
	    	% ========================
			case 'hdf5'
				if isempty(thisFrameList)
					hinfo = hdf5info(thisMoviePath);
					% hReadInfo = hinfo.GroupHierarchy.Datasets(1);
					% datasetNames = {hinfo.GroupHierarchy.Datasets.Name};
					% thisDatasetName = strmatch(inputDatasetName,datasetNames);
					% hReadInfo = hinfo.GroupHierarchy.Datasets(thisDatasetName);
					hReadInfo = getHdf5Info();
					% read in the file
	                % hReadInfo.Attributes
					tmpMovie = hdf5read(hReadInfo);
					if isempty(options.loadSpecificImgClass)
					else
						tmpMovie = cast(tmpMovie,imgClass);
					end
				else
					% xxxx = 1;
					if nanmin(diff(thisFrameList)==1)==1
						% if contiguous segment of HDF5 file, read that in as one block to save time if user specifies as such
						inputFilePath = thisMoviePath;
						framesToGrab = thisFrameList;
						hinfo = hdf5info(inputFilePath);
						% hReadInfo = hinfo.GroupHierarchy.Datasets(1);
						[hReadInfo thisDatasetName] = getHdf5Info();
						xDim = hReadInfo.Dims(1);
						yDim = hReadInfo.Dims(2);

						subfxnDisplay(['loading movie as contiguous chunk: ' num2str([0 0 framesToGrab(1)-1]) ' | ' num2str([xDim yDim length(framesToGrab)])],options);
						% tmpMovie = readHDF5Subset(inputFilePath,[0 0 framesToGrab(1)-1],[xDim yDim length(framesToGrab)],'datasetName',options.inputDatasetName);
						% size(tmpMovie)
						% size(outputMovie)
						if ischar(options.inputDatasetName)
							tmpMovie = readHDF5Subset(inputFilePath,[0 0 framesToGrab(1)-1],[xDim yDim length(framesToGrab)],'datasetName',options.inputDatasetName,'displayInfo',options.displayInfo);
						else
							tmpMovie = readHDF5Subset(inputFilePath,[0 0 framesToGrab(1)-1],[xDim yDim length(framesToGrab)],'datasetName',thisDatasetName,'displayInfo',options.displayInfo);
						end
					else
						% read frame-by-frame to save space
						inputFilePath = thisMoviePath;
						hinfo = hdf5info(inputFilePath);
						% hReadInfo = hinfo.GroupHierarchy.Datasets(1);
						[hReadInfo thisDatasetName] = getHdf5Info();
						xDim = hReadInfo.Dims(1);
						yDim = hReadInfo.Dims(2);
						% tmpMovie = readHDF5Subset(inputFilePath,[0 0 thisFrameList(1)],[xDim yDim length(thisFrameList)],'datasetName',options.inputDatasetName);
						framesToGrab = thisFrameList;
						nFrames = length(framesToGrab);
						reverseStr = '';
						for iframe = 1:nFrames
							readFrame = framesToGrab(iframe);

							if ischar(options.inputDatasetName)
								thisFrame = readHDF5Subset(inputFilePath,[0 0 readFrame-1],[xDim yDim 1],'datasetName',options.inputDatasetName,'displayInfo',0);
							else
								thisFrame = readHDF5Subset(inputFilePath,[0 0 readFrame-1],[xDim yDim 1],'datasetName',thisDatasetName,'displayInfo',0);
							end
							if isempty(options.loadSpecificImgClass)
								tmpMovie(:,:,iframe) = thisFrame;
							else
						    	% assume 3D movies with [x y frames] as dimensions
							    if(iMovie==1)
									outputMovie(1:dims.x(iMovie),1:dims.y(iMovie),iframe) = cast(thisFrame,imgClass);
							    else
							    	zOffset = sum(dims.z(1:iMovie-1));
							    	outputMovie(1:dims.x(iMovie),1:dims.y(iMovie),(zOffset+iframe)) = cast(thisFrame,imgClass);
							    end
							end
							if options.displayInfo==1
								reverseStr = cmdWaitbar(iframe,nFrames,reverseStr,'inputStr','loading hdf5','waitbarOn',options.waitbarOn,'displayEvery',50);
							end
						end
					end
				end
			% ========================
			case 'avi'
				xyloObj = VideoReader(thisMoviePath);

				if isempty(thisFrameList)
					nFrames = xyloObj.NumberOfFrames;
					framesToGrab = 1:nFrames;
				else
					nFrames = length(thisFrameList);
					framesToGrab = thisFrameList;
				end
				vidHeight = xyloObj.Height;
				vidWidth = xyloObj.Width;

				% Preallocate movie structure.
				tmpMovie = zeros(vidHeight, vidWidth, nFrames, 'uint8');

				% Read one frame at a time.
				reverseStr = '';
				iframe = 1;
				nFrames = length(framesToGrab);
				for iframe = 1:nFrames
					readFrame = framesToGrab(iframe);
					tmpAviFrame = read(xyloObj, readFrame);
					% check if frame is RGB or grayscale, if RGB only take one channel (since they will be identical for RGB grayscale)
					if size(tmpAviFrame,3)==3
						tmpAviFrame = squeeze(tmpAviFrame(:,:,1));
					end
				    tmpMovie(:,:,iframe) = tmpAviFrame;
		            % reduce waitbar access
		            if options.displayInfo==1
		    			reverseStr = cmdWaitbar(iframe,nFrames,reverseStr,'inputStr','loading avi','waitbarOn',options.waitbarOn,'displayEvery',50);
		    		end
		    		iframe = iframe + 1;
				end
			% ========================
			otherwise
				% let's just not deal with this for now
				return;
		end
		if exist('tmpMovie','var')
		    if(iMovie==1)
				outputMovie(1:dims.x(iMovie),1:dims.y(iMovie),1:dims.z(iMovie)) = tmpMovie;
		        % outputMovie(:,:,:) = tmpMovie;
		    else
		    	% assume 3D movies with [x y frames] as dimensions
		    	zOffset = sum(dims.z(1:iMovie-1));
		    	outputMovie(1:dims.x(iMovie),1:dims.y(iMovie),(zOffset+1):(zOffset+dims.z(iMovie))) = tmpMovie;
		        % outputMovie(:,:,end+1:end+size(tmpMovie,3)) = tmpMovie;
		    end
	    	clear tmpMovie;
		else

		end
	end

	% hinfo = hdf5info('A:\shared\concatenated_2013_07_05_p62_m728_MAG1.h5');
	% DFOF = hdf5read(hinfo.GroupHierarchy.Datasets(1));
	% get size of movie
	% DFOFsize = size(DFOF.Movie);
	movieDims = size(outputMovie);
	nPixels = movieDims(1)*movieDims(2);
    if length(movieDims)==2
        nFrames = 1;
    else
        nFrames = movieDims(3);
    end
	if options.waitbarOn==1
	    subfxnDisplay(['movie class: ' class(outputMovie)],options);
	    subfxnDisplay(['movie size: ' num2str(size(outputMovie))],options);
	    subfxnDisplay(['x-dims: ' num2str(dims.x)],options);
	    subfxnDisplay(['y-dims: ' num2str(dims.y)],options);
	    subfxnDisplay(['z-dims: ' num2str(dims.z)],options);
	end
	j = whos('outputMovie');j.bytes=j.bytes*9.53674e-7;
	subfxnDisplay(['movie size: ' num2str(j.bytes) 'Mb | ' num2str(j.size) ' | ' j.class],options);
    % display(dims);
	% Convert the movie to single
	% DFOF=single(DFOF);
	if options.convertToDouble==1
	    subfxnDisplay('converting to double...',options);
		outputMovie=double(outputMovie);
	end
	function [hReadInfo thisDatasetName] = getHdf5Info();
		if ischar(options.inputDatasetName)
			try
			    datasetNames = {hinfo.GroupHierarchy.Datasets.Name};
% 			    thisDatasetName = strmatch(options.inputDatasetName,datasetNames);%CUSTOM_CODE!
                thisDatasetName = datasetNames{1};%CUSTOM_CODE!
                try%CUSTOM_CODE!
                   hReadInfo = hinfo.GroupHierarchy.Datasets(datasetNames);%CUSTOM_CODE!
                catch%CUSTOM_CODE!
                    hReadInfo = hinfo.GroupHierarchy.Datasets;%CUSTOM_CODE!
                end%CUSTOM_CODE!
	        catch
	            try
	                datasetNames = {hinfo.GroupHierarchy.Groups.Datasets.Name};
	                thisDatasetName = strmatch(options.inputDatasetName,datasetNames);
	                hReadInfo = hinfo.GroupHierarchy.Groups.Datasets(thisDatasetName);
	            catch
	                nGroups = length(hinfo.GroupHierarchy.Groups);
	                datasetNames = {};
	                for groupNo = 1:nGroups
	                    datasetNames{groupNo} = hinfo.GroupHierarchy.Groups(groupNo).Datasets.Name;
	                end
	                thisDatasetName = strmatch(options.inputDatasetName,datasetNames);
	                thisGroupNo = strmatch(options.inputDatasetName,datasetNames);
	                % hReadInfo = hinfo.GroupHierarchy.Groups(thisDatasetName).Datasets;
	                % hinfo.GroupHierarchy.Groups(thisGroupNo).Datasets.Name
	                thisDatasetNo = strmatch(options.inputDatasetName,{hinfo.GroupHierarchy.Groups(thisGroupNo).Datasets.Name});
	                hReadInfo = hinfo.GroupHierarchy.Groups(thisGroupNo).Datasets(thisDatasetNo);
	            end
			end
		else
			hReadInfo = hinfo.GroupHierarchy.Datasets(options.inputDatasetName);
			thisDatasetName = hinfo.GroupHierarchy.Datasets(options.inputDatasetName).Name;
		end
	end
end

function [movieType supported] = getMovieFileType(thisMoviePath)
    % determine how to load movie, don't assume every movie in list is of the same type
	supported = 1;
    try
		[pathstr,name,ext] = fileparts(thisMoviePath);
	catch
		movieType = '';
		supported = 0;
		return;
	end
	% files are assumed to be named correctly (lying does no one any good)
	if strcmp(ext,'.h5')|strcmp(ext,'.hdf5')
		movieType = 'hdf5';
	elseif strcmp(ext,'.tif')|strcmp(ext,'.tiff')
		movieType = 'tiff';
	elseif strcmp(ext,'.avi')
		movieType = 'avi';
	else
		movieType = '';
		supported = 0;
	end
end
function subfxnDisplay(str,options)
	if options.displayInfo==1
		disp(str)
	end
end