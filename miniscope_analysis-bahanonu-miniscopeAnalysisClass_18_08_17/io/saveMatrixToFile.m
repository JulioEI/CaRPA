function [success] = saveMatrixToFile(inputMatrix,savePath,varargin)
	% save 3D matrix to arbitrary file type (HDF5, TIF, AVI for now)
	% biafra ahanonu
	% started: 2016.01.12 [11:09:53]
	% inputs
		% inputMatrix - [x y frame] movie matrix
		% savePath - character string of path to file with extension included,
	% outputs
		%

	% changelog
		%
	% TODO
		% Add checking of data size so tiff can be automatically switched

	%========================
	% default or force a save type
	options.saveType = 'avi';
	% how to save AVI, e.g. 'Motion JPEG AVI', see https://www.mathworks.com/help/matlab/ref/videowriter.html#inputarg_profile
	options.aviSaveType = 'Grayscale AVI';
	% whether to have the waitbar enabled
	options.waitbarOn = 1;
	% hierarchy name in hdf5 where movie is
	options.inputDatasetName = '/1';
	% frame rate, e.g. for AVI
	options.saveFPS = 20;
	% get options
	options = getOptions(options,varargin);
	% display(options)
	% unpack options into current workspace
	% fn=fieldnames(options);
	% for i=1:length(fn)
	% 	eval([fn{i} '=options.' fn{i} ';']);
	% end
	%========================

	try
		success = 0;
		[options.saveType supported] = getMovieFileType(savePath);
		if supported==0
			display('Unsupported save type, supported types are tif, tiff, h5, hdf5, and avi.')
			return;
		end
		startTime = tic;
		movieClass = class(inputMatrix);
		switch options.saveType
			case 'avi'
				fprintf('Saving to: %s',savePath);
				%
				nFrames = size(inputMatrix,3);
				writerObj = VideoWriter(savePath,options.aviSaveType);
				writerObj.FrameRate = options.saveFPS;
				open(writerObj);
				switch movieClass
					case 'single'
						inputMatrix = normalizeVector(inputMatrix,'normRange','zeroToOne');
					case 'double'
						inputMatrix = normalizeVector(inputMatrix,'normRange','zeroToOne');
					otherwise
						% do nothing
				end
				% same each frame to AVI file
				reverseStr = '';
				for frameNo = 1:nFrames
					thisFrame = inputMatrix(:,:,frameNo);
					writeVideo(writerObj,thisFrame);
					reverseStr = cmdWaitbar(frameNo,nFrames,reverseStr,'inputStr','saving avi','waitbarOn',options.waitbarOn,'displayEvery',50);
				end
				close(writerObj);
			case 'tiff'
				tiffOptions.comp = 'no';
				tiffOptions.overwrite = true;
				fprintf('Saving to: %s',savePath);
				saveastiff(inputMatrix, savePath, tiffOptions);
			case 'hdf5'
				fprintf('Saving to: %s',savePath);
				[output] = writeHDF5Data(inputMatrix,savePath,'datasetname',options.inputDatasetName);
			otherwise
				%
		end
		endTime = toc(startTime);
		fprintf('Done! Time elapsed: %0.1f seconds | %0.3f minutes \n',endTime,endTime/60);
		success = 1;
	catch err
		display(repmat('@',1,7))
		disp(getReport(err,'extended','hyperlinks','on'));
		display(repmat('@',1,7))
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