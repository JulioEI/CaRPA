function [frameVector frameVectorDecision frameVectorReset] = getSubsampleFrameMatrix(nFrames,nIterations,percentFramesPerIteration,varargin)
	% creates a vector of either randomly sampled or semi-randomly sampled frames
	% biafra ahanonu
	% started: 2017.04.01
	% inputs
		% nFrames - integer, number of frames to run
		% nIterations - integer, number of iterations
		% percentFramesPerIteration - between 0 and 1, fraction of frames to subsample per iteration
	% outputs
		%

	% changelog
		%
	% TODO
		%

	%========================
	% random, resampleRemaining
	options.subsampleMethod = 'random';
	%
	options.displayPlots = 0;
	% 0 to 1, for subsampleMethod=resampleRemaining, percent of remaining to subsample
	options.percentRemainingSubsample = 0.5;
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
		frameVector = zeros([nIterations nFrames]);
		frameVectorReset = zeros([nIterations nFrames]);
		frameVectorDecision = zeros([nIterations nFrames]);
		nSubsampled = round(nFrames*percentFramesPerIteration);
		% z = binocdf(1:nIterations,nFrames,percentFramesPerIteration);
		for iterNo = 1:nIterations
			switch options.subsampleMethod
				case 'random'
					selectedFrames = randsample(nFrames, nSubsampled);
					frameVector(iterNo,selectedFrames) = 1;
		    		frameVectorDecision(iterNo,selectedFrames) = 1;
					if iterNo~=1
						 frameVector(iterNo,:) = frameVector(iterNo,:)|frameVector(iterNo-1,:);
					end
				case 'resampleRemaining'
					if iterNo~=1
						fracFramesPrev = sum(frameVectorReset(iterNo-1,:))/nFrames;
					else
						fracFramesPrev = 0;
					end
					% fracFramesPrev
					if iterNo==1|fracFramesPrev==1
						selectedFrames = randsample(nFrames, nSubsampled);
					else
						unchosenFrames = find(frameVectorReset(iterNo-1,:)==0);
						% fprintf('value: %d\n',length(unchosenFrames))
						% unchosenFrames = randsample(unchosenFrames, round(length(unchosenFrames)/3));
						if length(unchosenFrames)>nSubsampled
							remainFramesToSubsample = round(nSubsampled*options.percentRemainingSubsample);
							unchosenFrames = randsample(unchosenFrames, remainFramesToSubsample);
						end
						selectedFrames = randsample(setdiff(1:nFrames,unchosenFrames), nSubsampled-length(unchosenFrames));
						selectedFrames = [selectedFrames unchosenFrames];
						% selectedFrames =	randsample(unchosenFrames, nSubsampled);
						fprintf('value: %d\n',length(selectedFrames))
					end
					frameVector(iterNo,selectedFrames) = 1;
		    		frameVectorDecision(iterNo,selectedFrames) = 1;
		    		frameVectorReset(iterNo,selectedFrames) = 1;
					if iterNo~=1
			    		if fracFramesPrev~=1
			    			frameVectorReset(iterNo,:) = frameVectorReset(iterNo,:)|frameVectorReset(iterNo-1,:);
			    		end
						 frameVector(iterNo,:) = frameVector(iterNo,:)|frameVector(iterNo-1,:);
					end
				otherwise
					% do nothing
			end

		end

		if options.displayPlots==1
			figure
				imagesc(frameVectorDecision)
				title(sprintf('fraction is %d',percentFramesPerIteration*100))
			figure(2)
				% subplot(1,2,memorySwitch+1)
				%plot(sum(frameVector,2)/nFrames)
				switch options.subsampleMethod
					case 'random'
						semilogx(1:nIterations,sum(frameVector,2)/nFrames)

					case 'resampleRemaining'
						semilogx(1:nIterations,sum(frameVector,2)/nFrames,'-.')
					otherwise
						% nothing
				end
				hold on
				drawnow
		end
	catch err
		display(repmat('@',1,7))
		disp(getReport(err,'extended','hyperlinks','on'));
		display(repmat('@',1,7))
	end
end