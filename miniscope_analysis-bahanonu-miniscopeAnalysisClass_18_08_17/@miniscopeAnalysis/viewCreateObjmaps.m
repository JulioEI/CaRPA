function obj = viewCreateObjmaps(obj)
	% creates obj maps and plots of high-SNR example signals
	% biafra ahanonu
	% branched from controllerAnalysis: 2014.08.01 [16:09:16]
	% inputs
		%
	% outputs
		%

	% changelog
		% 2017.01.14 [20:06:04] - support switched from [nSignals x y] to [x y nSignals]
	% TODO
		%

	[fileIdxArray idNumIdxArray nFilesToAnalyze nFiles] = obj.getAnalysisSubsetsToAnalyze();

	if obj.guiEnabled==1
		movieSettings = inputdlg({...
				'directory to save pictures: ',...
				'video file filter',...
				'video frames',...
				'image threshold (0-1)'
			},...
			'view movie settings',1,...
			{...
				obj.picsSavePath,...
				obj.fileFilterRegexp,...
				'1:500',...
				'0.5'
			}...
		);
		obj.picsSavePath = movieSettings{1};
		obj.fileFilterRegexp = movieSettings{2};
		userVideoFrames = str2num(movieSettings{3});
		userThreshold = str2num(movieSettings{4});
	else
		% obj.picsSavePath
		% obj.fileFilterRegexp
		userThreshold = 0.5;
		userVideoFrames = 1:500;
	end

	[figHandle figNo] = openFigure(959, '');
	[figHandle figNo] = openFigure(969, '');
	[figHandle figNo] = openFigure(970, '');
	[figHandle figNo] = openFigure(971, '');
	[figHandle figNo] = openFigure(972, '');

	for thisFileNumIdx = 1:nFilesToAnalyze
		try
			thisFileNum = fileIdxArray(thisFileNumIdx);
			obj.fileNum = thisFileNum;
			display(repmat('=',1,21))
			display([num2str(thisFileNum) '/' num2str(nFiles) ': ' obj.fileIDNameArray{obj.fileNum}]);
			% =====================
			% for backwards compatibility, will be removed in the future.
			nIDs = length(obj.stimulusNameArray);
			%
			nameArray = obj.stimulusNameArray;
			idArray = obj.stimulusIdArray;
			%
			% [inputSignals inputImages signalPeaks signalPeakIdx] = modelGetSignalsImages(obj,'returnType','raw');
			[inputSignals inputImages signalPeaks signalPeakIdx] = modelGetSignalsImages(obj,'returnType','filtered');
			if isempty(inputSignals);display('no input signals');continue;end
			% size(signalPeakIdx)
			% return
			% [inputSignals inputImages signalPeaks signalPeakIdx] = modelGetSignalsImages(obj);
			nIDs = length(obj.stimulusNameArray);
			nSignals = size(inputSignals,1);
			nFrames = size(inputSignals,2);
			%
			options.dfofAnalysis = obj.dfofAnalysis;
			timeSeq = obj.timeSequence;
			% subject = obj.subjectNum{obj.fileNum};
			subject = obj.subjectStr{obj.fileNum};
			assay = obj.assay{obj.fileNum};
			%
			framesPerSecond = obj.FRAMES_PER_SECOND;
			subjAssayIDStr = obj.fileIDNameArray{obj.fileNum};
			%
			figNoAll = obj.figNoAll;
			figNo = obj.figNo;
			figNames = obj.figNames;
			% magic numbers!
			% amount of time to make object maps before/after a stimulus
			prepostTime = 20;
			%
			picsSavePath = [obj.picsSavePath filesep 'cellmaps' filesep];
			fileFilterRegexp = obj.fileFilterRegexp;
			% =====================


	    	% thisFileID = obj.fileIDNameArray{obj.fileNum};
	    	thisFileID = obj.fileIDArray{obj.fileNum};

	    	normalFigs = 1;
	    	if normalFigs==1

	    		% [figHandle figNo] = openFigure(85, '');
	    		% 	[inputSignals inputImages signalPeaks signalPeakIdx] = modelGetSignalsImages(obj,'returnType','raw');

	    		% 	[signalSnr a] = computeSignalSnr(inputSignals,'testpeaks',signalPeaks,'testpeaksArray',signalPeakIdx);
	    		% 	[signalSnr sortedIdx] = sort(signalSnr,'descend');

	    		% 	inputImages2 = inputImages;
	    		% 	for cellNo = 1:size(inputImages,3)
	    		% 		inputImages2(:,:,cellNo) = normalizeVector(inputImages2(:,:,cellNo),'normRange','zeroToOne');
	    		% 	end
	    		% 	emImages2 = zeros([size(inputImages2,1) size(inputImages2,2) 3])
	    		% 	for i=1:3
	    		% 		randColorVector = 1*rand([size(inputImages2,3) 1]);
	    		% 		randColorVector = randColorVector.*matchAcross;
	    		% 		% randColorVector = randColorVector+0.05;
	    		% 		randColorVector(randColorVector==0) = 0.2;
	    		% 		% size(randColorVector)
	    		% 		% figure;plot(randColorVector)
	    		% 		emImages2(:,:,i) = nanmax(groupImagesByColor(inputImages2,randColorVector,'thresholdImages',0),[],3);
	    		% 	end

	    		% 	clear inputImages2
	    		% 	imagesc(emImages2);drawnow
	    		% 	axis off

	    		% 	markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};

	    		% 	legend()

	    		% 	title([subject ' | ' assay ' | overlap map | ' num2str(size(signalPeaks,1)) ' cells'],'fontsize',20)
	    		% 	set(figHandle,'PaperUnits','inches','PaperPosition',[0 0 9 9])
	    		% 	% figure(figHandle)
	    		% 	obj.modelSaveImgToFile([],'cellmapObjColor_','current',[]);

    			[figHandle figNo] = openFigure(85, '');
    				inputImages2 = inputImages;
    				for cellNo = 1:size(inputImages,3)
    					inputImages2(:,:,cellNo) = normalizeVector(inputImages2(:,:,cellNo),'normRange','zeroToOne');
    				end
    				emImages2 = zeros([size(inputImages2,1) size(inputImages2,2) 3]);
    				for i=1:3
    					emImages2(:,:,i) = nanmax(groupImagesByColor(inputImages2,1*rand([size(inputImages2,3) 1]),'thresholdImages',0),[],3);
    				end
    				clear inputImages2
    				imagesc(emImages2);drawnow
    				axis off
    				title([subject ' | ' assay ' | overlap map | ' num2str(size(signalPeaks,1)) ' cells'],'fontsize',20)
    				set(figHandle,'PaperUnits','inches','PaperPosition',[0 0 9 9])
    				% figure(figHandle)
    				obj.modelSaveImgToFile([],'cellmapObjColor_','current',[]);



	    			markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'}

	    		[figHandle figNo] = openFigure(7989, '');
	    			thres = thresholdImages(inputImages,'binary',1,'threshold',userThreshold);
	    			thisCellmap = createObjMap(thres,'mapType','sum');
	    			imagesc(thisCellmap);colorbar;
	    			% colormap(obj.colormap);
	    			colormap([0 0 0;jet(nanmax(thisCellmap(:)))])
	    			title([subject ' | ' assay ' | overlap map | ' num2str(size(signalPeaks,1)) ' cells'],'fontsize',20)
	    			set(figHandle,'PaperUnits','inches','PaperPosition',[0 0 9 9])
	    			% figure(figHandle)
	    			obj.modelSaveImgToFile([],'cellmapObjOverlap_','current',[]);
			    [figHandle figNo] = openFigure(969, '');
				    s1 = subplot(1,2,1);
					    % coloredObjs = groupImagesByColor(thresholdImages(inputImages),[]);
					    % thisCellmap = createObjMap(coloredObjs);
					    % firing rate grouped images
					    display(['signalPeaks: ' num2str(size(signalPeaks))])
					    numPeakEvents = sum(signalPeaks,2);
					    numPeakEvents = numPeakEvents/size(signalPeaks,2)*framesPerSecond;
					    display(['inputImages: ' num2str(size(inputImages))])
					    display(['numPeakEvents: ' num2str(size(numPeakEvents))])
					    thres = thresholdImages(inputImages,'threshold',userThreshold);
					    [groupedImagesRates] = groupImagesByColor(inputImages,numPeakEvents);
					    thisCellmap = createObjMap(groupedImagesRates);

					    % if fileNum==1
					    %     fig1 = figure(32);
					    %     % colormap gray;
					    % end
						% thisCellmap = createObjMap([thisDirSaveStr options.rawICfiltersSaveStr]);
						% subplot(round(nFiles/4),4,fileNum);
						plotBinaryCellMapFigure();
						title([subject ' | ' assay ' | firing rate map | ' num2str(size(signalPeaks,1)) ' cells'],'fontsize',20)
						hold on;

					[signalSnr a] = computeSignalSnr(inputSignals,'testpeaks',signalPeaks,'testpeaksArray',signalPeakIdx);
				[figHandle figNo] = openFigure(972, '');
					clf
					subplot(2,3,1);
						plotBinaryCellMapFigure();
				[figHandle figNo] = openFigure(969, '');
					s2 = subplot(1,2,2);
						if nSignals>1
							[signalSnr sortedIdx] = sort(signalSnr,'descend');
                            sortedinputSignals = signalPeaks.*inputSignals;
                            % [signalSnr sortedIdx] = sort(sum(sortedinputSignals,2),'descend');
                            [signalSnr sortedIdx] = sort(max(sortedinputSignals,[],2),'descend');
							sortedinputSignals = inputSignals(sortedIdx,:);
							display('==============')
							display(['signalPeakIdx: ' num2str(size(signalPeakIdx))])
							display(['sortedIdx: ' num2str(size(sortedIdx))])
							display('==============')
							signalPeakIdx = {signalPeakIdx{sortedIdx}};
							cutLength = 200;
							if cutLength*2>size(inputSignals,2);cutLength=floor(size(inputSignals,2)/2.2);end
							cutLength
							nSignalsShow = 20;
							if nSignalsShow>length(signalPeakIdx);nSignalsShow=length(signalPeakIdx);end
							sortedinputSignalsCut = zeros([nSignalsShow cutLength*2+1]);
							% display(['sortedinputSignalsCut: ' num2str(size(sortedinputSignalsCut))])
							shiftVector = round(linspace(round(cutLength/10),round(cutLength*0.9),nSignalsShow));
							shiftVector = shiftVector(randperm(length(shiftVector)));
							for i=1:nSignalsShow
								spikeIdx = signalPeakIdx{i};
								spikeIdxValues = sortedinputSignals(i,spikeIdx);
								[k tmpIdx] = max(spikeIdxValues);
								if isempty(tmpIdx)
									continue;
								end
								spikeIdx = spikeIdx(tmpIdx(1));
                                spikeIdx = spikeIdx(:);
								spikeIdx = spikeIdx-(round(cutLength/2)-shiftVector(i));
								% spikeIdx
								% cutLength
								nPoints = size(inputSignals,2);
								try
									if (spikeIdx-cutLength)<0
										beginDiff = abs(spikeIdx-cutLength);
										cutIdx = bsxfun(@plus,spikeIdx,-(cutLength-beginDiff-1):(cutLength+beginDiff+1));
										cutIdx = 1:(cutLength*2+1);
									elseif (spikeIdx+cutLength)>nPoints
										endDiff = abs(-spikeIdx);
										cutIdx = bsxfun(@plus,spikeIdx,-(cutLength+endDiff+1):(cutLength-endDiff-1));
										cutIdx = (nPoints-(cutLength*2)):nPoints;
                                    else
										cutIdx = bsxfun(@plus,spikeIdx,-cutLength:cutLength);
									end
                                catch err
                                    display(repmat('@',1,7))
                                    disp(getReport(err,'extended','hyperlinks','on'));
                                    display(repmat('@',1,7))
									cutIdx = [];
                                    cutIdx = -cutLength:cutLength;
								end
								if ~isempty(cutIdx)
									sortedinputSignalsCut(i,:) = sortedinputSignals(i,cutIdx(:)');
								end
                            end
                            %size(sortedinputSignalsCut)
                            %imagesc(sortedinputSignals)
                            % sortedinputSignalsCut = sortedinputSignals(1:nSignalsShow,:);
							% sortedinputSignalsCut = flip(sortedinputSignalsCut,1);
							display(['sortedinputSignalsCut: ' num2str(size(sortedinputSignalsCut))])
							% sortedinputSignalsCut = sortedinputSignals(1:7,:);
							% sortedinputSignalsCut = inputSignals(1:7,:);
							plotTracesFigure();
						else
							plot(inputSignals);
							xlabel('frames','fontsize',20);ylabel('\Delta F/F','fontsize',20);
							box off;
							title('example traces','fontsize',20);
						end

					d=0.02; %distance between images
					set(s1,'position',[d 0.15 0.5-2*d 0.8])
			     	set(s2,'position',[0.5+d 0.15 0.5-2*d 0.8])
				    saveFile = char(strrep(strcat(picsSavePath,'cellmap_',thisFileID,''),'/',''));
				    set(figHandle,'PaperUnits','inches','PaperPosition',[0 0 16 9])
				    % figure(figHandle)
				    obj.modelSaveImgToFile([],'cellmapObj_','current',[]);
				    % print('-dpng','-r200',saveFile)
				    % print('-dmeta','-r200',saveFile)
				    % saveas(gcf,saveFile);
					drawnow;

					[figHandle figNo] = openFigure(972, '');
						subplot(2,3,2);
							if nSignals>1
								plotTracesFigure();
							else
								plot(inputSignals);
								xlabel('frames','fontsize',20);ylabel('\Delta F/F','fontsize',20);
								box off;
								% title('example traces','fontsize',20);
							end

				[figHandle figNo] = openFigure(970, '');
					% timeVector = (1:size(sortedinputSignalsCut,2))/framesPerSecond;
					if nSignals>1
						plotSignalsGraph(sortedinputSignalsCut,'LineWidth',2.5);
						nTicks = 10;
						set(gca,'XTick',round(linspace(1,size(sortedinputSignalsCut,2),nTicks)))
						labelVector = round(linspace(1,size(sortedinputSignalsCut,2)/framesPerSecond,nTicks));
						set(gca,'XTickLabel',labelVector);
						xlabel('seconds','fontsize',20);ylabel('\Delta F/F','fontsize',20);
						box off;
						% axis off;
						% title('example traces');
					else
						plot(inputSignals);
						xlabel('frames','fontsize',20);ylabel('\Delta F/F','fontsize',20);
						box off;
						title('example traces','fontsize',20);
					end
					title([subject ' | ' assay ' | example traces'],'fontsize',20)
				    saveFile = char(strrep(strcat(picsSavePath,'traces_',thisFileID,''),'/',''));
				    saveFile
				    set(figHandle,'PaperUnits','inches','PaperPosition',[0 0 16 9])
				    % figure(figHandle)
				    obj.modelSaveImgToFile([],'cellmapTraces_','current',[]);
				    % print('-dpng','-r200',saveFile)
				    % print('-dmeta','-r200',saveFile)
				    % saveas(gcf,saveFile);
					drawnow;
			end
			movieList = getFileList(obj.inputFolders{obj.fileNum}, 'concat');
			if isempty(movieList)
				movieList = getFileList(obj.inputFolders{obj.fileNum}, fileFilterRegexp);
			end
			if ~isempty(movieList)
				movieFrame = loadMovieList(movieList{1},'convertToDouble',0,'frameList',1:100,'inputDatasetName',obj.inputDatasetName);
				% movieFrame = squeeze(movieFrame(:,:,1));
				movieFrame = squeeze(max(movieFrame,[],3));
				% imagesc(imadjust(movieFrame));
				% imagesc(movieFrame);
				% imshow(movieFrame);
				% axis off; colormap gray;
				% title([subject ' | ' assay ' | blue>green>red percentile rank']);
				% hold on;
				% imcontrast
				% continue
				% inputImagesThresholded = thresholdImages(inputImages,'binary',0)/3;
				% inputImagesThresholded = inputImages;
				% icaQ = quantile(numPeakEvents,[0.3 0.6]);
				% colorObjMaps{1} = createObjMap(inputImagesThresholded(numPeakEvents<icaQ(1),:,:));
				% colorObjMaps{2} = createObjMap(inputImagesThresholded(numPeakEvents>icaQ(1)&numPeakEvents<icaQ(2),:,:));
				% colorObjMaps{3} = createObjMap(inputImagesThresholded(numPeakEvents>icaQ(2),:,:));

				[inputSignals inputImages signalPeaks signalPeakIdx] = modelGetSignalsImages(obj,'returnType','raw');
				nSawSignals = size(inputSignals,1);

				try obj.valid{obj.fileNum}.(obj.signalExtractionMethod).manual;check.manual=1; catch check.manual=0; end
				if check.manual==1
					display(['using valid.' obj.signalExtractionMethod '.manual identifications...']);
					validAuto = obj.valid{obj.fileNum}.(obj.signalExtractionMethod).manual;
				else
					validAuto = obj.valid{obj.fileNum}.(obj.signalExtractionMethod).auto;
				end
				validAuto = logical(validAuto);
				% inputImages = inputImages(:,:,validAuto);
				% validAuto = logical(obj.validManual{obj.fileNum});
				display('==============')
				if isempty(obj.validRegionMod)
					% validRegionMod = ones(size(validAuto));
					validRegionMod = [];
				else
					validRegionMod = obj.validRegionMod{obj.fileNum};
					% validRegionMod = validRegionMod(logical(validAuto));
					% validRegionMod = validRegionMod(validAuto);
				end
				validRegionMod = logical(validRegionMod);
				% class(validAuto)
				% validRegionMod
				display(['validRegionMod: ' num2str(size(validRegionMod))])
				display(['validAuto: ' num2str(size(validAuto))])
				% inputImagesThresholded = thresholdImages(inputImages(validAuto,:,:),'binary',0)/3;
				inputImagesThresholded = thresholdImages(inputImages,'binary',0,'threshold',userThreshold)/3;
				% inputImagesThresholded = inputImagesThresholded(validAuto);
				display(['inputImagesThresholded: ' num2str(size(inputImagesThresholded))])
				if isempty(validRegionMod)
					colorObjMaps{1} = createObjMap(inputImagesThresholded(:,:,validAuto==0));
					colorObjMaps{2} = createObjMap(inputImagesThresholded(:,:,validAuto==1));
				else
					colorObjMaps{1} = createObjMap(inputImagesThresholded(:,:,validRegionMod==0));
					colorObjMaps{2} = createObjMap(inputImagesThresholded(:,:,validRegionMod==1));
				end
				display(['colorObjMaps{1}: ' num2str(size(colorObjMaps{1}))])
				display(['colorObjMaps{2}: ' num2str(size(colorObjMaps{2}))])
				E = normalizeVector(double(movieFrame),'normRange','zeroToOne')/2;
				% E = E';
				display(['E: ' num2str(size(E))])
				if isempty(colorObjMaps{1})
					Comb(:,:,1) = E;
				else
					Comb(:,:,1) = E+normalizeVector(double(colorObjMaps{1}),'normRange','zeroToOne')/4; % red
				end
				Comb(:,:,2) = E+normalizeVector(double(colorObjMaps{2}),'normRange','zeroToOne')/4; % green
				Comb(:,:,3) = E; % blue
				% Comb(:,:,3) = E+normalizeVector(double(colorObjMaps{1}),'normRange','zeroToOne')/4; % blue
				[figHandle figNo] = openFigure(971, '');
					imagesc(Comb)
					% clear Comb
					axis off; colormap gray;
					title([subject ' | ' assay ' | blue-green-red percentile rank | cells=' num2str(nSignals)]);

					[nanmax(movieFrame(:)) nanmin(movieFrame(:))]
					[nanmax(colorObjMaps{1}(:)) nanmin(colorObjMaps{1}(:))]

					% zeroMap = zeros(size(movieFrame));
					% oneMap = ones(size(movieFrame));
					% green = cat(3, zeroMap, oneMap, zeroMap);
					% blue = cat(3, zeroMap, zeroMap, oneMap);
					% red = cat(3, oneMap, zeroMap, zeroMap);
					% warning off
					% blueOverlay = imshow(blue);
					% greenOverlay = imshow(green);
					% redOverlay = imshow(red);
					% warning on
					% set(redOverlay, 'AlphaData', colorObjMaps{1});
					% set(greenOverlay, 'AlphaData', colorObjMaps{2});
					% set(blueOverlay, 'AlphaData', colorObjMaps{3});
					set(gca, 'LooseInset', get(gca,'TightInset'))
					hold off;
					saveFile = char(strrep(strcat(picsSavePath,'cellmap_overlay_',thisFileID,''),'/',''));
					saveFile
					set(figHandle,'PaperUnits','inches','PaperPosition',[0 0 16 9])
					% figure(figHandle)
					obj.modelSaveImgToFile([],'cellmapObjOverlay_','current',[]);
					% print('-dpng','-r200',saveFile)
					% print('-dmeta','-r200',saveFile)
					% saveas(gcf,saveFile);
					% pause

				[figHandle figNo] = openFigure(972, '');
					subplot(2,3,3);
						imagesc(movieFrame)
						title('raw movie')

				[figHandle figNo] = openFigure(972, '');
					subplot(2,3,4);
						imagesc(Comb)
						clear Comb
						axis off; colormap gray;
						title('detected cells in green')
						% title([subject ' | ' assay ' | blue-green-red percentile rank | cells=' num2str(nSignals)]);
						% num2str(size(signalPeaks,1))
				movieList2 = getFileList(obj.inputFolders{obj.fileNum}, fileFilterRegexp);
				if ~isempty(movieList2)
                    movieDims = loadMovieList(movieList2{1},'getMovieDims',1,'inputDatasetName',obj.inputDatasetName);
                    if movieDims.z<max(userVideoFrames);nFramesMax=1:movieDims.z;else nFramesMax=userVideoFrames;end
					movieFrame = loadMovieList(movieList2{1},'convertToDouble',0,'frameList',nFramesMax,'inputDatasetName',obj.inputDatasetName);
					movieFrame = squeeze(max(movieFrame,[],3));
					E = normalizeVector(double(movieFrame),'normRange','zeroToOne')/2;
					display(['E: ' num2str(size(E))])
					if isempty(colorObjMaps{1})
						Comb(:,:,1) = E;
					else
						Comb(:,:,1) = E;
						% Comb(:,:,1) = E+normalizeVector(double(colorObjMaps{1}),'normRange','zeroToOne')/4; % red
					end
					Comb(:,:,2) = E+normalizeVector(double(colorObjMaps{2}),'normRange','zeroToOne')/4; % green
					Comb(:,:,3) = E; % blue

					[figHandle figNo] = openFigure(972, '');
						subplot(2,3,5);
							imagesc(movieFrame)
							title(sprintf('processed movie (%d:%d frames max)',min(userVideoFrames),max(userVideoFrames)))

					[figHandle figNo] = openFigure(972, '');
						subplot(2,3,6);
							imagesc(Comb)
							clear Comb
							axis off; colormap gray;
							title('detected cells in green')
							% title([subject ' | ' assay ' | blue-green-red percentile rank | cells=' num2str(nSignals)]);
							% num2str(size(signalPeaks,1))
				end

				suptitle([subject ' | ' assay ' | firing rate map | ' num2str(sum(validAuto)) ' (' num2str(nSawSignals) ') cells'])
				% titleAxes = axes('Position', [0, 0.95, 1, 0.05],'Visible','off');
				% axes(titleAxes)
				% cla(ax,'reset')
				% set(titleAxes, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' );
				% set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' );
				% text(0.5, 0, [subject ' | ' assay ' | firing rate map | ' num2str(sum(validAuto)) ' cells'], 'FontSize', 14', 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' );
                set(figHandle,'PaperUnits','inches','PaperPosition',[0 0 22 16])
				obj.modelSaveImgToFile([],'objMapAll_','current',[]);

				% inputImagesThresholded = thresholdImages(inputImages,'binary',0);
				% saveFile = char(strcat(thisDirSaveStr,'cellmap_thresholded.h5'));
				% thisObjMap = createObjMap(inputImagesThresholded);
				% movieSaved = writeHDF5Data(thisObjMap,saveFile)
				% inputImagesThresholded = thresholdImages(inputImages,'binary',1);
				% saveFile = char(strcat(thisDirSaveStr,'cellmap_thresholded_binary.h5'));
				% thisObjMap = createObjMap(inputImagesThresholded);
				% movieSaved = writeHDF5Data(thisObjMap,saveFile)
			end
		catch err
			display(repmat('@',1,7))
			disp(getReport(err,'extended','hyperlinks','on'));
			display(repmat('@',1,7))
		end
	end
	function plotBinaryCellMapFigure()
		imagesc(thisCellmap);
		colormap(obj.colormap);
		s2Pos = get(gca,'position');
		cb = colorbar('location','southoutside'); ylabel(cb, 'Hz');
		set(gca,'position',s2Pos);
	    % colormap hot; colorbar;
		% title(regexp(thisDir,'m\d+', 'match'));
		box off;
		axis tight;
		axis off;
		set(gca, 'LooseInset', get(gca,'TightInset'))

	end
	function plotTracesFigure()
		if size(sortedinputSignalsCut,1)==1
			plot(sortedinputSignalsCut)
        else
			plotSignalsGraph(sortedinputSignalsCut,'LineWidth',2.5);
		end
		nTicks = 10;
		set(gca,'XTick',round(linspace(1,size(sortedinputSignalsCut,2),nTicks)))
		labelVector = round(linspace(1,size(sortedinputSignalsCut,2),nTicks)/framesPerSecond);
		set(gca,'XTickLabel',labelVector);
		xlabel('seconds','fontsize',20);ylabel('\Delta F/F','fontsize',20);
		box off;
		% title('example traces','fontsize',20);
	end
end