function [inputSnr inputMse] = computeSignalSnr(inputSignals,varargin)
	% obtains an approximate SNR for an input signal
	% biafra ahanonu
	% started: 2013.11.04 [11:54:09]
	% inputs
		% inputSignals: [nSignals frame] matrix
	% outputs
		% inputSnr: [1 nSignals] vector of calculated SNR. NaN used where SNR is not calculated.
		% inputMse: [1 nSignals] vector of MSE. NaN used where MSE is not calculated.
	% options
		% % type of SNR to calculate
		% options.SNRtype = 'mean(signal)/std(noise)';
		% % frames around which to remove the signal for noise estimation
		% options.timeSeq = [-10:10];
		% % show the waitbar
		% options.waitbarOn = 1;
		% % save time if already computed peaks
		% options.testpeaks = [];
		% options.testpeaksArray = [];
	% changelog
		% 2013.12.08 now uses RMS to calculate the SNR after removing the signal to get an estimated noise trace.


	%========================
	% type of SNR to calculate
	options.SNRtype = 'mean(signal)/std(noise)';
	% frames around which to remove the signal for noise estimation
	options.timeSeq = [-10:10];
	% show the waitbar
	options.waitbarOn = 1;
	% save time if already computed peaks
	options.testpeaks = [];
	options.testpeaksArray = [];
	% get options
	options = getOptions(options,varargin);
	%========================
	try
		display(['SNR type: ' options.SNRtype])
		% to later calculate the signal idx
		% outerFun = @(x,y) x+y;
		nSignals = size(inputSignals,1);
		reverseStr = '';
		% calculate peak locations
		if isempty(options.testpeaks)
			[testpeaks testpeaksArray] = computeSignalPeaks(inputSignals,'makeSummaryPlots',0,'waitbarOn',options.waitbarOn);
		else
			testpeaks = options.testpeaks;
			testpeaksArray = options.testpeaksArray;
		end
		inputSnr = NaN([1 nSignals]);
		inputMse = NaN([1 nSignals]);
		for signalNo=1:nSignals
			loopSignal=inputSignals(signalNo,:);
			testpeaks = testpeaksArray{signalNo};
			switch options.SNRtype
				case 'mean(signal)/std(noise)'
					% Xapp=zeros(1,length(X));
					if ~isempty(testpeaks)
						% peakIdx = bsxfun(outerFun,options.timeSeq',testpeaks);
						peakIdx = bsxfun(@plus,options.timeSeq',testpeaks);
						peakIdx = unique(peakIdx(:));
						if ~isempty(peakIdx)
							% remove peaks outside range of signal
							peakIdx(find(peakIdx>length(loopSignal)))=[];
							peakIdx(find(peakIdx<=0))=[];

							% remove signal then add back in noise based on signal statistics
							noiseSignal = loopSignal;
							noiseSignal(peakIdx) = NaN;
							% noiseSignal(peakIdx) = normrnd(nanmean(noiseSignal),nanstd(noiseSignal),[1 length(peakIdx)]);

							% remove noise from signal vector
							xtmp = zeros([1 length(loopSignal)]);
							xtmp(peakIdx) = 1;
							loopSignal(~logical(xtmp)) = NaN;

							% compute SNR
							% x_snr = (rootMeanSquare(loopSignal)/rootMeanSquare(noiseSignal))^2;
							x_snr = nanmean(loopSignal)/nanstd(noiseSignal);
							xRms = rootMeanSquare(loopSignal);
						else
							x_snr = NaN;
						end
					else
						x_snr = NaN;
						xRms = NaN;
					end

					inputSnr(signalNo)=x_snr;
					inputMse(signalNo)=xRms;
					% Xapp(testpeaks)=X(testpeaks);
					% [psnr,mse,maxerr,L2rat] = measerr(X,Xapp);
					% IcaSnr(signalNo)=psnr;
					% IcaMse(signalNo)=mse;
				otherwise
					display('incorrect SNR type')
			end

			% print progress
			if options.waitbarOn==1
				reverseStr = cmdWaitbar(signalNo,nSignals,reverseStr,'inputStr','obtaining SNR','waitbarOn',options.waitbarOn,'displayEvery',50);
			end
		end
	catch err
        display(repmat('@',1,7))
        disp(getReport(err,'extended','hyperlinks','on'));
        display(repmat('@',1,7))
    end
end
function [RMS] = rootMeanSquare(x)
	RMS = sqrt(nanmean(x.^2));
end