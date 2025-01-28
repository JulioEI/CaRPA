function [success speedup poolSizeVector] = parallelWorkerSpeedTest(numHands,numPlayers,varargin)
	% test speed-up on specific hardware
	% based on https://www.mathworks.com/help/distcomp/examples/simple-benchmarking-of-parfor-using-blackjack.html
	% biafra ahanonu
	% started: 2016.10.31
	% inputs
		% numHands = number of hands
		% numPlayers = number of players, most important for parallelization
	% outputs
		% success - 1=run no errors, 0 = run with errors
		% speedup - speedup for each # of workers, [1 nWorkers]
		% poolSizeVector - number of workers used in each run, [1 nWorkers]

	% changelog
		%
	% TODO
		%

	%========================
	% options.exampleOption = '';
	% execute parallel workers
	% options.parallel = 1;

	% options.numHands = 1e4;
	% options.numPlayers = 100; % most important for parallelization

	% get options
	% options = getOptions(options,varargin);
	% display(options)
	% unpack options into current workspace
	% fn=fieldnames(options);
	% for i=1:length(fn)
	% 	eval([fn{i} '=options.' fn{i} ';']);
	% end
	%========================

	try
		% Check number of inputs.
		if nargin == 0
			display('loading default arguments')
			numHands = 1e4;
			numPlayers = 100;
		end

		% Open parallel pool, can also run manageParallelWorkers.m
			% check maximum number of cores available
			numWorkersToOpen = java.lang.Runtime.getRuntime().availableProcessors;

			% check that local matlabpool configuration is correct
			myCluster = parcluster('local');
			if myCluster.NumWorkers~=numWorkersToOpen
				myCluster.NumWorkers = numWorkersToOpen; % 'Modified' property now TRUE
				saveProfile(myCluster);   % 'local' profile now updated
			end

			% open max workers if pool not already available
			if ~isempty(gcp('nocreate'))
			else
				% matlabpool('open',maxCores-1);
				parpool('local',numWorkersToOpen,'IdleTimeout', Inf);
			end

			success = 0;
			p = gcp;
			if isempty(p)
			    error('pctexample:backslashbench:poolClosed', ...
			        ['This example requires a parallel pool. ' ...
			         'Manually start a pool using the parpool command or set ' ...
			         'your parallel preferences to automatically start a pool.']);
			end
			poolSize = p.NumWorkers;

		fprintf('Simulating each player playing %d hands.\n', numHands);
		t1 = zeros(1, poolSize);
		poolSizeVector = 1:(poolSize+1);
		for n = poolSizeVector
		    tic;
		        % pctdemo_aux_parforbenchLocal(numHands, n*numPlayers, n);
		        pctdemo_aux_parforbenchLocal(numHands, numPlayers, n-1);
		    t1(n) = toc;
		    fprintf('%d workers simulated %d players in %3.2f seconds.\n', ...
		            n-1, numPlayers, t1(n));
		end
		t1
		% speedup = (1:poolSize).*t1(1)./t1;
		speedup = t1(2)./t1;
		fig = pctdemo_setup_blackjack(1.0);
		fig.Visible = 'on';
		ax = axes('parent', fig);
		poolSizeVector = poolSizeVector-1;
		x = plot(ax, poolSizeVector, poolSizeVector, '--', ...
		    poolSizeVector, speedup, 's', 'MarkerFaceColor', 'b');
		t = ax.XTick;
		t(t ~= round(t)) = []; % Remove all non-integer x-axis ticks.
		ax.XTick = t;
		legend(x, 'Linear Speedup', 'Measured Speedup', 'Location', 'NorthWest');
		xlabel(ax, ['Number of MATLAB workers participating in computations' 10 '(0 = parfor OFF)']);
		ylabel(ax, 'Speedup');
		success = 1;
	catch err
		success = 0;
		display(repmat('@',1,7))
		disp(getReport(err,'extended','hyperlinks','on'));
		display(repmat('@',1,7))
	end
end

function S = pctdemo_aux_parforbenchLocal(numHands, numPlayers, n)
	% PCTDEMO_AUX_PARFORBENCH Use parfor to play blackjack.
	%   S = pctdemo_aux_parforbench(numHands, numPlayers, n) plays
	%   numHands hands of blackjack numPlayers times, and uses no
	%   more than n MATLAB(R) workers for the computations.

	%   Copyright 2007-2009 The MathWorks, Inc.

	S = zeros(numHands, numPlayers);
	parfor (i = 1:numPlayers, n)
		S(:, i) = pctdemo_task_blackjack(numHands, 1);
	end
end