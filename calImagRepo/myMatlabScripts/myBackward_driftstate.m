function [n, Ffit, F0, LL, xest, yfit] = myBackward_driftstate(F,par) %%

%% THis is mine (including two more outputa)
% function [n, Ffit, F0, LL, xest, yfit, cestn, xsaturation] = myBackward_driftstate(F,par)

DEBUG = eval('false');

% Input
if isempty(F), [n Ffit LL xest yfit] = deal([],[],0,[],[]); return, end
if ~isvector(F), error argument, end

% Physiological parameters
if length(par.F0)~=2, error 'when estimating a drift, F0 should be an interval', end
switch par.drift.effect
    case 'additive'
        F0 = 1;
        if ~(F0>par.F0(1) && F0<par.F0(2))
            error 'additive drifts: calcium signals must be already normalized by F0'
        else
            %disp 'additive drifts: calcium signals supposed already normalized by F0'
        end
    case 'multiplicative'
        F0 = mean(par.F0);
end
baselineinterval = par.F0 / F0;
a = par.a;
decay = exp(-par.dt/par.tau);
sat = par.saturation;
pnonlin = par.pnonlin;
if ~isempty(pnonlin) && sat~=0, error 'saturation and nonlinearity cannot be applied simultaneously', end
hill = par.hill;
spikerate = par.finetune.spikerate;
if isempty(par.finetune.sigma), error programming, end
sigmay = par.finetune.sigma/F0;
sigmab = par.drift.parameter/F0 * sqrt(par.dt);
if sigmab==0, error programming, end

% Algo parameters
estimate = par.algo.estimate;
doMAP = strcmpi(estimate,'MAP');
doproba = strcmp(estimate,'proba');
dosample = ismember(estimate,{'sample' 'samples'});
interpmode = fn_switch(doMAP,par.algo.interpmode,'linear'); % spline interpolation can yield negative weights, which is not acceptable for probabilities
nsample = fn_switch(dosample,par.algo.nsample,1);

% GPU implementation?
dogpu = par.algo.dogpu;
% (function for gpuArray<->array conversion only if requested)
if dogpu
    mygpu = @gpuArray;
    mygather = @gather;
else
    mygpu = @(x)x;
    mygather = @(x)x;
end

% Special: non-integer spikes
nonintegerspike = par.special.nonintegerspike_minamp;
if nonintegerspike && ~doMAP, error 'noninteger spike are available only for MAP estimations', end

% Get the normalized fluorescence
y = F/F0;
y = double(y(:));

% Algorithm parameters
nc = par.algo.nc; if isempty(nc), nc = par.algo.nc_norise; end % discretization
cmax = par.algo.cmax;
dc = cmax/(nc-1);
cc = (0:nc-1)'*dc; % column vector
nb = par.algo.nb; if isempty(nb), nb = par.algo.nb_driftstate; end
db = diff(baselineinterval)/(nb-1);
bb = linspace(baselineinterval(1),baselineinterval(2),nb); % row vector

% Sizes
T = length(y);

% Precomputations for the interpolation function x <- x*decay + n
nspikemax = par.algo.nspikemax;
if nonintegerspike==0
    MM = cell(1,1+nspikemax);
    for i=0:nspikemax
        % (value before to look at if there was 0, 1, 2, 3 spikes)
        if i==0
            cci = cc*decay;
        else
            cci = min(cci + 1,cmax);
        end
        % (interpolation matrixc)
        Mi = interp1(cc,eye(nc),cci,interpmode);
        MM{1+i} = Mi;
    end
else
    if ~doMAP, error 'nonintegerspike handled only for MAP estimates', end
    M0 = interp1(cc,eye(nc),cc*decay,interpmode);
    minjump = ceil(nonintegerspike/dc); % minimal calcium jump of an event
end

% Precomputations for the spike likelihood
% p(n) = exp(-rate*dt) (rate*dt)^n/n!, then take the negative log
if spikerate
    if par.special.burstcostsone
        nspikcost = [0 ones(1,nspikemax)];
    else
        nspikcost = 0:nspikemax;
    end
    lspike = +spikerate*par.dt +log(factorial(nspikcost)) -nspikcost*log(spikerate*par.dt);
    pspike = exp(-lspike)/sum(exp(-lspike)); % make the sum 1
    lspike = -log(pspike);
else
    % no a priori on spikes!
    lspike = zeros(1,1+nspikemax);
end

% Precomputation for probability update
% f1(c) = sum_n p(n) f(c*decay+n)
if doMAP && ~nonintegerspike
    MM = cat(1,MM{:});
elseif ~doMAP
    MS = 0;
    for i=0:nspikemax, MS = MS + pspike(1+i)*MM{1+i}; end
end

% Precomputation for the baseline drift
if doMAP
    % time update will involve finding the baseline drift that maximizes
    % probability; the maximum on the discretization grid will first be
    % located, then interpolation will be used to find the maximum with a
    % finer resolution
    % (drifting matrix)
    maxdrift = max(2,ceil(3*sigmab/db));
    DD = zeros(nb,2*maxdrift+1);
    for i=1:2*maxdrift+1
        DD(:,i) = max(min((1:nb)+(i-1-maxdrift),nb),1);
    end
    % ldrift = -log(1/(sqrt(2*pi)*sigmab)) + ((-maxdrift:maxdrift)*db).^2/(2*sigmab^2);
    ldrift = ((-maxdrift:maxdrift)*db).^2/(2*sigmab^2);
    ldrift = repmat(shiftdim(ldrift,-1),[nc nb]);
    % (quadratic interpolation of a triplet of points)
    tmp = eye(3);
    QQ = zeros(3);
    for i=1:3, QQ(:,i) = polyfit([-1 0 1],tmp(i,:),2); end
    QQ = QQ'; % operation on columns
else
    % time update will involve averaging probabilities accross possible
    % baseline drifts, which are described by a continuous (rather than
    % discrete) probability
    % we can construct a matrix multiplication that will realize the
    % interpolation and averaging at once; this matrix is obtained by first
    % replacing the continuous distribution by a fine-grain discrete
    % distribution
    
    % define bins: central bins have all equal probabilities, while some
    % extreme bins with lower probabilities are added
    sides = [-Inf -100 -50 -30 -20 -10 -7 -5 -3 norminv(.04:.04:.96) 3 5 7 10 20 30 50 100 Inf];
    ndrift = length(sides)-1;
    ndhalf = (ndrift-1)/2;
    pdrift = diff(normcdf(sides)); % probability of each bin
    pdrift(ndhalf+2:end) = fliplr(pdrift(1:ndhalf));  % correct numerical error for the last one
    discretesteps = norminv(cumsum(pdrift)-pdrift/2); % representative element of each bin has central probability
    discretesteps(ndhalf+2:end) = -fliplr(discretesteps(1:ndhalf));   % correct numerical error for the last one   
    pdrift = pdrift/sum(pdrift);
        
    % matrix for baseline time update
    bb1 = fn_add((1:nb)',discretesteps*(sigmab/db));
    BB = interp1(eye(nb),bb1(:),'linear',NaN); % (nb*ndrift)*nb
    BB = reshape(BB,[nb ndrift nb]);
    pdriftc = repmat(pdrift,[nb 1 nb]);
    pdriftc(isnan(BB)) = 0; pdriftc = fn_div(pdriftc,sum(pdriftc,2));
    BB(isnan(BB)) = 0;
    BB = squeeze(sum(BB.*pdriftc,2)); % nb*nb  
    BB = BB'; % will operate on columns
    
    % for the forward step only
    ldrift = -log(pdrift);
    ldrift([1 end]) = 100^2/2+log(50);  % an approximative value, instead of Inf
    ldrift([2 end-1]) = 50^2/2+log(50); % an approximative value, instead of Inf
end

% Precomputation for the measure after saturation/nonlinearity
ccn = (par.c0+cc).^hill-par.c0^hill;
if isempty(pnonlin)
    % p(y|x) = 1/(sqrt(2*pi)*sigma) exp(-(y-a*x/(1+sat*x))^2/2*sigma^2)
    dye = 1 + a * ccn./(1+sat*ccn);
else
    dye = 1 + a * polyval([fliplr(pnonlin) 1-sum(pnonlin) 0],ccn);
end
switch par.drift.effect
    case 'additive'
        xxmeasure = fn_add(dye-1,bb);
    case 'multiplicative'
        xxmeasure = fn_mult(dye,bb);
    otherwise
        error flag
end
lmeasure = -log(1/(sqrt(2*pi)*sigmay));

% Precomputation for the a priori probability of calcium c(1)
% m = spikerate*par.dt/(1-decay);
% v = spikerate*par.dt/(1-decay^2);
% pcalcium = 1/(sqrt(2*pi*v))*exp(-(cc-m).^2/(2*v));
% pcalcium = pcalcium / sum(pcalcium); % re-normalize
% lcalcium = -log(pcalcium);
% lcalcium = repmat(lcalcium,[1 nb]);
lcalcium = zeros(nc,nb);

% Debug display
if DEBUG
    tt = (0:T-1)*par.dt;
    figure(429), clf
    hda=subplot(321); 
    plot(tt,y*F0,'parent',hda)
    hx = line(0,0,'linestyle','none','marker','*','color','k','parent',hda);
    hdb=subplot(322);
    hdc=subplot(323); hdd=subplot(324); 
    hde=subplot(325); 
    ime = imagesc(cc,bb*F0,zeros(nb,nc),'parent',hde,[0 1e-3]);
    xlabel(hde,'calcium'), ylabel(hde,'baseline'), set(hde,'ydir','normal')
    hdf=subplot(326);
end

% Backward sweep
% L(x,t) remembers what is the best log-likelihood with xt=x
% L(x,t) = min_{x(t+1),..,x(T)} -log(p(x(t+1),..,x(T),y(t),..,y(T)|x(t)=x))
% while N(x,t+1) and D(x,t+1) remember respectively the number of spikes
% between t and t+1 and the baseline drift that give this best likelihood
% N(c,b,t) = argmin_n(t+1) min_{x(t+2),..,x(T)}        -log(p(n(t+1),x(t+2),..,x(T),y(t+1),..,y(T)|c(t)=c,b(t+1)=b))
% D(c,b,t) = argmin_b(t+1) min_{n(t+1),x(t+2),..,x(T)} -log(p(x(t+1),x(t+2),..,x(T),y(t+1),..,y(T)|c(t)=c,b(t)=b))
if ~doMAP
    try
        L = zeros(nc,nb,T); 
    catch
        L = zeros(nc,nb,T,'single'); 
    end
end
if doMAP
    D = zeros(nc,nb,T,'single');
    if nonintegerspike==0
        N = zeros(nc,nb,T,'uint8');
    else
        N = zeros(nc,nb,T,'single');
    end
end
for t=T:-1:1
    % L(x,t) = min_n(t+1) -log(p(n(t+1)) + min_b(t+1) -log(p(b(t+1)|b(t))) + L(x(t+1),t+1)   <- time update (minimize first over the drift in b, then over the number of spikes)
    %          - log(p(y(t)|x(t))                                                            <- measure update
    
    if DEBUG, set(hx,'xdata',tt(t),'ydata',y(t)*F0), end
    
    % Time update (find the best n(t+1))
    if t==T
        % initialization with 'empty probability' p([])
        lt = zeros(nc,nb);
    else
        % calcium time update
        if doMAP && ~nonintegerspike
            % what is the best number of spikes
            %             lt1 = lspike_interp + reshape(MM*lt,nc,nspikemax+1,nb);
            lt1 = fn_add(lspike, reshape(MM*lt,nc,nspikemax+1,nb));
            [lt , n1] = min(lt1,[],2);
            lt = squeeze(lt); n1 = squeeze(n1-1);
            N(:,:,t+1) = n1;
        elseif doMAP && nonintegerspike
            % decay
            lt1_noevent = lspike(1)+M0*lt; % lspike(1) is the cost of zero spike
            % jump
            lt1_event = lt1_noevent;
            jumps = zeros(nc,nb);
            % (initialize with hypothetical jumps (nc-minevent)->nc+1 of
            % infinite cost)
            ltk = Inf;
            jumpk = minjump+1;
            for k=nc:-1:(1+minjump)
                ltk1 = lt(k,:);
                smaller = (ltk1<=ltk); % does a jump of amplitude minjump reach a value which is more interesting than the current minimum?
                ltk(smaller) = ltk1(smaller);
                jumpk(smaller) = minjump;
                jumpk(~smaller) = jumpk(~smaller)+1;
                lt1_event(k-minjump,:) = ltk;
                jumps(k-minjump,:) = jumpk;
            end
            lt1_event = lspike(2)+lt1_event; % lspike(2) is the cost of 1 spike (here, 1 'event') 
            % what is better between the decay and the optimal jump
            lt = lt1_noevent;
            smaller = (lt1_event<lt1_noevent);
            lt(smaller) = lt1_event(smaller);
            N = reshape(N,[nc*nb T]);
            N(smaller,t) = jumps(smaller)*dc;
            N = reshape(N,[nc nb T]);
        elseif ~doMAP
            lt = logmultexp(MS,lt); % interpolation of probabilities rather than of log-probabilities
        end
        
        if DEBUG
            set(ime,'cdata',log2proba(lt)')
            set(hde,'ydir','normal')
            if doMAP
                imagesc(cc,bb,N(:,:,t+1)','parent',hdf,[0 3])
                set(hdf,'ydir','normal')
            end
            drawnow
        end
        
        % baseline time update
        if doMAP
            % what is the optimal baseline drift
            % get: lt = min_b(t+1) -log(p(b(t+1)|b(t))) + L(x(t+1),t+1)
            lt1 = ldrift + reshape(lt(:,DD),[nc nb 2*maxdrift+1]);
            [lt idrift] = min(lt1,[],3);
            
            % find a over-sampling minimum using a quadratic interpolation when
            % drifting values are not on the sides defined by the maximum
            % allowed
            oksides = ~(idrift==1 | idrift==2*maxdrift+1);
            lt1 = reshape(lt1,[nc*nb 2*maxdrift+1]);
            lt1ok = lt1(oksides,:);
            idriftok = idrift(oksides);
            nok = sum(oksides(:));
            indices3 = fn_add((1:nok)'+nok*(idriftok-1),nok*[-1 0 1]);
            values3 = lt1ok(indices3);
            qq = values3 * QQ;
            idriftmin = -qq(:,2)./(2*qq(:,1)); % (q(x) = ax^2 + bx + c -> the min is -b/2a)
            idrift(oksides) = idrift(oksides) + idriftmin;
            lt(oksides) = (qq(:,1).*idriftmin + qq(:,2)).*idriftmin + qq(:,3);
            
            D(:,:,t+1) = (idrift-1-maxdrift)*db;
            
            if DEBUG
                imagesc(cc,bb,log2proba(lt)','parent',hdc,[0 1e-3])
                imagesc(cc,bb,D(:,:,t+1)','parent',hdd,[-1 1]*maxdrift*db)
                set([hdc hdd],'ydir','normal')
            end
        else
            lt = logmultexp_column(BB,lt); % interpolation of probabilities
        end
    end
    
    % Measure update
    lt = lt + (lmeasure+(y(t)-xxmeasure).^2/(2*sigmay^2));
    if ~doMAP, L(:,:,t) = lt; end
    
    % A priori on calcium concentration at t=1
    if t==1
        lt = lt + lcalcium;
    end
    
    if DEBUG && doMAP
        imagesc(cc,bb,log2proba(lt)','parent',hdb,[0 1e-3])
        set(hdb,'ydir','normal')
        drawnow
    end
end

% Precomputations for forward sweep
if doproba
    % Precomputations for the interpolation function f1(c) = f((c-n)/decay)
    % and for probability update
    MS = 0; NS = 0;
    for i=0:nspikemax
        % (value before to look at if there was 0, 1, 2, 3 spikes)
        cci = (cc-i)/decay;
        % (interpolation matrices: interpolated probabilities will be zero where cc<nspike)
        Mi = interp1(cc,eye(nc),cci,interpmode,0);
        if i==0
            Mi(1,1:2) = [decay 1-decay]; % part of the probability in the 'zero calcium' bin at time t should be interpolated from the first 'non-zero calcium' bin at time t-1
        end
        % f1(c) = sum_n p(n) f((c-n)/decay)
        MS = MS + pspike(1+i)*Mi;
        % f1(c) = sum_n n p(n) f((c-n)/decay)
        NS = NS + i*pspike(1+i)*Mi;
    end
elseif dosample
    lspike_drift = fn_add(column(lspike),row(ldrift));
end

% Forward collecting/sampling/smoothing step
if doproba
    n = zeros(T,nsample,'single');
else
    n = zeros(T,nsample,'uint8');
end
doxest = ~doproba || (nargout>=2) || par.dographsummary;
if doxest, xest = zeros(T,2,nsample,'single'); end
if dosample && strcmp(par.display,'steps'), fn_progress('sampling',T), end
for t=1:T
    if dosample && strcmp(par.display,'steps'), fn_progress(t), end
    if t==1
        if doMAP
            if par.drift.baselinestart
                % impose that the initial calcium level is baseline
                cidx = 1;
                ystart = mean(y(1:ceil(0.1/par.dt))); % average over 100ms to get the start value
                [dum bidx] = min(abs(ystart-xxmeasure(cidx,:))); 
                LL = lt(cidx,bidx);
            else
                % LL is the minimum negative log likelihood
                [LL cidx] = min(lt,[],1);
                [LL bidx] = min(LL,[],2);
                cidx = cidx(bidx);
            end
            xest(t,:) = [cc(cidx) bb(bidx)];
        elseif dosample
            % initiate samples
            LL = []; % would be quite useless to compute log likelihoods, isn't it?
            [cidx bidx] = logsample(lt,nsample);
            xest(t,1,:) = cc(cidx);
            xest(t,2,:) = bb(bidx);
        elseif doproba
            LL = logsumexp(lt(:));
            pt = log2proba(lt);
            if doxest
                xest(t,1) = sum(row(fn_mult(cc,pt)));
                xest(t,2) = sum(row(fn_mult(bb,pt)));
            end
            
            % now lt represents -log p(xt|y1,..,yt), so we use the time
            % zero prior and perform a single measure update
            lt = lcalcium + (lmeasure+(y(t)-xxmeasure).^2/(2*sigmay^2));
        end
    else
        if doMAP
            xest(t,2) = fn_coerce(xest(t-1,2) + D(cidx,bidx,t),baselineinterval);
            bidx = 1+round((xest(t,2)-bb(1))/db);
            n(t) = N(cidx,bidx,t);
            xest(t,1) = min(xest(t-1,1)*decay + double(n(t)),cmax);
            cidx = 1+round(xest(t,1)/dc);
        elseif dosample
            % draw calcium and baseline evolutions at once
            % too difficult this time to do all particles at once 
            % -> use a for loop
            nspike = column(0:nspikemax);                        % putative number of spikes
            ct = fn_add(xest(t-1,1,:)*decay, nspike);           % corresponding putative calcium values [(1+nspikemax)*1*nsample]
            ct1 = repmat(ct,[1 ndrift 1]);                      % idem [(1+nspikemax)*ndrift*nsample]
            bt = fn_add(xest(t-1,2,:), discretesteps*sigmab);   % putative baseline values [1*ndrift*nsample]
            bt1 = repmat(bt,[1+nspikemax 1 1]);                  % idem [(1+nspikemax)*ndrift*nsample]
            Lt = L(:,:,t);                                      % -log p(yt,..,yT|xt) [nc*nb]
            Lt = mygpu(Lt); % the computation in the next line seems to be the only one where GPU is profitable!
            ltk0 = interpn(Lt,1+ct1/dc,1+(bt1-bb(1))/db,'linear',Inf); % -log p(yt,..,yT|xt)   [(1+nspikemax)*ndrift*nsample]
            ltk0 = mygather(ltk0);
            ltk = bsxfun(@plus,lspike_drift,ltk0);                   % ~ -log p(xt|x(t-1),yt,..,yT) [(1+nspikemax)*ndrift*nsample]
            [cidx bidx] = logsample(ltk,'2D');                  % [nsample]
            n(t,:) = cidx-1;
            xest(t,1,:) = ct(sub2ind([1+nspikemax nsample],cidx,1:nsample));
            xest(t,2,:) = bt(sub2ind([ndrift nsample],bidx,1:nsample));
            badsample = all(all(isinf(ltk))); % some samples ran out uncharted low-proba territory: put them back in the max-proba position
            if any(badsample)
                xest(t,1,badsample) = 0;
                [~, bidx] = min(Lt(1,:));
                xest(t,2,badsample) = bb(bidx);
            end
            %             for ksample = 1:nsample
            %                 ct = xest(t-1,1,ksample)*decay + nspike;    % corresponding putative calcium values
            %                 bt = xest(t-1,2,ksample) + discretesteps*sigmab;        % putative baseline values
            %                 ltk0 = interpn(Lt,1+ct/dc,1+(bt-bb(1))/db,'linear',Inf); % -log p(yt,..,yT|ct,Bt) [size (1+nspikemax)*nb]
            %                 ltk = lspike_drift + ltk0;                   % ~ -log p(xt|x(t-1),yt,..,yT) [size (1+nspikemax)*ndrift]
            %                 [cidx bidx] = logsample(ltk);
            %                 n(t,ksample) = cidx-1;
            %                 xest(t,1,ksample) = ct(cidx);
            %                 xest(t,2,ksample) = bt(bidx);
            %             end
        elseif doproba
            if eval('true')
                % Below is the implementation described in the paper.
                % But this does not seem stable enough, in particular,
                % because it combines the past and future of each t, 
                % the result can be inconsistant (due to numerical
                % approximations probably) between time t-1 and time t.
                % Furthermore, the initialization for this computation is
                % incorrect, as it should only take y1 into account.
                
                % time update
                % . lt represents information from the past only, it will be
                %   updated from -log p(x(t-1)|y1,..,y(t-1)) to -log p(xt|y1,..,yt)
                % . L(:,:,t) is -log p(yt,..,yT|xt), i.e. combines information
                %   from past and future
                % note that interpolations must occur in the 'proba' rather
                % than 'log proba' to be accurate
                
                % lt time update
                lt1 = lt;               % -log p(x(t-1)|y1,..,y(t-1))
                lmin = min(lt1(:));
                pt1 = exp(lmin-lt1);    % ~ p(x(t-1)|y1,..,y(t-1))
                pt1b = pt1*BB;
                pt = MS*pt1b;% ~ p(xt|y1,..,y(t-1))
                lt = lmin-log(pt);      % -log p(xt|y1,..,y(t-1))
                
                % update L(:,:,t), i.e. combine lt and previous L(:,:,t)
                lty = lt + L(:,:,t);    % ~ -log p(xt|y)
                L(:,:,t) = lty;
                pty = log2proba(lty);   % ~ p(xt|y)
                
                % expectancy for number of spikes
                nt = (NS*pt1b)./pt; nt(pt==0) = 0;      % E(nt|xt,y1,..,y(t-1))
                n(t) = sum(row(nt.*pty));               % E(nt|y)
                if doxest
                    xest(t,1) = sum(row(fn_mult(cc,pty)));  % E(ct|y)
                    xest(t,2) = sum(row(fn_mult(bb,pty)));  % E(bt|y)
                end
                
                % lt measure update
                lt = lt + (lmeasure+(y(t)-xxmeasure).^2/(2*sigmay^2));
            else
                
                % Well... an alternative that would produce results more
                % similar to the 'samples' mode is too difficult to write,
                % in particular because it might involve square matrices
                % with side the total number of states
                
                error 'not implemented'
                %                 lt1 = lt;               % -log p(x(t-1)|y1,..,y(t-1))
                %                 lmin = min(lt1(:));
                %                 pt1 = exp(lmin-lt1);    % p(x(t-1)|y1,..,y(t-1))
                %
                %                 ltfuture = L(:,:,t);    % -log p(xt|yt,..,yT)
                %                 ptfuture = exp(min(ltfuture(:))-ltfuture);
                %                 ptfuture = ptfuture/sum(ptfuture(:));   % p(xt|yt,..,yT)
                
                
            end
        end
    end
end

% We can stop here if we want only spikes
if ~doxest, return, end

% Saturation/nonlinearity and scaling
cestn = squeeze(par.c0+xest(:,1,:)).^hill-par.c0^hill;
if isempty(pnonlin)
    xsaturation = a * cestn./(1+sat*cestn);
else
    xsaturation = a * polyval([fliplr(pnonlin) 1-sum(pnonlin) 0],cestn);
end

% Predicted measure (taking drifts into account)
switch par.drift.effect
    case 'additive'
        yfit = xsaturation + squeeze(xest(:,2,:));
    case 'multiplicative'
        yfit = (1+xsaturation).*squeeze(xest(:,2,:));
    otherwise
        error flag
end

% Back from normalized to data scale
Ffit = yfit*F0;

% Reajust F0 and xest(:,2) to make the mean of xest(:,2) 1
avgb = mean(row(xest(:,2,:)));
F0 = F0 * avgb;
xest(:,2,:) = xest(:,2,:) / avgb;

%-------------------------------------------------------------------------%