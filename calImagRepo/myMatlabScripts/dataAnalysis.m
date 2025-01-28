classdef dataAnalysis < handle
    %This is a parent class giving basic utility to general analysis 
    % functions and general analysis parameters. To be inherited by more 
    % specific classes.
    
    properties
        %
        animalList = [];
        
        %Analysis parameters
        dtCamera = 0.05; %period of miniscope (seconds/frame)
        dT = 0.05; %Integration parameter (seconds)
        binN = [20,1]; %Number of bins in X Y
        minVel = 4;%[4, 0]; %Frames under this speed will be discarded (cm/s)
        pxPerCm = [6.5, 6.5]; %How many px a centimeter is (px/cm)
        traceField = 'rawProb'; %Default field to extract traces
        spikeField = 'spikeDeconv'; %Default field to extract spikes
        filterPF = false; %Remove cells with non significant place fields
    end
    
    properties (Constant)
        %File regexp
        tracesEventsRegexp = {'TracesAndEvents'};
        analysisFileRegexp = {'emAnalysis_?\d*.mat$','pcaicaAnalysis_?\d*.mat$'};
        decisionsRegexp = {'Sorted_?\d*.mat$','decisions_?\d*.mat$'}

        %Time parsers
        dateParser = 'yyyymmdd'; %Reads and prints time with this format
        timeParser = 'HHMMSS';
    end
     
    methods (Static)
        
        function [rTB,xTB] = integrateInTime(r,x,dT,dtCamera)
            %Takes a response and position matrix sampled at dtCamera and
            %integrates in time with a moving window of size dT. Sums
            %responses and averages positions.
            
            timeBin = dT/dtCamera;
            sizeTB = floor(length(x)/timeBin);
            rTB = zeros([sizeTB,size(r,2)]);
            xTB = zeros([sizeTB,2]);
            for k = 1:sizeTB
                lb = round((k-1)*timeBin + 1);
                ub = round((k)*timeBin);
                if ~isempty(r)
                    rTB(k,:) = sum(r(lb:ub,:),1);
                else
                    rTB = r;
                end
                xTB(k,:) = mean(x(lb:ub,:),1);
            end            
        end
        
        function [posB,binCounts] = binPosition(x,binN)
            %Takes a matrix [time * 2] of positions and returns the binned 
            %matrix from 1 to binN.
            
            %Normalize position
            posN = (x - min(x)) ./ (max(x) - min(x));
            
            %Find bin edges
            [binCounts,binEdgesX,binEdgesY] = histcounts2(posN(:,1),posN(:,2),binN);
            
            %Bin positions
            posB = zeros(size(posN));
            for dim = 1:2
                if dim == 1;dimBinEdges = binEdgesX;else;dimBinEdges = binEdgesY;end
                for k = 1:(length(dimBinEdges)-1)
                    posInBin = find(posN(:,dim)>=dimBinEdges(k) & posN(:,dim)<=dimBinEdges(k+1));
                    posB(posInBin,dim) = k;
                end
            end
        end
        
        function matHot = ind2OneHot(idx,totalSize)
            %Converts cells of indices into one hot encoding matrix
            if iscell(idx)
                matHot = zeros([totalSize,length(idx)]);
                for cell = 1:length(idx)
                    matHot(idx{cell},cell) = 1;
                end
            else
                matHot = idx;
            end
        end

        function [xTB,rTB,binCounts,cmPerBin] = curateXandR(x,r,param)
            
            if isfield(param,'v')
                v = param.v;
            else
                %v = abs(diff(x(:,1))/param.dtCamera);%When v not specified we only take the x velocity
                v = abs(diff(x(:,1)./param.pxPerCm(1))/param.dtCamera);%When v not specified we only take the x velocity

            end
            %Looks for outliers, removes slow periods, integrates in time and bins x
            
            %Get parameters
            dT = param.dT; %Temporal bin size (s)
            binN = param.binN; %Number of bins
            dtCamera = param.dtCamera; %Period of camera (s)
            pxPerCm = param.pxPerCm; %Pixels/cm ratio
            minVel = param.minVel; %Minimal speed (cm/s)
            
            %Check for outliers
            outliers = x>(mean(x)+3*std(x));
            if sum(outliers(:,1)) ~= 0 %only for x
                warning('OUTLIERS DETECTED, MANUAL ACTION REQUIERED')
                %pause
                %x = x(100:end,:);
                %r = r(100:end,:);
            end
            
            %Put pos in cm units and remove slow frames
            x = x./pxPerCm(1:sum(size(x)>1));
            [x,r] = dataAnalysis.removeFramesUnderMinVel(x,r,v,minVel);
                       
            %Integrate in time
            [rTB,xTB] = dataAnalysis.integrateInTime(r,x,dT,dtCamera);
%             figure;plot(x,'-*');hold on;plot(linspace(5,length(x)+5,length(xTB)),xTB,'o');title('>4 cm/s');plot(repmat(1:length(x),[20,1])',(ones([20,length(x)]).*linspace(10,110,20)')','k--')
%             figure;plot(linspace(5,length(x)+5,length(xTB)),dataAnalysis.binPosition(xTB,binN),'-o');hold on;plot(dataAnalysis.binPosition(x,binN),'-*');plot(repmat(1:length(x),[20,1])',(ones([20,length(x)]).*linspace(1,20,20)')','k--')
            
            %Bin positions
            [xTB,binCounts] = dataAnalysis.binPosition(xTB,binN);
            
            %Get cm per bin
            cmPerBin = range(x)./binN;
            
            %Remove non needed dimensions
            xTB = xTB(:,1:size(x,2));
            
        end
        
        function [x,r] = removeFramesUnderMinVel(x,r,v,minVel)
            %Get absolute velocity
            %v = abs(diff(x)/dtCamera);
            
            %Get frames under min velocity
            framesUnderMinVel = find(v < minVel);%unique([find(v(:,1) < minVel(1)); find(v(:,2) < minVel(2))]);
            
            %Remove frames under min velocity
            x(framesUnderMinVel,:) = [];
            if ~isempty(r)
                r(framesUnderMinVel,:) = [];
            end
        end
        
        function [xR,xL,rR,rL] = separateLR1D(r,x)
            %Separates data into left right frames
            %Get only directional PF (linear track only)
            framesL = find(diff(x(:,1))<0);
            xL = x(framesL,:);
            rL = r(framesL,:);
            framesR = find(diff(x(:,1))>0);
            xR = x(framesR,:);
            rR = r(framesR,:);
        end
        
    end
    
    methods 
        function param = getAnalysisParameters(obj)
            param = obj.getDefaultParams({'binN','minVel','pxPerCm','dtCamera','dT','filterPF'},{});
        end
    end
    
    methods (Access = protected)
                
        function param = getDefaultParams(obj,myFields,input)
            for field = myFields
                param.(field{1}) = dataAnalysis.parseInput({input,field{1},obj.(field{1})});
            end
        end
        
    end
    
    methods(Static, Access = protected) 
         
        function match = checkRegexp(sourceText,regexpList)
            match = zeros(size(sourceText));
            k = 1;
            for text = sourceText
                regexpMatch = 0;
                for regexpSeq = regexpList
                    if regexp(text{1},regexpSeq{1})
                        regexpMatch = regexpMatch + 1;
                    end
                end
                if regexpMatch > 0
                    match(k) = 1;
                end
                k = k + 1;
            end
        end  
        
        function myOutput = selectOutput(output,message)
            if nargin < 2
                message = 'Select output:';
            end
            myFields = fields(output);
            s = listdlg('PromptString',message,...
                'SelectionMode','single',...
                'ListString',myFields);
            myOutput = myFields{s};
        end
        
          function myOutput = selectOutputMulti(output,message)
            if nargin < 2
                message = 'Select output:';
            end
            myFields = fields(output);
            s = listdlg('PromptString',message,...
                'SelectionMode','multiple',...
                'ListString',myFields);
            myOutput = myFields(s);
          end
         
        function [varargout] = parseInput(varargin)
            %{{inputs},optional(arg1 default1, arg2 default2, ...)}
            tmpVar = varargin{1};
            if iscell(tmpVar{1})
                argumentList = tmpVar(2:2:end);
                defaultVal = tmpVar(3:2:end);
                input = tmpVar{1};
                inputArg = input(1:2:end);
                inputVal = input(2:2:end);
                for k = 1:length(argumentList)
                    idx = strcmp(argumentList{k},inputArg);
                    warning('off')
                    if sum(idx) == 0
                        varargout{k} = defaultVal{k}();
                    else
                        varargout{k} = inputVal{idx};
                    end
                    warning('on')
                end                
            else
                for k = 2:2:length(tmpVar)
                    varargout{k} = tmpVar{k};
                end
            end
        end
        
     end
end

