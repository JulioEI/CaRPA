function [ behavFix,caFix,behavIdx,caIdx, mode] = fixBehaviorCalciumLength(behav,ca,mode)
    
    param.dtCamera = 0.05; %period of miniscope (seconds/frame)
    param.dT = 0.5; %Integration parameter (seconds)
    param.binN = [20,1]; %Number of bins in X Y
    param.minVel = [4, 0]; %Frames under this speed will be discarded (cm/s)
    param.pxPerCm = [6.5, 6.5]; %How many px a centimeter is (px/cm)

    if nargin < 3
        mode = 'interpol'%'best';
    end
    behavLen = size(behav,1);
    caLen = size(ca,1);
    if behavLen > caLen
        if behavLen - caLen < 4 %Less than 4 frames difference
            behavFix = behav(1:caLen,:);
            behavIdx = 1:caLen;
            caFix = ca;
            caIdx = 1:caLen;
        else
            [caFix,behavFix,caIdx,behavIdx] = longerLenFix(ca,behav,'calcium','behavior',mode,param);
        end
    elseif caLen > behavLen
        if caLen - behavLen < 4 %Less than 4 frames difference
            caFix = ca(1:behavLen,:);
            caIdx = 1:behavLen;
            behavFix = behav;
            behavIdx = 1:behavLen;
        else
            [behavFix,caFix,behavIdx,caIdx] = longerLenFix(behav,ca,'behavior','calcium',mode,param);
        end
    end
end

function [sVecFix,lVecFix,sVecIdx,lVecIdx] = longerLenFix(sVec,lVec,sVecStr,lVecStr,mode,param)
    sVecLen = size(sVec,1);
    lVecLen = size(lVec,1);
    s0 = ['Warning, ',sVecStr,' is ',num2str(sVecLen),' frames, but ',lVecStr,' is ',num2str(lVecLen),' frames.'];
    s1 = ['Interpolate ', sVecStr,' uniformly'];
    s2 = ['Shift ',sVecStr,' in block to minimize decoder error and shorten ', lVecStr];
    
    switch mode
        case 'manual'
            choice = questdlg(s0, 'MANUAL ACTION NEEDED',s1,s2,'Abort','Abort');
        case 'best'
            choice = s2;
        case 'interpol'
            choice = s1;
        case 'shift'
            choice = s2;
    end

    switch choice
        case s1
            sVecFix = fixVecInterpolate(sVec,lVecLen);
            sVecIdx = 1:lVecLen;
            lVecFix = lVec; 
            lVecIdx = 1:lVecLen;
        case s2 
            score = zeros([1,(lVecLen-sVecLen)]);
            for k = 1:(lVecLen-sVecLen)
                lVecCut = lVec(k:(end-(lVecLen-sVecLen-k+1)),:);
                if strcmp(lVecStr,'calcium')
                    x = sVec(:,1);
                    r = lVecCut;
                    %score(k) = mean(quickForestDecoder(lVecCut,sVec(:,1),10)); %we only predict the x component
                elseif strcmp(lVecStr,'behavior')
                    x = lVecCut(:,1);
                    r = sVec;
                    %score(k) = mean(quickForestDecoder(sVec,lVecCut(:,1),10)); %we only predict the x component
                end
                [x,r] = dataAnalysis.curateXandR(x,r,param); %Bin, remove slow frames
                score(k) = mean(decoderAnalysis.decodePosKFoldParallel(x,r,10,'SVM',1,0,0));
            end
            
            bestK = find(score==min(score),1);
            
            lVecFix = lVec(bestK:(sVecLen+bestK-1),:);
            lVecIdx = bestK:(sVecLen+bestK-1);
            sVecFix = sVec;
            sVecIdx = 1:sVecLen;
            
        case 'Abort'
            error('User aborted operation')
    end
    
    if strcmp(mode,'best') %Try the best of the two approaches
        sVecInterp = fixVecInterpolate(sVec,lVecLen);
        lVecInterp = lVec;
        if strcmp(lVecStr,'calcium')
            x = sVecInterp(:,1);
            r = lVecInterp;
            %scoreInterp = mean(quickBayesDecoder(lVecInterp,sVecInterp(:,1),10));
        elseif strcmp(lVecStr,'behavior')
            x = lVecInterp(:,1);
            r = sVecInterp;
            %scoreInterp = mean(quickBayesDecoder(sVecInterp,lVecInterp(:,1),10));
        end
        
        [x,r] = dataAnalysis.curateXandR(x,r,param); %Bin, remove slow frames
        scoreInterp = mean(decoderAnalysis.decodePosKFoldParallel(x,r,10,'SVM',1,0,0));
  
        if score(bestK) > scoreInterp
            sVecFix = sVecInterp;
            sVecIdx = 1:lVecLen;
            lVecFix = lVecInterp;       
            lVecIdx = 1:lVecLen;
        end
    end
    
end

