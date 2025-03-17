% Junsu Jang (junsu.jang94@gmail.com)
% 2025/01/20
function bMLAllMethods = estimateBetaAllMethods(Z,r,f,refFreqIdx,...
                                                betas,minNumF,isTones,...
                                                methodREs,methodIdxs,MM,y2Hyps,noiseVar)
% In this function, the loglikelihood for each betas per method are computed.
% First, the likelihood for each beta is obtained, and then in the last
% for-loop, the maximum likelihood estimates are obtained per method.
% 
% Since the nonlinear transformation, i.e., the interpolation, took the
% longest time but was commonly needed for multiple methods, we perform the
% interpolation first. 

% see ../RangeComputation/estimateRangeAllMethods.m for further comments

%% Parameters
numMethods = numel(methodREs);
numBetas = numel(betas);

% Identify the indicies of the neighboring broadband signals
toneIdx = find(isTones);
tonalNeighborIdx = [toneIdx-1; toneIdx + 1];
numTones = sum(isTones);
numComb = sum(1:(numTones-1));

% Initialize variables
logLikelihoodAll = nan(numBetas,size(Z,1),6);
bMLAllMethods = nan(methodIdxs(end),1);
xcorrCoeff = nan(numBetas,numComb);

striationIdxStartAll = nan(numBetas,1);
I = abs(Z).^2;

for bIdx=1:numBetas
    curB = betas(bIdx);
    [IInterp,striationIdx] = computeNLTData(f,curB,r,I,refFreqIdx);

    validData = ~isnan(IInterp);
    numValidTonesPerStriation = sum(validData,2);
    striationIdxStartAll(bIdx) = striationIdx;

    validStriations = numValidTonesPerStriation >= minNumF;
    IValid = IInterp(validStriations,:);

    avgMValid = mean(sqrt(IValid));
    [maxM,~] = max(avgMValid(~isTones));
    MRatio = avgMValid/(maxM);
    IRatio = MRatio.^2;
    IRatio(isTones) = 1;
    
    for m=1:numMethods
        curMethodIdx = methodIdxs(m);
        curMethodRE = methodREs(m);
                if strcmp(curMethodRE,"G")
            
            [jPDFAll,jPDFTonalOnly] = computeStriationPDFGeneralized2(IValid,...
                                            isTones,tonalNeighborIdx,MM,noiseVar);
            logLikelihoodAll(bIdx,validStriations,curMethodIdx) = jPDFAll;
            logLikelihoodAll(bIdx,validStriations,curMethodIdx+1) = jPDFTonalOnly;
        elseif strcmp(curMethodRE,"BT")
            [jPDFAll,jPDFTonalOnly] = computeStriationPDFBroadTonalHorizontal(IValid,...
                                            IRatio,isTones,tonalNeighborIdx);
            logLikelihoodAll(bIdx,validStriations,curMethodIdx) = jPDFAll;
            logLikelihoodAll(bIdx,validStriations,curMethodIdx+1) = jPDFTonalOnly;
        elseif strcmp(curMethodRE,"BB")
            [jPDFAll] = computeStriationPDFBroadOnly(IValid,IRatio,isTones,MM);
            logLikelihoodAll(bIdx,validStriations,curMethodIdx) = jPDFAll;
        elseif strcmp(curMethodRE,"TonalFast")
            [jPDFAll] = computeStriationPDFTonalOnlyFast(IValid,toneIdx,numTones,isTones,MM);
            logLikelihoodAll(bIdx,validStriations,curMethodIdx) = jPDFAll;
        elseif strcmp(curMethodRE,"XCORR")
            df = 1;
            xcorrCoeff(bIdx,:) = computeCorrCoeffRangeBroadband(IValid(:,isTones),numTones,df,numComb);    

        elseif strcmp(curMethodRE,"TonalOG")
            jPDFAll = computeStriationPDFTonalOnlyOriginal(IValid,toneIdx,numTones,isTones,y2Hyps);
            logLikelihoodAll(bIdx,validStriations,curMethodIdx) = jPDFAll;
        end
    end
end

for m=1:numMethods
    curMethodIdx = methodIdxs(m);
    curMethodRE = methodREs(m);
    tic
    if strcmp(curMethodRE,"G") || strcmp(curMethodRE,"BT")
        ll1 = squeeze(logLikelihoodAll(:,:,curMethodIdx));
        ll2 = squeeze(logLikelihoodAll(:,:,curMethodIdx+1));
        ll1(ll1 == 0) = nan; % ar
        ll2(ll2 == 0) = nan;
        % startIdx = round(size(ll1,2)-(striationIdxStartAll(1)+numel(rPotential)));
        llAligned1 = nan(size(ll1));  %,1),size(ll1,2)+startIdx);
        llAligned2 = nan(size(ll2)); % ,1),size(ll1,2)+startIdx);
        minSNum = min(striationIdxStartAll);
        for i=1:numBetas
            ds1 = striationIdxStartAll(i)-minSNum; % striationIdxStartAll(1);
            llAligned1(i,1:(end-ds1)) = ll1(i,ds1+1:end);
            llAligned2(i,1:(end-ds1)) = ll2(i,ds1+1:end);
        end
        validStriations = sum(isnan(llAligned1)) == 0;
        llValid1 = llAligned1(:,validStriations);                
        llValid2 = llAligned2(:,validStriations);

        LLAll1 = sum(llValid1,2);
        LLAll2 = sum(llValid2,2);

        [~,bEstIdx1] = max(LLAll1);
        [~,bEstIdx2] = max(LLAll2);
        bMLAllMethods(curMethodIdx) = betas(bEstIdx1);
        bMLAllMethods(curMethodIdx+1) = betas(bEstIdx2);
    elseif strcmp(curMethodRE,"BB") || strcmp(curMethodRE,"TonalFast")  || strcmp(curMethodRE,"TonalOG")
        ll1 = logLikelihoodAll(:,:,curMethodIdx);
        ll1(ll1 == 0) = nan; % ar
        llAligned1 = nan(size(ll1));  %,1),size(ll1,2)+startIdx);
        minSNum = min(striationIdxStartAll);

        for i=1:numBetas
            ds1 = striationIdxStartAll(i)-minSNum; % striationIdxStartAll(1);
            llAligned1(i,1:(end-ds1)) = ll1(i,ds1+1:end);
        end
        validStriations = sum(isnan(llAligned1)) == 0;
        llValid1 = llAligned1(:,validStriations);       

        LLAll = sum(llValid1,2);
        [~,bEstIdx1] = max(LLAll);
        bMLAllMethods(curMethodIdx) = betas(bEstIdx1);
    elseif strcmp(curMethodRE,"XCORR")
        xcorrCoeffmean = mean(xcorrCoeff,2,'omitnan');
        [~,bEstIdxCorr] = max(xcorrCoeffmean);
        bMLAllMethods(curMethodIdx) = betas(bEstIdxCorr);
    end

end

end


