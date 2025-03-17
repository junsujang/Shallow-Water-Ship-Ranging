% Junsu Jang (junsu.jang94@gmail.com)
% 2025/01/20
function [rMLAllMethods,tMethods,LLAll2] = estimateRangeAllMethods(Z,rPotential,rHyp,f,refFreqIdx,...
                                                betaML,minNumF,isTones,...
                                                methodREs,methodIdxs,MM,y2Hyps,noiseVar)
% In this function, the loglikelihood for each rHyp per method are computed.
% First, the likelihood for each rHyp is obtained, and then in the last
% for-loop, the maximum likelihood estimates are obtained per method.
% 
% Since the nonlinear transformation, i.e., the interpolation, took the
% longest time but was commonly needed for multiple methods, we perform the
% interpolation first. 

%% Parameters
numMethods = numel(methodREs);
numRHyp = numel(rHyp);
% Identify the indicies of the neighboring broadband signals
toneIdx = find(isTones);
tonalNeighborIdx = [toneIdx-1; toneIdx + 1];
numTones = sum(isTones);
numComb = nchoosek(numTones,2); % number of pair-wise combinations expected for xcorr

% Initialize variables
logLikelihoodAll = nan(numRHyp,size(Z,1),6);
rMLAllMethods = nan(methodIdxs(end),1);
tMethods = zeros(methodIdxs(end)+1,1);  % also measure the time elapsed for each method
xcorrCoeff = nan(numRHyp,numComb); % This is the objective function per pair of tonals in XCORR
% The index of the striations after performing the nonlinear transformation
% does not actually match the index of the spectrogram due to the way it's
% been optimized. However, to utilize the common striations, we would like
% to align the striations. Therefor, we keep track of the index of the
% striation that corresponds to the first spectrogram snapshot.
striationIdxStartAll = nan(numRHyp,1);

I = abs(Z).^2;

%% Compute the Loglikelihood or objective function per rHyp
for rShiftIdx=1:numRHyp
    % rPotential is the rr*t
    % the guess of the range for the given spectrogram
    curR = rPotential+rHyp(rShiftIdx);
    
    curI = I;
    tic
    % Perform the nonlinear transformation
    [IInterp,striationIdx] = computeNLTData(f,betaML,curR,curI,refFreqIdx);
    % Identify the valid striations
    validEntries = ~isnan(IInterp);
    numValidFreqPerStriation = sum(validEntries,2);
    validStriations = numValidFreqPerStriation >= minNumF;
    IValid = IInterp(validStriations,:);

    % For scaling and normalization of the intensities across frequencies
    % (needed by BT and B)
    avgMValid = mean(sqrt(IValid));
    [maxM,~] = max(avgMValid(~isTones));
    MRatio = avgMValid/(maxM);
    IRatio = MRatio.^2;
    IRatio(isTones) = 1;
    tMethods = tMethods+toc;
    
    % Keep track of the idx corresponding to the first snapshot of the
    % spectrogram
    striationIdxStartAll(rShiftIdx) = striationIdx;
    tMethods(end) = toc;
    %% Perform loglikelihood or objective function per method
    for m=1:numMethods
        curMethodIdx = methodIdxs(m);
        curMethodRE = methodREs(m);
        tic
        if strcmp(curMethodRE,"G")
            
            [jLLAll,jLLTonalOnly] = computeStriationPDFGeneralized2(IValid,...
                                            isTones,tonalNeighborIdx,MM,noiseVar);
            logLikelihoodAll(rShiftIdx,validStriations,curMethodIdx) = jLLAll;
            logLikelihoodAll(rShiftIdx,validStriations,curMethodIdx+1) = jLLTonalOnly;
        elseif strcmp(curMethodRE,"BT")
            [jLLAll,jLLTonalOnly] = computeStriationPDFBroadTonalHorizontal(IValid,...
                                            IRatio,isTones,tonalNeighborIdx);
            logLikelihoodAll(rShiftIdx,validStriations,curMethodIdx) = jLLAll;
            logLikelihoodAll(rShiftIdx,validStriations,curMethodIdx+1) = jLLTonalOnly;
        elseif strcmp(curMethodRE,"BB")
            [jLLAll] = computeStriationPDFBroadOnly(IValid,IRatio,isTones,MM);
            logLikelihoodAll(rShiftIdx,validStriations,curMethodIdx) = jLLAll;
        elseif strcmp(curMethodRE,"TonalFast")
            [jLLAll] = computeStriationPDFTonalOnlyFast(IValid,toneIdx,numTones,isTones,MM);
            logLikelihoodAll(rShiftIdx,validStriations,curMethodIdx) = jLLAll;
        elseif strcmp(curMethodRE,"XCORR")
            df = 1;
            xcorrCoeff(rShiftIdx,:) = computeCorrCoeffRangeBroadband(IValid(:,isTones),numTones,df,numComb);    

        elseif strcmp(curMethodRE,"TonalOG")
            jLLAll = computeStriationPDFTonalOnlyOriginal(IValid,toneIdx,numTones,isTones,y2Hyps);
            logLikelihoodAll(rShiftIdx,validStriations,curMethodIdx) = jLLAll;
        end
        tMethods(curMethodIdx) = tMethods(curMethodIdx)+toc;
    end
end

%% For each method, compute the joint loglikelihood and determine the range
for m=1:numMethods
    curMethodIdx = methodIdxs(m);
    curMethodRE = methodREs(m);
    tic
    if strcmp(curMethodRE,"G") || strcmp(curMethodRE,"BT")
        % 1 and 2 here mean two different types of loglikelihoods (both
        % broadband and tonal / tonal only)
        ll1 = squeeze(logLikelihoodAll(:,:,curMethodIdx));
        ll2 = squeeze(logLikelihoodAll(:,:,curMethodIdx+1));
        % ll1(ll1 == 0) = nan;
        % ll2(ll2 == 0) = nan;
        % Align the LL based on the idx of the first snapshot of the
        % spectrogram.
        llAligned1 = nan(size(ll1));
        llAligned2 = nan(size(ll2));
        minSNum = min(striationIdxStartAll);
        for i=1:numRHyp
            ds1 = striationIdxStartAll(i)-minSNum;
            llAligned1(i,1:(end-ds1)) = ll1(i,ds1+1:end);
            llAligned2(i,1:(end-ds1)) = ll2(i,ds1+1:end);
        end
        % Striations that have valid entries across all rHyps are used for
        % unbiased loglikelihood computation
        validStriations = sum(isnan(llAligned1)) == 0;
        llValid1 = llAligned1(:,validStriations);                
        llValid2 = llAligned2(:,validStriations);

        % The joint LL across all striations
        LLAll1 = sum(llValid1,2);
        LLAll2 = sum(llValid2,2);

        % Obatin the ML estimates
        [~,rEstIdx1] = max(LLAll1);
        [~,rEstIdx2] = max(LLAll2);
        rMLAllMethods(curMethodIdx) = rHyp(rEstIdx1);
        rMLAllMethods(curMethodIdx+1) = rHyp(rEstIdx2);
    elseif strcmp(curMethodRE,"BB") || strcmp(curMethodRE,"TonalFast")  || strcmp(curMethodRE,"TonalOG")
        ll = logLikelihoodAll(:,:,curMethodIdx);
        ll(ll == 0) = nan; 
        llAligned = nan(size(ll)); 
        minSNum = min(striationIdxStartAll);

        for i=1:numRHyp
            ds1 = striationIdxStartAll(i)-minSNum; % striationIdxStartAll(1);
            llAligned(i,1:(end-ds1)) = ll(i,ds1+1:end);
        end
        validStriations = sum(isnan(llAligned)) == 0;
        llValid1 = llAligned(:,validStriations);       

        LLAll = sum(llValid1,2);
        [~,rEstIdx1] = max(LLAll);
        rMLAllMethods(curMethodIdx) = rHyp(rEstIdx1);
    elseif strcmp(curMethodRE,"XCORR")
        % Compute the joint xcorr objective function and obtain the range
        xcorrCoeffmean = mean(xcorrCoeff,2,'omitnan');
        [~,rEstIdxCorr] = max(xcorrCoeffmean);
        rMLAllMethods(curMethodIdx) = rHyp(rEstIdxCorr);
    end
    tMethods(curMethodIdx) = tMethods(curMethodIdx)+toc;

end

end


