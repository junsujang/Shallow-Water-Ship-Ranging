% mainSlopeBased_SBCEX17.m
% 
% Process the SBCEX17 acoustic data to performing ranging using the 2D-DFT
% Radon transform method. The approach utilizes the slope of the striation.
% This correponds to method (S) in the paper Jang and Meyer (2025).
% 
% Based on the methods presented in Cockrell and Schmidt (2010) and 
% Yao et al. (2021)
%
% author: Junsu Jang (junsu.jang94@gmail.com) 
% date: 2025/01/25
clear; close all;

addpath('../Libraries/SBCEX17');
addpath('../Libraries/RangeComputation')

dataDir = '../AcousticData/dataByMinute2';
resultsDir = '../Results/RealData';
% resultName = sprintf("%s/resultRangeDFT_allTime_b118_20s_half.mat",resultsDir);
resultName = sprintf("%s/PairModesResult20250218.mat",resultsDir);

%% Parameters
numDataPoints = 61;
MIN2SEC = 60; % Convert minute to seconds
rrAvg = 10.2; % Average range rate used
numTones = 5; % number of tones in the data
betaML = 1.18;% assumed WI value 
% frequency extraction parameter - include the bins with leakage to ensure
% even frequency bins
isEvenFreq = true; removeLeakage = false;

cMin = 1450; % minimum sound speed feasible
cMax = 1750; % maximum soudn speed feasible

% search angles for slope calculation in the Radon domain
dTheta = 0.01;
thetas = -90:dTheta:(90-dTheta);


%% Filtered Acoustic Data information parameter setup
yr = 2017; mm = 3; dd = 24;
startHrs = [18 18 19 19];
startMins = [30 45 00 15];
dataTime = NaT(4,1);
for dataIdx=1:4
    hr = startHrs(dataIdx);
    curMin = startMins(dataIdx);
    dataTime(dataIdx) = datetime(yr,mm,dd,hr,curMin,0);
end

tElapsed = nan(numDataPoints,1);
%% Process
rTrueRef = nan(numDataPoints,1);
rS = nan(numDataPoints,1);
tAll = NaT(numDataPoints,1);
for n=26:numDataPoints
    %% Load the acoustic data
    refTimePassedSec = (n-1+15)*MIN2SEC;
    endTime = dataTime(1) + seconds(refTimePassedSec);
    % filename = sprintf("%s/DataK_%02d%02d.mat",dataDir,hour(endTime),minute(endTime));

    filename = sprintf("%s/DataK_%02d%02d_20s_half.mat",dataDir,hour(endTime),minute(endTime));
    snapshotInterval = 2.5;
    if contains(filename,"10s")
        snapshotInterval = 5;
    elseif contains(filename,"16s")
        snapshotInterval = 8;
    elseif contains(filename,"20s")
        snapshotInterval = 10;
    elseif contains(filename,"30s")
        snapshotInterval = 15;
    end

    load(filename,'f','f1','f2','Z','rTrue');
    % f2 = f2 + 0.2;
    % [~,isTones,tonalF,~] = getTonalFrequencies(f,f1,f2,isEvenFreq,removeLeakage,Z);
    fFilter = f >= f1 & f <= f2;
    fInterest = f(fFilter);
    curZ = Z(fFilter,:);
    stepsizeR = rrAvg*snapshotInterval; %(rPotential(2)-rPotential(1));
    stepSizeF = fInterest(2)-fInterest(1);
    rSearchMax = rTrue(end); % This is needed for the filtering in the 2D DFT domain
    
    tic
    %% Window the data
    % Compute the window size (Cockrell and Schmidt, 2010)
    % We ignore the frequency window size recommendation and use the full 7Hz
    % window
    kr_mn = 2*pi*fInterest(1)*(1/cMin-1/cMax)*1/3;
    win_r = 3*2*pi/kr_mn;
    win_r_px = ceil(win_r/stepsizeR);

    % Make sure that the windows are all odd numbers so that the radon
    % transform happens around the center
    if mod(win_r_px,2) == 0; win_r_px = win_r_px+1; end
    win_r = win_r_px*stepsizeR;

    rIdx2 = size(curZ,2);
    rIdx1 = rIdx2 - win_r_px + 1;
    curZ = curZ(:,rIdx1:rIdx2).';
    curI = abs(curZ).^2;

    % The signals are 'normalized' for less biased comparison
    ISum = sum(curI);
    IRatio = ISum./max(ISum);
    IScaled = (curI./repmat(IRatio,[size(curI,1),1])).';

    % Undo the following if testing 'unnormalized' (i.e., raw) intensity
    % signal
    % IScaled = curI; 

    [NFFT_F,NFFT_R] = size(IScaled);
    midIdx = floor(size(IScaled)./2); % Identify the center idx of 2D DFT
    midFreq = fInterest(midIdx(1)); 

    %% Compute the range
    rObjFn = computeRangeDFTRadonObjFn(IScaled,rSearchMax,fInterest,...
                                    cMin,cMax,thetas,stepsizeR,stepSizeF);
    % The slope of interest is that of the complement angle
    slopes = tand(90-thetas); % 
    % Apply the ratio of the range and frequency windows to get meaningful
    % slope.
    r0 = -betaML*midFreq*(NFFT_R*stepsizeR)/(NFFT_F*stepSizeF).*(slopes);
    r0 = r0 + stepsizeR*NFFT_R/2; % make the estimation w.r.t the middle of the spectrogram
    % Only perform the estimation within the search window of other
    % methods.
    positiveSlopeFilter = r0 > rSearchMax*0.6 & r0 < rSearchMax*1.4;
    r0 = r0(positiveSlopeFilter);
    rObjFn = rObjFn(positiveSlopeFilter);
    [~,maxIdx] = max(rObjFn);
    tElapsed(n) = toc;
    rS(n) = r0(maxIdx);
    rTrueRef(n) = rTrue(end);
    tAll(n) = endTime;
end
save(resultName,'rS','rTrueRef','tAll'); 
