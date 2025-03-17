% mainWI_SBCEX17Kalamata_AllMethods.m
% 
% Process the SBCEX17 acoustic data to performing WI estimation using the
% following methods:
%   G:          The proposed generalized approach (Jang et al., 2025)
%   BT:         Legacy code that only considers the scenario without the background
%               noise (includes the tonals)
%   BB:         Broadband based approach (Jang and Meyer, 2024) 
%               This is referred to as (B) in Jang and Meyer (2025). 
%   TonalFast:  Faster version of the method proposed in Young et. al.
%               (2020). Only uses the tonal signal and assumes high background noise.
%   XCORR:      cross-correlation coefficient approach (Jang and Meyer, 2023)
%               This is referred to as (C) in Jang and Meyer (2025).
%   TonalOG:    The original version of Young et. al. (2020)
%               This is referred to as (T) in Jang and Meyer (2025).
% Unlike in mainRE, the range is known and the values of hypothetical WI
% are searched over.
%
% author: Junsu Jang (junsu.jang94@gmail.com) 
% date: 2025/01/25
clear; rng(7);

addpath('../Libraries/WIComputation');
addpath('../Libraries/ScaleTime_20201020'); % Required for fast interpolation
addpath('../Libraries/SBCEX17'); 
addpath('../Libraries/AcousticProcessing');
addpath('../Libraries/NonlinearTransformation');

dataDir = '../AcousticData/dataByMinute';

% Include methodRE "TonalOG" and methodIdx 8 to test TonalOG
methodREs = ["G"; "BT"; "BB"; "TonalFast"; "XCORR"];%; "TonalOG"];
% For testing purpose, G and BT also compute likelihoods using (1) the
% broadband and tonal component and (2) only the tonal component. 
% There are additional methodIdx 2 and 4, respectively, as a result.
methodIdxs = [1 3 5 6 7];% 8];
numMethods = numel(methodREs);

%% Parameters
betas = 0.8:0.01:1.3;
numSpectrograms = 61;
MIN2SEC = 60;
snapshotStepSize = 10.2*10; % range rate (m/s) * step size in time (s)
numTones = 5;
MM = 3;             % movingmean sample length
isEvenFreq = true; removeLeakage = true;
% The hypothetical search range for the noncentral parameter 
nc2ParamHyp = 0:0.05:15;
% background noise variance in linear intensity, measured at a different time
backgroundNoiseVar = 10^(82/20); % linear intensity
resultName = sprintf("../Results/RealData2/result_WI_allTime_MM%d_noAY.mat",MM);

%% Filtered Acoustic Data information
yr = 2017; mm = 3; dd = 24;
startHrs = [18 18 19 19];
startMins = [30 45 00 15];
dataTime = NaT(4,1);
for dataIdx=1
    hr = startHrs(dataIdx);
    curMin = startMins(dataIdx);
    dataTime(dataIdx) = datetime(yr,mm,dd,hr,curMin,0);
end


%% Process
rTrueRef = nan(numSpectrograms,1);
% WI computed from the range based on the AIS
bML_AIS = nan(methodIdxs(end),numSpectrograms);
% WI computed from the range based on the average range rate 
bML_AVG = nan(methodIdxs(end),numSpectrograms);


for n=3:numSpectrograms
    fprintf("Time step, n=%d\n",n);

    refTimePassedSec = (n-1+15)*MIN2SEC;
    endTime = dataTime(1) + seconds(refTimePassedSec);
    filename = sprintf("%s2/DataK_%02d%02d_20s_half.mat",dataDir,hour(endTime),minute(endTime));
    
    load(filename,'f','f1','f2','Z','rTrue');
    
    [fInterest,isTones,tonalF,curZ] = getTonalFrequencies(f,f1,f2,isEvenFreq,removeLeakage,Z);
    % minimum number of frequencies that are required to have valid intensities
    minNumFreqProjections = numel(fInterest); 
    % the reference frequency for nonlinear transformation
    refFreqIdx = floor(numel(fInterest)/2);     
    
    % curZ = Z(fFilter,:).';
    
    bMLAllMethodsAIS = estimateBetaAllMethods(curZ,rTrue,fInterest,refFreqIdx,...
                                            betas,minNumFreqProjections,isTones,...
                                            methodREs,methodIdxs,MM,nc2ParamHyp,backgroundNoiseVar);

    bML_AIS(:,n) = bMLAllMethodsAIS;

    % rAvg = rTrue(end) + (-snapshotStepSize*(numel(rTrue)-1):snapshotStepSize:0);
    % bMLAllMethodsAvg = estimateBetaAllMethods(curZ,rAvg,fInterest,refFreqIdx,...
    %                                         betas,minNumFreqProjections,isTones,...
    %                                         methodREs,methodIdxs,MM,nc2ParamHyp,backgroundNoiseVar);
    % 
    % bML_AVG(:,n) = bMLAllMethodsAvg;

    rTrueRef(n) = rTrue(end);
    
end


save(resultName,'bML_AIS','bML_AVG','rTrueRef','methodREs'); 
