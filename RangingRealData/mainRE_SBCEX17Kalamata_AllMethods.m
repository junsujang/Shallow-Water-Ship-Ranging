% mainRE_SBCEX17Kalamata_AllMethods.m
% 
% Process the SBCEX17 acoustic data to performing ranging using the
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
% author: Junsu Jang (junsu.jang94@gmail.com) 
% date: 2025/01/25
clear; rng(7);

addpath('../Libraries/RangeComputation');
addpath('../Libraries/ScaleTime_20201020'); % Required for fast interpolation
addpath('../Libraries/SBCEX17'); 
addpath('../Libraries/AcousticProcessing');
addpath('../Libraries/NonlinearTransformation');

dataDir = '../AcousticData/dataByMinute2';

% Include methodRE "TonalOG" and methodIdx 8 to test TonalOG
methodREs = ["G"; "BT"; "BB"; "TonalFast"; "XCORR"; "TonalOG"];
% For testing purpose, G and BT also compute likelihoods using (1) the
% broadband and tonal component and (2) only the tonal component. 
% There are additional methodIdx 2 and 4, respectively, as a result.
methodIdxs = [1 3 5 6 7 8];
numMethods = methodIdxs(end); % numel(methodREs);
MName = "20s_half"; % "20s_half"; % "onethird"; % "half" for five min or nothing for 10 min worht of striations
if contains(MName,"10s")
    snapshotInterval = 5;
elseif contains(MName,"16s")
    snapshotInterval = 8;
elseif contains(MName,"20s")
    snapshotInterval = 10;
elseif contains(MName,"30s")
    snapshotInterval = 15;
end


%% Parameters
numSpectrograms = 61;
MIN2SEC = 60;
rrAvg = 10.2;
numTones = 5;
rHypSpan = 0.4;     % The range of search range (+/-0.4 of the true r value)
drHyp = 10;         % search range resolution
MM = 3;             % movingmean sample number
isEvenFreq = true; % isEven=true removes broadband frequencies that are not multiples of 0.2 Hz
removeLeakage = true;
betaML = 1.18;
% The hypothetical search range for the noncentral parameter 
nc2ParamHyp = 0:0.05:100;            
% background noise variance in linear intensity, measured at a different time
backgroundNoiseVar = 10^(82/20);    


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


%% If we want to save the log likelihood data cross all data, 
% identify the smalles and largest range
refTimePassedSec1 = (26-1+15)*MIN2SEC;
endTime1 = dataTime(1) + seconds(refTimePassedSec1);
filename1 = sprintf("%s/DataK_%02d%02d_%s.mat",dataDir,hour(endTime1),minute(endTime1),MName);
load(filename1,'rTrue','Z');
rS = (rTrue(end)*(1-rHypSpan))-mod((rTrue(end)*(1-rHypSpan)),drHyp);

refTimePassedSec2 = (numSpectrograms-1+15)*MIN2SEC;
endTime2 = dataTime(1) + seconds(refTimePassedSec2);
filename2 = sprintf("%s/DataK_%02d%02d_%s.mat",dataDir,hour(endTime2),minute(endTime2),MName);
load(filename2,'rTrue','Z');
rE = (rTrue(end)*(1+rHypSpan))-mod((rTrue(end)*(1+rHypSpan)),drHyp);

totalRHyp = rS:drHyp:rE;
allTime = NaT(numSpectrograms,1); % for saving the timestamp of the data


%% Process
rTrueRef = nan(numSpectrograms,1);
rML = nan(methodIdxs(end),numSpectrograms);
LLAll = nan(numel(totalRHyp),numSpectrograms);
resultName = sprintf("../Results/RealData2/result_RE_allTimeNew_MM%d_b%d_%s_3BB.mat",MM,betaML*100,MName);
fInterestAll = [];
tMAll = nan(numMethods+1,numSpectrograms); % Last index for method is the NLT duration
for n=26:numSpectrograms  % 26:numSpectrograms
    refTimePassedSec = (n-1+15)*MIN2SEC;
    endTime = dataTime(1) + seconds(refTimePassedSec);
    filename = sprintf("%s/DataK_%02d%02d_%s.mat",dataDir,hour(endTime),minute(endTime),MName);
    % filename = sprintf("%s/DataK_%02d%02d.mat",dataDir,hour(endTime),minute(endTime));

    load(filename,'f','f1','f2','Z','rTrue','rr');
    rPotential = (rrAvg*snapshotInterval)*(-(numel(rTrue)-1):0);
    % set the rHyp for this data
    rS = (rTrue(end)*(1-rHypSpan))-mod((rTrue(end)*(1-rHypSpan)),10);
    rE = (rTrue(end)*(1+rHypSpan))-mod((rTrue(end)*(1+rHypSpan)),10);
    rHyp = rS:drHyp:rE; 

    [fInterest,isTones,tonalF,curZ] = getTonalFrequencies(f,f1,f2,isEvenFreq,removeLeakage,Z);
    minNumFreqProjections = numel(fInterest);   % minimum number of frequencies that are required to have valid intensities
    refFreqIdx = floor(numel(fInterest)/2);     % the reference frequency for nonlinear transformation
    % fInterest
    % curZ = Z(fFilter,:).';
    fprintf("Time step, n=%d\n",n);
    % rHyp(1) = rTrue(end);
    [rMLAllMethods,tM,LLAll2] = estimateRangeAllMethods(curZ,rPotential,rHyp,fInterest,refFreqIdx,...
                                            betaML,minNumFreqProjections,isTones,...
                                            methodREs,methodIdxs,MM,nc2ParamHyp,backgroundNoiseVar);
    rML(:,n) = rMLAllMethods;
    % tonalF
    rTrueRef(n) = rTrue(end);
    
    % Align the current rShift with the global rShift range.
    [~,rMinIdx] = min(abs(totalRHyp - rHyp(1)));
    LLAll(rMinIdx:(rMinIdx+numel(rHyp)-1),n) = LLAll2;
    allTime(n) = endTime;
    fInterestAll = [fInterestAll; fInterest];
    tMAll(:,n) = tM;
end

save(resultName,'rML','rTrueRef','methodREs','LLAll','allTime','totalRHyp','fInterestAll','tonalF','tMAll'); 


%% Figures to visualize the results 
errAll = rML-rTrueRef.'; 
errAllP = 100*errAll([1,5,6,7,8],:)./rTrueRef.';
validData = 26:numSpectrograms;
figure; plot(1:numSpectrograms,errAllP,'LineWidth',2); % movmean(rML([2,5,6,7],:).'./1e3-rTrueRef./1e3,1),'Linewidth',2); 

% figure; plot(rTrueRef./1e3,errAllP,'LineWidth',2); % movmean(rML([2,5,6,7],:).'./1e3-rTrueRef./1e3,1),'Linewidth',2); 
legend(methodREs([1 3 4 5 6 7]),'Location','northwest'); 
xlabel("Range (km)"); ylabel("Error (percent)"); 
grid on; box on; 
set(gca,'FontSize',20); 
title("b=1.18,movmean=1");

figure; plot(rTrueRef./1e3,movmean(abs(rML([2,5,6,7],:).'./1e3-rTrueRef./1e3),5),'Linewidth',2); 
legend(methodREs([1 3 4 5]),'Location','northwest'); 
xlabel("Range (km)"); ylabel("Error (km)"); 
grid on; box on; 
xlim([rTrueRef(validData(1)),rTrueRef(end)]./1e3); 
set(gca,'FontSize',20); 
title("b=1.18,movmean=5");