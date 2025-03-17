% mainWIPairAnalysis.m
% 
% Extract the pair-wise mode contribution to the spectrogram based on the
% method presented in Le Gall and Bonnel (2013)
%
% author: Junsu Jang (junsu.jang94@gmail.com) 
% date: 2025/01/25
clear; % close all;

dataDir = '../AcousticData/dataByMinute2';
resultsDir = '../Results/RealData2';
filePostfix = "20s_15minLong";
isSaveFigure = true;

%% Parameters
rrAvg = 10.2;       % Average range rate
dt = 10;           % Snapshot step size in seconds
dr = rrAvg*dt;
numSpectrograms = 60;
MIN2SEC = 60;

% b is the hypothetical WI values correpsonding to the slope angles
b = nan(3,numSpectrograms);
% m is the Radon transform magnitude
m = nan(3,numSpectrograms);
rTrueRef = nan(1,numSpectrograms);

for n=3:numSpectrograms
    %% Get the relevant information from the data file
    refTimePassedSec = (n-1+15)*MIN2SEC;
    yr = 2017; mm = 3; dd = 24; hr = 18; min = 30;
    allDataStartTime = datetime(yr,mm,dd,hr,min,0);
    endTime = allDataStartTime + seconds(refTimePassedSec);
    filename = sprintf("%s/DataK_%02d%02d_%s.mat",dataDir,hour(endTime),minute(endTime),filePostfix);
    A = load(filename);
    % Only use 15 minutes worth of spectrogram
    % A.Z = A.Z(:,(end-90+1):end);
    % A.rTrue = A.rTrue((end-90+1):end);
    % A.t = A.t((end-90+1):end);

    snapshotStepsize = A.t(2)-A.t(1);
    f=A.f; f1 = A.f1; f2=A.f2;
    fFilter = f >= f1 & f <= f2;
    broadbandF = f(fFilter);    
    
    %% Pre-processing of the spectrogram
    curZ = A.Z(fFilter,:).';
    curM = abs(curZ);
    curIo = curM.^2;
    rNew = A.rTrue(1):(rrAvg*dt):A.rTrue(end);
    % interpolate the spectrogram to an even range value
    curI = interp1(A.rTrue,curIo,rNew);
    
    magSum = sum(curI);
    magRatio = magSum./max(magSum);
    
    IScaled = curI./repmat(magRatio,[size(curI,1),1]);
    
    % x-axis: range, y-axis: frequency 
    % Also remove the effect of the decay due to range
    rScale = repmat(rNew.',[1,size(IScaled,2)]);
    IScaled = rScale.*IScaled;
    
    %% Compute the bins of the FFT of the range (i.e., wavenumber)
    Fs_R = 1/dr;
    NFFT_R = size(IScaled,1);
    kr = (0:(NFFT_R-1))*Fs_R/NFFT_R;
    demeanedI = IScaled-repmat(mean(IScaled),[size(IScaled,1),1]);
    windowedI = demeanedI.*repmat(hamming(size(IScaled,1)),[1,size(IScaled,2)]);
    TFI = abs(fft(windowedI)).';
    TFI(:,kr >= Fs_R/2) = [];
    % TFI(:,1:2) = [];
    kr = 2*pi*kr(kr<Fs_R/2);
    % kr(1:2) = [];
    
    if isSaveFigure
        if ismember(n,[11,26,61])
            % kr vs f figure
            mFilename = convertStringsToChars(sprintf("%s/RDC_%d_%s",resultsDir,n,filePostfix));
            save(mFilename,'kr','broadbandF','TFI');

        end
    end


    %% Sample in the log domain    
    [KR,F] = meshgrid(kr,broadbandF);
    ILogInterp = scatteredInterpolant(KR(:),F(:),TFI(:),'linear','none');
    
    samp_dkr=0.01;
    samp_df =0.002;
    log_kr = log10(0.005):samp_dkr:log10(kr(end));
    log_f = log10(broadbandF(1)):samp_df:log10(broadbandF(end));
    
    samp_kr = 10.^(log_kr);
    samp_f = 10.^(log_f);
    [sT,sK] = meshgrid(samp_kr,samp_f);
    ILog = ILogInterp(sT,sK);

    % figure; imagesc(log_kr,log_f,10*log10(ILog)); colorbar
    % figure; imagesc(ILog);
    
    %% Radon transform
    dTheta = 0.1;
    % thetas = -180:dTheta:-5;
    % [R,xp] = radon(ILog.',thetas);
    % % Select the three mode pairs in the Radon domain
    % betaTheta = -tand(thetas)*(samp_df/samp_dkr);
    thetas = 91:dTheta:179; % -180:dTheta:-5;
    [R,xp] = radon(ILog.',thetas);

    % Select the three mode pairs in the Radon domain
    % delta kr is negative valued so multiply by -1 in front of samp_dkr
    % subtract 90 deg from thetas because that is the projection angle not
    % the actual angle of the slope relative to f
    slopekf = tand(thetas-90)*(-samp_dkr/samp_df); % = -1/beta
    betaTheta = -1./slopekf;



    tFilt = betaTheta >= 0.5 & betaTheta < 1.5;
    betaTheta = betaTheta(tFilt);
    % Three modes appear within the following windows
    % xpFilt1 = xp >= -10 & xp <= 10;
    % xpFilt2 = xp >= 23 & xp <= 33;
    % xpFilt3 = xp >= 43 & xp <= 53;
    xpFilt1 = xp >= -40 & xp < -25;
    xpFilt2 = xp >= -25 & xp < 0;
    xpFilt3 = xp >= 0 & xp < 40;
    % figure; imagesc(thetas(tFilt)-90,xp,10*log10(R(:,tFilt))); colorbar; clim([170 185]);
    % hold on;
    % yline(-40,'--r','Linewidth',2);
    % yline(-25,'--r','Linewidth',2);
    % yline(-25,'--w','Linewidth',2);
    % yline(0,'--w','Linewidth',2);
    % yline(0,'--y','Linewidth',2);
    % yline(40,'--y','Linewidth',2);


    R1 = R(xpFilt1,tFilt);
    R2 = R(xpFilt2,tFilt);
    R3 = R(xpFilt3,tFilt);
    xp1 = xp(xpFilt1);
    xp2 = xp(xpFilt2);
    
    % Find the xp that yields the peak in the filtered data
    [~,mIdx1] = max(max(R1.'));
    [~,mIdx2] = max(max(R2.'));
    [~,mIdx3] = max(max(R3.'));
    R1 = R1(mIdx1,:);
    R2 = R2(mIdx2,:);
    R3 = R3(mIdx3,:);

    [m1,b1] = max(R1);
    [m2,b2] = max(R2);
    [m3,b3] = max(R3);
    b(:,n) = betaTheta([b1,b2,b3]);
    m(:,n) = [m1;m2;m3];
    % wiObjFnFilename = sprintf('%s/WI_objFn/objFn_no%d.mat',resultsDir,n);
    % save(wiObjFnFilename,'m','betaTheta','R1','R2','R3');
 
    rTrueRef(n) = A.rTrue(end);
end
mm = m./repmat(sum(m,'omitnan'),[3,1]);
bm = sum(b.*mm,'omitnan');
% 
% figure; plot(movmean(b.',1));
% figure; plot(movmean(bm,5));
% legend(["R1: mode12"; "R2: mode23"; "R3: mode23"])

save(sprintf('%s/PairModesResult20250306.mat',resultsDir),'m','mm','b','bm','rTrueRef');