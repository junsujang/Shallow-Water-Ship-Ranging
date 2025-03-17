% Junsu Jang (junsu.jang94@gmail.com)
% 2025/01/20
function [fInterest,isTones,tonalF,curZ] = getTonalFrequencies(f,f1,f2,isEvenFreq,removeLeakage,Z)

%% Select the broadband at intervals outside of the tonal narrowband
% isEvenFreq and removeLeakage are legacy variables
fFilter = f >= f1 & f <= f2;
f = f(fFilter);
curZ = Z(fFilter,:);
PSD = mean(abs(curZ).^2,2);
% Find peaks of the PSD as the tonal central frequencies
[~,tonalF] = findpeaks(PSD,f,'NPeaks',5,'MinPeakDistance',1);

% Compute the broadband frequencies (three bins between two tonals)
td1 = diff(tonalF);
td2 = (1:3).'*(td1/4); 
beforeTonals = tonalF(1) - td2(:,1);
midTonal = tonalF(1:4) + td2;
afterTonals = tonalF(5) + td2(:,1);
broadbandF = [beforeTonals; midTonal(:); afterTonals].';
broadbandF(broadbandF < f1 | broadbandF > f2) = [];
fInterest = sort([broadbandF, tonalF]);

[~,fIdx] = min(abs(f - fInterest.'),[],2);
fIdx = uniquetol(fIdx);
fInterest=(f(fIdx)); % discretize

curZ = curZ(fIdx, :).';
isTones = ismember(fInterest,tonalF).';

end

