% Junsu Jang (junsu.jang94@gmail.com)
% 2025/01/20
function [weightsAll,weightsTonalOnly] = computeStriationPDFGeneralized2(IValid,...
                                isTones,tonalNeighborIdx,MM,noiseVar)
% Implementation of the generalized approach proposed in Jang et al. (2025)

isBroadband = ~isTones;
IValidBB = IValid(:,isBroadband); 
ISampMean = mean(IValid);
% the total noiseVar is 2*STD^2 and this value is assumed to be known
SignalOnly = ISampMean - noiseVar;
SignalOnly(SignalOnly <= 0) = 0;
maxI = max(SignalOnly(~isTones));
% Compute the propagated broadband scaling factor \alpha (Eq. 15 in Jang
% et al. (2025)
SignalRatio = SignalOnly/maxI;
SignalRatio(isTones) = 1;% Scale the rest (broadband) signals

%% Broadband processing (Exponential)
% Compute the propagated broadband signal variance v (Eq. 16)
striationSampMean = mean(IValidBB,2,'omitnan');
striationSignalParam = (striationSampMean - noiseVar)./mean(SignalRatio(~isTones));
striationSignalParam(isnan(striationSignalParam) | striationSignalParam < 0) = 0;
striationSignalParam = movmean(striationSignalParam,MM);

% Compute the exponential parameter (Eq. 17)
expParam = striationSignalParam*SignalRatio + noiseVar;

% Compute the exponential PDF in case needed
IScaledBB = IValidBB./expParam(:,~isTones);
pdfBB = exp(-IScaledBB)./expParam(:,~isTones);


%% Tonal processing 
IValidBT = IValid(:,isTones);
tonalVarEst = squeeze(mean(reshape(expParam(:,tonalNeighborIdx),size(expParam,1),2,[]),2));
% The variance of each of the real and imaginary parts of the acoustic
% field is half of the exponential parameter
ISBR = IValidBT./(tonalVarEst/2); 

% Eq. 19 in Jang et al. (2025) to obtain the NC2Param (lambda)
Denom = mean(ISBR,'all') - 2;
Nume1 = mean(ISBR) - 2;
Nume2 = mean(ISBR,2) - 2;
NC2Param = (Nume2*Nume1)/Denom;
NC2Param(isnan(NC2Param) | NC2Param < 0) = 0;
% Apply moving mean striation wise to smooth out the noise
NC2Param2 = movmean(NC2Param,MM);
% [~,pkIdx] = findpeaks(-NC2Param2(:,1),'NPeaks',10);
nonCentralChi2PDF= ncx2pdf(ISBR, 2, NC2Param2);

% % Compute the PDF of the observed data based on the estimated NC2Param
% besselInput = sqrt(NC2Param.*ISBR);
% I0= besseli(0,besselInput); % Modified Bessel fn of the first kind of 0th order
% % At high SNR, the bessel function requires approximation
% bigI0idx = NC2Param > 100;
% I0Big = exp(besselInput)./sqrt(2*pi*besselInput);
% I0(bigI0idx) = I0Big(bigI0idx);
% 
% nonCentralChi2PDF = 1/2*exp(-(ISBR+NC2Param)./2).*I0;

%% Combine the PDFs
% Broadband + tonal PDF
pdfAllFreq = [pdfBB, nonCentralChi2PDF];
logpdf2 = log(pdfAllFreq);
weightsAll = sum(logpdf2,2);

% Tonal only componenet (used in the paper)
pdfTonalsOnly = nonCentralChi2PDF;
logpdf1 = log(pdfTonalsOnly);
weightsTonalOnly = sum(logpdf1,2); % (:,[1 2 3 5]),2);
% weightsTonalOnly(~ismember(1:size(weightsTonalOnly,1),pkIdx),:) = 0;
end