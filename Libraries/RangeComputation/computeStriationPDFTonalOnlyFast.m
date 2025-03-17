% Junsu Jang (junsu.jang94@gmail.com)
% 2025/01/20
function weights = computeStriationPDFTonalOnlyFast(IValid,toneIdx,numTones,isTones,MM)
% Implementation of the fast version of the tonal-only approach
noiseVar = estimateNoiseMLwithI(numTones,1,toneIdx,IValid);
Itones = IValid(:,isTones);
% instantaenous signal-to-broadband component ratio
ISBR = Itones./repmat(noiseVar.',[size(Itones,1),1]);

% Eq. 19 in Jang et al. (2025) to obtain the NC2Param (lambda)
Denom = mean(ISBR,'all') - 2;
Nume1 = mean(ISBR) - 2;
Nume2 = mean(ISBR,2) - 2;
NC2Param = (Nume2*Nume1)/Denom;
NC2Param(isnan(NC2Param) | NC2Param < 0) = 0;
% Apply moving mean striation wise to smooth out the noise
NC2Param = movmean(NC2Param,MM);

% Compute the PDF of the observed data based on the estimated NC2Param
besselInput = sqrt(NC2Param.*ISBR);
I0= besseli(0,besselInput); % Modified Bessel fn of the first kind of 0th order
% At high SNR, the bessel function requires approximation
bigI0idx = NC2Param > 100;
I0Big = exp(besselInput)./sqrt(2*pi*besselInput);
I0(bigI0idx) = I0Big(bigI0idx);

pdf = 1/2*exp(-(ISBR+NC2Param)./2).*I0;
logpdf = log(pdf);
weights = sum(logpdf,2,'omitnan');

end
