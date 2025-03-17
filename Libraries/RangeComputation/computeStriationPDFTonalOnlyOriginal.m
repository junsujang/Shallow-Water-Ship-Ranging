% Junsu Jang (junsu.jang94@gmail.com)
% 2025/01/20
function weights = computeStriationPDFTonalOnlyOriginal(IValid,toneIdx,numTones,isTones,y2Hyps)
% Implementation of the tonal-based statistical approach proposed by Young
% et al (2020). It requires 'searching' through the potential noncentrality
% parameters (y2hyp). Here, the y2hyp value that yields the highest PDF is
% selected.

noiseVar = estimateNoiseMLwithI(numTones,1,toneIdx,IValid);
Itones = IValid(:,isTones);
interp_y2 = Itones./repmat(noiseVar.',[size(Itones,1),1]);

ISum = sum(IValid);
MSum = sqrt(ISum(:,isTones));
ncParamRatio = (MSum/(max(MSum))).^2;
ncParam = (y2Hyps.')*ncParamRatio; % lambda_hyp (noncentrality parameter search)

numStriations = size(interp_y2,1);

ncParam = permute(repmat(ncParam,[1,1,numStriations]),[3,2,1]);
numHyp = size(y2Hyps,2);
sqrtAc2 = sqrt(ncParam);
y = sqrt(interp_y2);   % Rician requires y=sqrt(x1^2+x2^2), y2=y.^2
bigI0idx = sqrtAc2 > 10;

y = repmat(y,[1,1,numHyp]);
besselInput= sqrtAc2.*y;
I0 = besseli(0,besselInput); % Modified Bessel fn of the first kind of 0th order
I0Log = log(I0);
% I0Big = exp(besselInput)./sqrt(2*pi*besselInput);
I0BigLog = besselInput - .5*log(2*pi*besselInput);
% I0(bigI0idx) = I0Big(bigI0idx);
I0Log(bigI0idx) = I0BigLog(bigI0idx);
% pdfCandidates = .5*exp(-(interp_y2+ncParam)./2).*I0;
logPDF = log(.5)-.5*(interp_y2+ncParam)+I0Log;    
llMax = max(squeeze(sum(logPDF,2)),[],2);

weights = llMax;

end
