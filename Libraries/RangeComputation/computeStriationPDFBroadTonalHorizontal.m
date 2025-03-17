function[weightsAll,weightsTonalOnly] = computeStriationPDFBroadTonalHorizontal(IValid,...
                                IRatio,isTones,tonalNeighborIdx) 
% Identify the indicies of the neighboring broadband signals
isBroadband = ~isTones;  % member(1:numF,[toneIdx; toneIdx-1; toneIdx + 1].');
IValidBB = IValid(:,isBroadband);

% Scale the rest (broadband) signals
IScaled = IValidBB./repmat(IRatio(isBroadband),[size(IValidBB,1),1]);
MValidBB = sqrt(IScaled);

%% Broadband processing (Rayleigh)
% method of moment based estimator
rParamMLE = mean(MValidBB,2,'omitnan')/sqrt(pi/2);
% MLE estimate of the magnitude (known to be biased)
% rParamMLE = sqrt(mean((broadSignal.^2)./2,2,'omitnan'));
rParamMLE = repmat(rParamMLE,[1,size(MValidBB,2)]);

MValidBBNormalized = MValidBB./rParamMLE; % normalize the rParam so that it is all 1
IValidBBNormalized = MValidBBNormalized.^2; % this should turn the distribution to chi2 w/ 2DOF (or exponential with 1/2 param)
expParam = .5;
pdfBB = expParam*exp(-expParam*IValidBBNormalized); % chi2pdf(IValidBBNormalized,2); % get 2DOF chi2 PDF


% MValidBB = MValidBB./rParamMLE; % normalize the rParam so that it is all 1
% rayleighPDF = MValidBB.*exp(-(MValidBB.^2)./(2));

% % Obtain the PDF
% rayleighPDF = (broadbandMagnitude./(rParamMLE.^2)).*exp(-(broadbandMagnitude.^2)./(2*rParamMLE.^2));



%% Tonal processing (Rice)
tonalIntensity = IValid(:,isTones);
% IRatioAll = ones(numF,1); IRatioAll(isBroadband) = IRatio.';
IRatioAll = IRatio.';
IRatioTonal = mean(IRatioAll(tonalNeighborIdx));
tonalVarEst = IRatioTonal.*(rParamMLE(:,1).^2); %
tonalSNR = tonalIntensity./tonalVarEst; % movmean(tonalStdEst.^2,5); % (interp_M_all(:,isTones)./tonalStdEst).^2;

Ac2 = mean(tonalSNR,'omitnan')-2;
Ac2(Ac2 < 0) = 0;
Ac2 = repmat(Ac2, [size(tonalSNR,1),1]);
sqrtAc2 = sqrt(Ac2);

tonalMNR = sqrt(tonalSNR);   % Rician requires y=sqrt(x1^2+x2^2), y2=y.^2
besselInput = sqrtAc2.*tonalMNR;
I0= besseli(0,besselInput); % Modified Bessel fn of the first kind of 0th order
bigI0idx = sqrtAc2 > 10;
I0Big = exp(besselInput)./sqrt(2*pi*besselInput);
I0(bigI0idx) = I0Big(bigI0idx);

%% Non Central Chi^2 model
nonCentralChi2PDF = 1/2*exp(-(tonalSNR+Ac2)./2).*I0;
 

%% Combine the PDFs
pdfTonalsOnly = nonCentralChi2PDF;
pdfTonalsOnly(pdfTonalsOnly==0) = 1e-500;
logpdf1 = log(pdfTonalsOnly);
weightsTonalOnly = sum(logpdf1,2);

pdfAllFreq = [pdfBB, nonCentralChi2PDF];
% pdfAllFreqMin = nan(size(pdfAllFreq,1),minNumFreqProjections);
% for i=1:size(pdfAllFreq,1)
%     sIdx = find(~isnan(pdfAllFreq(i,:)));
%     randIdx = randperm(numel(sIdx),minNumFreqProjections);
%     pdfAllFreqMin(i,:) = pdfAllFreq(i,sIdx(randIdx));
% end
logpdf2 = log(pdfAllFreq);
weightsAll = sum(logpdf2,2);


end