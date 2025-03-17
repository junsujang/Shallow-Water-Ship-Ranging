% Junsu Jang (junsu.jang94@gmail.com)
% 2025/01/20
function weightsBB = computeStriationPDFBroadOnly(IValid,IRatio,isTones,MM)
% Implementation of the broadband only approach introduced in Jang and
% Meyer (2024)

isBroadband = ~isTones; 
MRatio = sqrt(IRatio);
MValidBB = sqrt(IValid(:,isBroadband));
MScaled = MValidBB./repmat(MRatio(isBroadband),[size(MValidBB,1),1]);

% MLE estimate of the magnitude (known to be biased, but OK so far)
% Also apply moving mean to reduce the effect of noise
rParamMLE = movmean(sqrt(mean((MScaled.^2)./2,2,'omitnan')),MM);
rParamMLE = repmat(rParamMLE,[1,size(rParamMLE,2)]);

% 'whiten' (or normalize) the data so it is less affected by the bias
% introduced at different values of the Rayleigh parameters
MUnitScale = MScaled./rParamMLE; % normalize the rParam so that it is all 1
rayleighPDF = MUnitScale.*exp(-(MUnitScale.^2)./(2));

%% Combine the PDFs
pdfBB = rayleighPDF;
logpdf = log(pdfBB);
weightsBB = sum(logpdf,2);

end