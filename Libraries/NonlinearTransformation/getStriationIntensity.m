% Junsu Jang (junsu.jang94@gmail.com)
% 2025/01/20
function IInterp = getStriationIntensity(f,betaML,r,X,refFIdx)
% Given the range and the reference frequency, the expected ranges along
% each striation are projected using the WI equation. Then, the intensities
% of the spectrogram are interpolated at those expected ranges to obtain
% the data in striation-frequency domain.
fIdx = 1:numel(f);
IInterp = nan(numel(r),size(X,2));
projected_r = r.'*((f(fIdx)./f(refFIdx)).^(1/betaML));
% The data at the reference frequency is left alone
fIdx(refFIdx) = [];
% For faster interpolation, we use ScaleTime cited below:
%   Jan (2025). ScaleTime 
%   (https://www.mathworks.com/matlabcentral/fileexchange/25463-scaletime), 
%   MATLAB Central File Exchange. Retrieved January 20, 2025.
for i=fIdx
    IInterp(:,i) = ScaleTimeNU(r,X(:,i),projected_r(:,i),'extrap');
end
IInterp(:,refFIdx) = X(:,refFIdx);

end