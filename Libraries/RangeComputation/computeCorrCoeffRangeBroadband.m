% Junsu Jang (junsu.jang94@gmail.com)
% 2025/01/20
function [xcorrCoeff] = computeCorrCoeffRangeBroadband(interp_I_valid,numF,df,numComb)
% Compute frequency pair-wise cross correlation coefficient (at lag=0)
% See Jang and Meyer (2023) for details
xcorrCoeff = zeros(numComb,1);
dIdx = 1;

for i=1:numF-df
    for j = (i+df):df:numF
        originalI = interp_I_valid(:,j);
        projectedI = interp_I_valid(:,i); 
        nanFilter = ~(isnan(originalI) | isnan(projectedI));
        xc = corrcoef(originalI(nanFilter),projectedI(nanFilter));
        xcorrCoeff(dIdx) = xc(1,2);
        dIdx = dIdx + 1;
    end
end
end