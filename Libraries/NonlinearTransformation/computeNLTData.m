% Junsu Jang (junsu.jang94@gmail.com)
% 2025/01/20
function [IInterp,striationIdxStart] = computeNLTData(f,betaML,r,I,refFreqIdx)
    % Project the first range at high frequency to ref freq to get the
    % lowest range needed to make a full projection from middle to high
    % freq. Similarly project the last range at the lowest freq to ref freq
    % to get the highest range needed.
    %
    % note also that when non-uniform range rate is utilized, it is simpler
    % to interpolate the data such that the striations are uniformly
    % spaced.
    dr = max(abs(diff(r))); % assumes constant range rate of the last 
    projected_r_low = r(1)*((f(refFreqIdx)./f(end)).^(1/betaML));
    r_low = max([projected_r_low - mod(projected_r_low,dr),dr]);
    projected_r_high = r(end)*((f(refFreqIdx)./f(1)).^(1/betaML));
    r_high = projected_r_high + (dr - mod(projected_r_high,dr));


    % The interpolation is done over all possible values of rRequired whose
    % minimum value may span less than rTrue(1) and maximum value spans 
    % more than rTrue(end).
    %
    % 'ZBig' is the spectrogram with nan values at ranges outside of r for
    % interpolation pruposes.
    %
    % 'striationIdxStart' is the index of the striation that corresponds to
    % the first snapshot of the initially given spectrogram.
    %
    % There are some edge cases to consider
    if (r(1) < dr)
        rRequired = r_low:dr:r_high;
    
        % indices for comparison of the same striations (the first snapshot)
        [~,striationIdxEnd] = min(abs(rRequired-r(end)));
        [~,striationIdxStart] = min(abs(dr-r));
        IBig = [I(striationIdxStart:end,:); nan(numel(rRequired)-striationIdxEnd,size(I,2))];
    else        
        r1 = r_low:dr:r(1);
        r2 = (r(end)+dr):dr:r_high;
        rRequired = [r1, r, r2];
        [~,striationIdxStart] = min(abs(rRequired-r1(end)));
        [~,striationIdxEnd] = min(abs(rRequired-r2(1)));

        IBig = [nan(striationIdxStart,size(I,2)); I; nan(numel(rRequired)-striationIdxEnd+1,size(I,2))];
    end
    % Perform the interpolation 
    % IBig = abs(ZBig).^2;
    IInterp = getStriationIntensity(f,betaML,rRequired,IBig,refFreqIdx);
end