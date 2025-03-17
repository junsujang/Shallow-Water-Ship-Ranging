% Junsu Jang (junsu.jang94@gmail.com)
% 2025/01/20
function rObjFn = computeRangeDFTRadonObjFn(I,rMax,broadbandF,cMin,cMax,thetas,stepSizeR,stepSizeF)
% Compute the slope of the striation in the Radon transform domain
% Yao et al. (2021)

[NFFT_F,NFFT_R] = size(I);

% demean the intensity to remove the DC component
I = I - mean(I,'all');          
% Compute the FFT bins in delay-wavenumber domain
Fs_R = 1/stepSizeR;
Fs_F = 1/stepSizeF;
dkr = Fs_R/NFFT_R; 
dkf = Fs_F/NFFT_F;

% For FFT shift
kr = (0:(NFFT_R-1))*dkr;
kr(kr >= Fs_R/2) = kr(kr >= Fs_R/2) - Fs_R;
kr = 2*pi*fftshift(kr);
kf = (0:(NFFT_F-1))*dkf;
kf(kf >= Fs_F/2) = kf(kf >= Fs_F/2) - Fs_F;
kf = 2*pi*fftshift(kf);

% Filter out the signal unnecessary in the delay-wavenumber domain
% This is rather essential because a little extra noise changes the result
% significantly.
kf_max = 2*pi*rMax*(1/cMin-1/cMax);
kr_max = 2*pi*broadbandF(end)*(1/cMin-1/cMax)/2; % 3 for 15 min segment
krValid = abs(kr) <= kr_max;
kfValid = abs(kf) <= kf_max;

% Take the 2D-DFT --> obtain X
I2 = fft2(I);
I2DT = abs(fftshift(I2));
% Filter out the window with slope-information
X = I2DT(kfValid,krValid);

[R,xp] = radon(X,thetas);
% Identify the slope of the 2DFT: Only consider the energy that passes 
% through the center of the 2D-DFT (the DC component)
rObjFn = R((xp == 0),:);

end