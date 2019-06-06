
function [amplitudeHT] = getAmplitudeHTCT(xSignal)

%% This software calculates Hilbert Transform Amplitude
%  ---inputs---
%  xSignal: Filtered signal
%  ---outputs---
%  amplitudeHT: Hilbert Transform Amplitude or envelope
%% Usage
% % [amplitudeHT] = getAmplitudeHT(xSignal);

%% Original Author: Cagdas Topcu, 2018 <Topcu.Cagdas@mayo.edu>
%% 
%% fft based Hilbert Transform and amplitude calculation
amplitudeHT = sqrt(imag(hilbert(xSignal)).^2+xSignal.^2); 

end
