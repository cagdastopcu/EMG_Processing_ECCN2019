
function [filterCoef] = FilterBPRawCT (cF1, cF2, Fs)

%% This software implements a bandpass butterworth filter
%  inputs:
%  cF1: First cutoff frequency
%  cF2: Second cutoff freq.
%  Fs: Sampling freq. in Hertz.
%  outputs:
%  filterCoef: Created filter coefficients, filterCoef.b and
%  filterCoef.a.
%% Usage
% Fs = 32000;
% cF1 = 65;
% cF2 = 110;
% [filterCoef] = FilterBPRaw (cF1, cF2, Fs);
% filteredSignal = filtfilt(filterCoef.b,filterCoef.a,rawSignal);

%% Original Author: Cagdas Topcu, 2018 <Topcu.Cagdas@mayo.edu>
%% 
% b  = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);
% Hd = dfilt.dffir(b);

N   = 4;     % Order after filtfilt and state that number in the publication

%% Calculate numerator b

[b,a] = butter(N/2, [cF1 cF2]/(Fs/2),'bandpass');

% Calculate the zpk values using the BUTTER function.
% [z,p,k] = butter(N/2, [cF1 cF2]/(Fs/2));
% 
% [sos_var,g] = zp2sos(z, p, k);
% Hd = dfilt.df2sos(sos_var, g);

filterCoef.a = a;

filterCoef.b = b;

end