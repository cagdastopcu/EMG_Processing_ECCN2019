function [onsetTime] = getOnsetTime(mostActiveChannelData)

nBP = 4;
Fs = 2000;
cF1 = 10;
cF2 = 100;

x = mostActiveChannelData;

[Mx Nx] = size(x);

[b,a] = butter(nBP/2,[cF1, cF2]/(Fs/2));

    xFtemp = filtfilt(b,a,repmat(x',3,1));
    xF = xFtemp(Mx+1:2*Mx);

dim = 2;
globaltolerance = std(x);
tau = 1;
winsize = 50;
wininc = winsize/50;
datawin = hamming(winsize);
dispstatus = 0;


feat = getfuzzyenfeat(x',globaltolerance,winsize,wininc,datawin,dispstatus);

dfuzzyEn = abs(diff(feat));
dfuzzyEn(dfuzzyEn<=0.01) = 0;
onset1 = find(dfuzzyEn>0);

onsetTime = onset1(1);

end