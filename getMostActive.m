function [mostActiveSTD] = getMostActive(dataAllChannels)

%% Band-pass filter parameters
nBP = 4;
Fs = 2000;
cF1 = 10;
cF2 = 100;

[N, M] = size(dataAllChannels);

[bLF,aLF] = butter(nBP/2,[cF1, cF2]/(Fs/2));

%% Spectrogram parameters
paramsLF.tapers = [1 1 1]; %freq res x time res (s) product at end ->
% controls how many tapers used in the analysis 2*1st*2nd - 3rd
paramsLF.Fs=2000; %sampling frequency
paramsLF.fpass=[10 100]; %freq band --> Tory chnage from 62 118 to 1 50
paramsLF.pad=-1; %DSP technique to smooth data before you do the fft
%-1 is switched off
paramsLF.err=[]; %useful for stats, put alpha value of test and which test you want
paramsLF.trialave = 0;   %if you arrange in matrices, will avg over columns


for m = 1:M

LFtemp(:,m) = filtfilt(bLF,aLF,repmat(dataAllChannels(:,m),3,1));

LF(:,m) = LFtemp(N+1:2*N,m);

[LFS(:,:,m),t1(m,:),f1(m,:)] = mtspecgramc(LF(:,m)-mean(LF(:,m)),[0.5 .05],paramsLF);

end


[KLF LLF MLF] = size(LFS);

meanLFS = mean(LFS(:,:,:),2);

meanLFS=reshape(meanLFS,KLF,MLF);

meanLFS = mean(meanLFS);

for i = 1:KLF
    
    stdLFS(i,:,:) = std(LFS(i,:,:),0,2);

end



[maxLF,mostActiveMean] = max(meanLFS)

[maxstdLF,mostActiveSTD] = max(meanLFS)

if ~isequal(mostActiveMean,mostActiveSTD)
   error('Most acrive channels are different for std andmean values')
end

mostActive.mostActiveMean = mostActiveMean;

mostActive.mostActiveSTD = mostActiveSTD;
end