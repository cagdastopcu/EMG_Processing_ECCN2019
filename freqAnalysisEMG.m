clear all;
close all;
% clc,
format long g
% C:\Arm paper\Data to process\pelinleymanData\pelinleyman1\pelinleyman28012015
% C:\Arm paper\Data to process\cihantopal
% H01 C:\Arm paper\Data to process\UpperLimbEMG19122014\cagdas27032014
% H02 C:\Arm paper\Data to process\UpperLimbEMG19122014\emredogan08042014
% H03 C:\Arm paper\Data to process\UpperLimbEMG19122014\fezacarlak09042014
% H04 C:\Arm paper\Data to process\UpperLimbEMG19122014\hakanaktas 02042014
% H05 C:\Arm paper\Data to process\UpperLimbEMG19122014\nurettin 24032014
% H06 C:\Arm paper\Data to process\UpperLimbEMG19122014\omercolak10042014
% H07 C:\Arm paper\Data to process\UpperLimbEMG19122014\ovuncpolat14032014
% H08 C:\Arm paper\Data to process\UpperLimbEMG19122014\refiksever 04042014
% H09 C:\Arm paper\Data to process\UpperLimbEMG19122014\seckinsen 28032014
% H10 C:\Arm paper\Data to process\UpperLimbEMG19122014\umitulusar26032014

folderpath = 'C:\Arm paper\Data to process\UpperLimbEMG19122014\umitulusar26032014';
folderpath = fullfile(folderpath, '**');    % What is the meaning of "/**/" ???
filelist   = dir(folderpath);
filenames_cell       = {filelist.name};
filenames_cell       = filenames_cell(~strncmp(filenames_cell, '.', 1))   % No files starting with '.'
% load('filenames_cell.mat')

% MovNo = 16;

for MovNo = 1:22

indvNames = {'ArmTransplant' 'Arm Replant' 'H01' 'H02' 'H03' ...
    'H04' 'H05' 'H06' 'H07' 'H08' 'H09' 'H10'};

dataDir = strcat('C:\Arm paper\Data to process\UpperLimbEMG19122014\umitulusar26032014', '\', filenames_cell{MovNo})

% dataDir = strcat('C:\Arm paper\Data to process\UpperLimbEMG19122014\cagdas27032014', '\', filenames_cell{MovNo})%H01

% dataDir = strcat('C:\Arm paper\Data to process\UpperLimbEMG19122014\emredogan08042014', '\', filenames_cell{MovNo})%H02


load(dataDir)

tic
for i = 1:24
% dataAllChannels(:,i) = data(datastart(i)+11000:dataend(i)-41000);
% chan2(:,i) = data(datastart(i)+31000:dataend(i)-21000);
dataAllChannels(:,i) = data(datastart(i)+4000:dataend(i));
end

mostActive= getMostActive(dataAllChannels);
dataMostActive = dataAllChannels(:,mostActive(19));
[onsetTime,feat] = getOnsetTime(dataMostActive');

clear dataAllChannels;
numPoints = 66000;
extractedTime = 80000-(datastart(1)+3800+onsetTime-1)-numPoints;
dataAllChannels = zeros(numPoints,24);

for i = 1:24
% dataAllChannels(:,i) = data(datastart(i)+11000:dataend(i)-41000);
% chan2(:,i) = data(datastart(i)+31000:dataend(i)-21000);
dataAllChannels(:,i) = data(datastart(i)+3800+onsetTime:dataend(i)-extractedTime);
end

[MData NData] = size(dataAllChannels);

nBP = 4;
Fs = 2000;
cF1 = 10;
cF2 = 50;
[bVLF,aVLF] = butter(nBP/2,[cF1, cF2]/(Fs/2));
cF1 = 50;
cF2 = 100;
[bLF,aLF] = butter(nBP/2,[cF1, cF2]/(Fs/2));
cF1 = 100;
cF2 = 150;
[bMF,aMF] = butter(nBP/2,[cF1, cF2]/(Fs/2));

cF1 = 150;
cF2 = 250;
[bHF,aHF] = butter(nBP/2,[cF1, cF2]/(Fs/2));

cF1 = 250;
cF2 = 450;
[bVHF,aVHF] = butter(nBP/2,[cF1, cF2]/(Fs/2));
%% spectrogram parameters
paramsVLF.tapers = [4 0.25 1]; %freq res x time res (s) product at end ->
% controls how many tapers used in the analysis 2*1st*2nd - 3rd
paramsVLF.Fs=2000; %sampling frequency
paramsVLF.fpass=[10 50]; %freq band --> Tory chnage from 62 118 to 1 50
paramsVLF.pad=-1; %DSP technique to smooth data before you do the fft
%-1 is switched off
paramsVLF.err=[]; %useful for stats, put alpha value of test and which test you want
paramsVLF.trialave = 0;   %if you arrange in matrices, will avg over columns

paramsLF = paramsVLF;
paramsLF.fpass = [50 100];
paramsMF = paramsVLF;
paramsMF.fpass = [100 150];
paramsHF = paramsVLF;
paramsHF.fpass = [150 250];
paramsVHF = paramsVLF;
paramsVHF.fpass = [250 450];

% ltht2(:,m) = filtfilt(theta,1,repmat(BIPOLAR{1,x(i)}(:,NR(m)),3,1));
%         LTHT2(:,m) = ltht2(1501:3000,m);
% dataAllChannelsF = filtfilt(b,a,dataAllChannels);% sifir faz

% %% Band durduran centik filtresi
% 
% [b,a] = butter(1,[49 51]/(Fs/2), 'stop');
% dataAllChannelsF = filtfilt(b,a,dataAllChannelsF);% sifir faz
% [b,a] = butter(1,[249 251]/(Fs/2), 'stop');
% dataAllChannelsF = filtfilt(b,a,dataAllChannelsF);% sifir faz
% [b,a] = butter(1,[349 351]/(Fs/2), 'stop');
% dataAllChannelsF = filtfilt(b,a,dataAllChannelsF);% sifir faz
% [b,a] = butter(1,[449 451]/(Fs/2), 'stop');
% dataAllChannelsF = filtfilt(b,a,dataAllChannelsF);% sifir faz

%% PSD


%     window = 400;

%     noverlap = round(window/2);
% 
%     freqs = 2:2:500;
% 
%     fs = 2000;
%     
%     [PSD,f] = pwelch(dataAllChannels-mean(dataAllChannels),...
%         window,noverlap,freqs,fs);
% 
%     plot(freqs,10*log10(PSD),'r','linewidth',2)
%     hold on
%     [PSD,f] = pwelch(dataAllChannelsF-mean(dataAllChannelsF),...
%         window,noverlap,freqs,fs);
% 
%     plot(freqs,10*log10(PSD),'k','linewidth',2)


dataTempStart = MData+1;
dataTempStop = 2*MData;

%% spectrogram
for m = 1:NData
    
    %% VLF: 5-50 Hz LF: 50-100 MF: 100-150
    %HF: 150-250 VHF:250-400 Hz
    
   [b,a] = butter(1,[49 51]/(Fs/2), 'stop');
   dataAllChannelsF1Temp(:,m) = filtfilt(b,a,repmat(dataAllChannels(:,m),3,1));% sifir faz
   dataAllChannelsF1(:,m) = dataAllChannelsF1Temp(dataTempStart:dataTempStop,m);
   [b,a] = butter(1,[249 251]/(Fs/2), 'stop');
   dataAllChannelsF2Temp(:,m) = filtfilt(b,a,repmat(dataAllChannelsF1(:,m),3,1));
   dataAllChannelsF2(:,m) = dataAllChannelsF2Temp(dataTempStart:dataTempStop,m);
   [b,a] = butter(1,[349 351]/(Fs/2), 'stop');
   dataAllChannelsF3Temp(:,m) = filtfilt(b,a,repmat(dataAllChannelsF2(:,m),3,1));
   dataAllChannelsF3(:,m) = dataAllChannelsF3Temp(dataTempStart:dataTempStop,m);
   
   VLFtemp(:,m) = filtfilt(bVLF,aVLF,repmat(dataAllChannelsF3(:,m),3,1));
   VLF(:,m) = VLFtemp(dataTempStart:dataTempStop,m);
   [VLFS(:,:,m),t1(m,:),f1(m,:)] = mtspecgramc(VLF(:,m)-mean(VLF(:,m)),[0.5 .05],paramsVLF);
   
   LFtemp(:,m) = filtfilt(bLF,aLF,repmat(dataAllChannelsF3(:,m),3,1));
   LF(:,m) = LFtemp(dataTempStart:dataTempStop,m);
   [LFS(:,:,m),t2(m,:),f2(m,:)] = mtspecgramc(LF(:,m)-mean(LF(:,m)),[0.5 .05],paramsLF);
   
   MFtemp(:,m) = filtfilt(bMF,aMF,repmat(dataAllChannelsF3(:,m),3,1));
   MF(:,m) = MFtemp(dataTempStart:dataTempStop,m);
   [MFS(:,:,m),t3(m,:),f3(m,:)] = mtspecgramc(MF(:,m)-mean(MF(:,m)),[0.5 .05],paramsMF);
   
   HFtemp(:,m) = filtfilt(bHF,aHF,repmat(dataAllChannelsF3(:,m),3,1));
   HF(:,m) = HFtemp(dataTempStart:dataTempStop,m);
   [HFS(:,:,m),t4(m,:),f4(m,:)] = mtspecgramc(HF(:,m)-mean(HF(:,m)),[0.5 .05],paramsHF);
   
   VHFtemp(:,m) = filtfilt(bVHF,aVHF,repmat(dataAllChannelsF3(:,m),3,1));
   VHF(:,m) = VHFtemp(dataTempStart:dataTempStop,m);
   [VHFS(:,:,m),t5(m,:),f5(m,:)] = mtspecgramc(VHF(:,m)-mean(VHF(:,m)),[0.5 .05],paramsVHF);
end

[K L M] = size(VLFS);
% VLFS = zscore(VLFS);
meanVLFS = mean(VLFS(:,:,:),2);

meanVLFS=reshape(meanVLFS,K,M);
meanVLFS =  mean(meanVLFS,2);

for i = 1:K
    
    stdVLFS(i,:,:) = std(VLFS(i,:,:),0,2);

end
meanstdVLFS = mean(stdVLFS,2);

[K L M] = size(LFS);
% VLFS = zscore(VLFS);
meanLFS = mean(LFS(:,:,:),2);

meanLFS=reshape(meanLFS,K,M);
meanLFS =  mean(meanLFS,2);
for i = 1:K
    

    stdLFS(i,:,:) = std(LFS(i,:,:),0,2);


end
meanstdLFS = mean(stdLFS,2);

[K L M] = size(MFS);
% VLFS = zscore(VLFS);
meanMFS = mean(MFS(:,:,:),2);

meanMFS=reshape(meanMFS,K,M);
meanMFS =  mean(meanMFS,2);
for i = 1:K
    

    stdMFS(i,:,:) = std(MFS(i,:,:),0,2);


end
meanstdMFS = mean(stdMFS,2);


[K L M] = size(HFS);
meanHFS = mean(HFS(:,:,:),2);

meanHFS=reshape(meanHFS,K,M);
meanHFS =  mean(meanHFS,2);

for i = 1:K
    

    stdHFS(i,:,:) = std(HFS(i,:,:),0,2);

end
meanstdHFS = mean(stdHFS,2);


[K L M] = size(VHFS);
meanVHFS = mean(VHFS(:,:,:),2);

meanVHFS=reshape(meanVHFS,K,M);
meanVHFS =  mean(meanVHFS,2);

for i = 1:K
    

    stdVHFS(i,:,:) = std(VHFS(i,:,:),0,2);

end
meanstdVHFS = mean(stdVHFS,2);
%%
% figure('units','normalized','outerposition',[0 0 1 1])

maxLFS = max(meanLFS(:));

maxstdLFS = max(meanstdLFS(:));

maxVLFS = max(meanVLFS(:));

maxstdVLFS = max(meanstdVLFS(:));

meanNormS = mean([maxLFS maxVLFS]);
meanNormstdS = mean([maxstdLFS maxstdVLFS]);

FigH = figure;

[linehandlesVLFS, patchhandles] = errorarea(t1(1,1:end),(meanVLFS(:)/meanNormS),(meanstdVLFS(:))/sqrt(14)/meanNormstdS);
hold on
[linehandlesLFS, patchhandles] = errorarea(t2(1,1:end),(meanLFS(:)/meanNormS),(meanstdLFS(:))/sqrt(14)/meanNormstdS);

[linehandlesMFS, patchhandles] = errorarea(t3(1,1:end),(meanMFS(:)/meanNormS),(meanstdMFS(:))/sqrt(14)/meanNormstdS);

[linehandlesHFS, patchhandles] = errorarea(t4(1,1:end),(meanHFS(:)/meanNormS),(meanstdHFS(:))/sqrt(14)/meanNormstdS);

[linehandlesVHFS, patchhandles] = errorarea(t5(1,1:end),(meanVHFS(:)/meanNormS),(meanstdVHFS(:))/sqrt(14)/meanNormstdS);

lgnd = legend([linehandlesVLFS linehandlesLFS linehandlesMFS linehandlesHFS linehandlesVHFS], 'VLF: 10-50Hz', 'LF: 50-100Hz', 'MF: 100-150Hz', 'HF: 150-250Hz', 'VHF: 250-450Hz');
lgnd.FontSize = 18;
axis([0 34 0 1.8])
xlabel('time (s)')
ylabel('normalized power spectrum')
set(gca,'fontsize',32)
 set(gca,'box','off')
set(gca,'linewidth',10)
legend boxoff     

% 
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1) -0.14;
bottom = outerpos(2) + ti(2) - 0.1;
ax_width = outerpos(3) - ti(1) - ti(3) + 0.13;
ax_height = outerpos(4) - ti(2) - ti(4) + 0.08;
ax.Position = [left bottom ax_width ax_height];
 h = findobj(gca,'Type','line');
 set(FigH,'WindowState','maximized')
% %  %% it works
% % set(FigH,'Units','Inches');
% % pos = get(FigH,'Position');
% % set(FigH,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % print(FigH,'filename','-dpdf','-r600')
% % %% %%%%%%%%%%%%%%%%%

%  print('ScreenSizeFigure','-dpng','-r0')
%  res = 300;
%  size = [255*ax_width 255*ax_height];
% set(gcf,'paperunits','inches','paperposition',[left bottom size/res]);
% print('resized.tiff','-dtiff',['-r' num2str(res)]);
%  WindowAPI(FigH, 'full');  % fill the current monitor
%  frame_h = get(handle(gcf),'JavaFrame');
% set(frame_h,'Maximized',1);
toc
% 
% pos = get(gca,'Position');
% set(ax,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(ax,'fistTest','-dpdf','-r0')

%% Save 

MovNumInCell = {'01' '02' '03' '04' '05'...
        '06' '07' '08' '09' '10' ...
        '11' '12' '13' '14' '15' '16'...
        '17' '18' '19' '20' ...
        '21' '22'};
    
    mkdir('H10')
    nameFile = [pwd '\H10\' MovNumInCell{MovNo} '-' indvNames{12} '-' filenames_cell{MovNo} '.png'];

print(nameFile,'-dpng','-r0')
clear dataAllChannels dataAllChannelsF1 dataAllChannelsF1Temp...
    dataAllChannelsF2 dataAllChannelsF2Temp...
    dataAllChannelsF3 dataAllChannelsF3Temp;
close all;
end