function feat = getrecewpcdb5featCT(x,winsize,wininc,datawin,dispstatus)

%% This software implements wavelet packet energy feature extraction
% Please cite: Topçu, Ça?da?, et al. 
% "Assessment of Emotional Expressions after Full-Face Transplantation."
% Neural plasticity 2017 (2017).
%% Original Author: Cagdas Topcu, 2015 <Topcu.Cagdas@mayo.edu>

if nargin < 5
    if nargin < 4
        if nargin < 3
            if nargin < 2
                winsize = size(x,1);
            end
            wininc = winsize;
        end
        datawin = ones(winsize,1);
    end
    dispstatus = 0;
end

% winsize = 500;
% wininc = 250;
% datawin = hamming(winsize);
datasize = size(x,1);
Nsignals = size(x,2);
numwin = floor((datasize - winsize)/wininc)+1;

if dispstatus
    h = waitbar(0,'Computing EWPC features...');
end

st = 1;
en = winsize;

    global wname nMax;
    nMax = 7;
    % allocate memory
    if nMax == 5
        columbFeat = 32;
    elseif nMax == 7
        columbFeat = 128;
    else
        columbFeat = [];
    end
feat = zeros(numwin,columbFeat);
    wname = 'db5';
% disp(numwin)
for i = 1:numwin
    
    curwin = x(st:en,:);
    [T] = wpdec(curwin ,nMax ,wname) ;
    tn = leaves(T,'s');
    nbtn = length(tn);
    E  = zeros(1,nbtn);
    for k=1:nbtn
    C = wprcoef(T,tn(k));
    E(k) = sum(C(:).^2);
    end
    feat(i,:) = E;
   
    st = st + wininc;
    en = en + wininc;
% for testing
%     disp(i)
end

if dispstatus
    close(h)
end