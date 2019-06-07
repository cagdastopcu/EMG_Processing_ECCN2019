%
% GETFUZZYENFEAT Gets the Fuzzy Entropy feature.
%
% feat = getfuzzyenfeat(x,globaltolerance,winsize,wininc,datawin,dispstatus)

% Inputs
%    x: 		columns of signals
%    winsize:	window size (length of x)
%    wininc:	spacing of the windows (winsize)
%    datawin:   window for data (e.g. Hamming, default rectangular)
%               must have dimensions of (winsize,1)
%    dispstatus:zero for no waitbar (default)

function feat = getfuzzyenfeat(x,globaltolerance,winsize,wininc,datawin,dispstatus)

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

datasize = size(x,1);
Nsignals = size(x,2);
numwin = floor((datasize - winsize)/wininc)+1;

% allocate memory
feat = zeros(numwin,Nsignals);

if dispstatus
    h = waitbar(0,'Computing HFD features...');
end

st = 1;
en = winsize;

%    if dispstatus
%        waitbar(i/numwin);
%    end%alttaki for dongusunun icine gomulebilir.

for i = 1:numwin

   curwin = x(st:en,:).*repmat(datawin,1,Nsignals);
%       feat(i,:) = sqrt(mean(curwin.^2));

   
   dim = 2;
   r = globaltolerance;
%    r = 0.2 * std(curwin);
   tau = 1;

   feat(i,:) = FuzzyEn( dim, r, curwin, tau );
   
   st = st + wininc;
   en = en + wininc;
end

if dispstatus
    close(h)
end
