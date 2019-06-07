function fuzzyEn = FuzzyEn( dim, r, data, tau )
% FuzzyEN Fuzzy Entropy
%   calculates the fuzzy entropy of a given time series data

%   dim     : embedded dimension
%   r       : tolerance (typically 0.2 * std)
%   data    : time-series data
%   tau     : delay time for downsampling (user can omit this, in which case
%             the default value is 1)
%
%---------------------------------------------------------------------
% implemented by Cagdas Topcu, Modified from Kijoon Lee's (2012) Sample Entropy function
% cagdastopcu@gmail.com
% 2015
% Cite: Topçu, Çağdaş, et al. 
% "Recovery of facial expressions using functional electrical stimulation after full-face transplantation." 
% Journal of neuroengineering and rehabilitation 15.1 (2018): 15.
%---------------------------------------------------------------------
% data =randn(1,1000);
% r = std(data);
% tau=1;
% dim=2;
if nargin < 4, tau = 1; end
if tau > 1, data = downsample(data, tau); end

N = length(data);
correl = zeros(1,2);
dataMat = zeros(dim+1,N-dim);
for i = 1:dim+1
    dataMat(i,:) = data(i:N-dim+i-1);
end

for m = dim:dim+1
    count = zeros(1,N-dim);
    tempMat = dataMat(1:m,:);
    
    for i = 1:N-m
        % calculate Chebyshev distance, excluding self-matching case
        dist = max(abs(tempMat(:,i+1:N-dim) - repmat(tempMat(:,i),1,N-dim-i)));
		
        D =  exp(-(dist.^2)/(0.25*r));
        
        count(i) = sum(D)/(N-dim);
    end
    
    correl(m-dim+1) = sum(count)/(N-dim);
end

fuzzyEn = log(correl(1)/correl(2));

end