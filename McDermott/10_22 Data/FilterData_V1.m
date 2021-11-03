function smoothDat = FilterData_V1(t,datasets,freqThreshold)

if length(t) ~= length(datasets)
    error('Time and datasets must have the same length.')
end
if nargin == 2
    freqThreshold = 0.1;
end

L = length(t);
if mod(L,2)
    L = L-1;
end
Fs = round(L/(t(end) - t(1)));
nSets = size(datasets,2);

% Fourier-transformed data. Computed following MATLAB Documentation for
% fft() built-in function.
f = Fs*(0:(L/2))/L;
trans = repmat(struct(),nSets,1);
for n3 = 1:nSets
    trans(n3).fourier = fft(datasets(:,n3));     % Fourier transform.
    P2 = abs(trans(n3).fourier/L);
    P1 = P2(1:L/2+1);
    trans(n3).sss = P1;
    trans(n3).sss(2:end-1) = 2*P1(2:end-1);   % Single-sided spectrum.
    trans(n3).maxFreq = find(((trans(n3).sss - freqThreshold*max(trans(n3).sss))<0),1,'first');
end

smoothDat = zeros(size(datasets));

%{
figure();
for n3 = 1:nSets
    smoothDat(:,n3) = lowpass(datasets(:,n3),trans(n3).maxFreq,Fs);
    subplot(1,nSets,n3);
    plot(smoothDat(:,n3));
end
%}
%figure();
for n3 = 1:nSets
    smoothDat(:,n3) = movmean(datasets(:,n3),(1/trans(n3).maxFreq)/(1/Fs));
    %subplot(1,nSets,n3);
    %plot(movmean(datasets(:,n3),(1/trans(n3).maxFreq)/(1/Fs)));
end