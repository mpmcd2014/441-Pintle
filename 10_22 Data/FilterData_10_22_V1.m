filename = 'Data\Pintle\Pintle_10-22_125-1.csv';
freqThreshold = 0.1;   % Include frequencies that contribute more than 10% of amplitude.

rawData = csvread(filename,1,0);

raw = struct();
raw.t = rawData(:,1);
raw.data.P = rawData(:,7);
raw.data.m = sum(rawData(:,16:19),2);
clear rawData

L = length(raw.t);
Fs = round(L/(raw.t(end) - raw.t(1)));
sets = fieldnames(raw.data);

% Fourier-transformed data. Computed following MATLAB Documentation for
% fft() built-in function.
f = Fs*(0:(L/2))/L;
trans = struct();
for n3 = 1:length(sets)
    trans.(sets{n3}) = struct();
    trans.(sets{n3}).fourier = fft(raw.data.(sets{n3}));     % Fourier transform.
    P2 = abs(trans.(sets{n3}).fourier/L);
    P1 = P2(1:L/2+1);
    trans.(sets{n3}).sss = P1;
    trans.(sets{n3}).sss(2:end-1) = 2*P1(2:end-1);   % Single-sided spectrum.
    trans.(sets{n3}).maxFreq = find(((trans.(sets{n3}).sss - freqThreshold*max(trans.(sets{n3}).sss))<0),1,'first');
end

figure();
filt = struct();
for n3 = 1:length(sets)
    filt.(sets{n3}) = lowpass(raw.data.(sets{n3}),trans.(sets{n3}).maxFreq,Fs);
    subplot(1,length(sets),n3);
    plot(filt.(sets{n3}));
end

figure();
for n3 = 1:length(sets)
    subplot(1,length(sets),n3);
    plot(movmean(raw.data.(sets{n3}),(1/trans.(sets{n3}).maxFreq)/(1/Fs)));
end

tmp = FilterData_V1(raw.t,[raw.data.P,raw.data.m],freqThreshold);