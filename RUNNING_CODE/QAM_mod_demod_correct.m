close all;
clear all;

addpath('..');

sampleFreq = 4e6;
sampleTime = 1/sampleFreq;
fc = 1.1e6;
resampleFactor = 32;

data_source = [1 5 3 2 8 9 2 4 5 6 3 9 0 10 3 11 15 2 0 2 4 3 12 13 15 3 2 7 10 2 5 7 8 1 10 11 2 15 13];
modulated_data = qammod(data_source, 16);

ofdm_signal = ifft(modulated_data);
ofdm_signal = resample(ofdm_signal,resampleFactor,1);   %resample to have more samples per symbol

tc = [0:sampleTime:(length(ofdm_signal)-1)*sampleTime];
mixedSignal = real(ofdm_signal .* exp(i*2*pi*fc*tc));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DEMODULATION
%   QPSK
%       - 128 carriers
%       - 13 bits for cyclic prefix
%       - 4 times upsampled

origSignal = mixedSignal .* exp(-i*2*pi*fc*tc);

figure;
plot(abs(fft(origSignal,1024)));

% Zeroing fft spectrum looks good. IT'S RUNNING ALSO WITHOUT IT! WHY?
filtSignal = fft(origSignal);
filtSignal(150:750) = 0;
filtSignal = 2*ifft(filtSignal);

figure;
plot(abs(fft(filtSignal,1024)));

filtSignal = resample(filtSignal,1,resampleFactor);
fftSignal = fft(filtSignal);
demodSignal = qamdemod(fftSignal,16);

figure;
stem(data_source,'-o'); hold on;
stem(demodSignal,'r--x'); hold off;


