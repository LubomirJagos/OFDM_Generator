close all;
clear all;

addpath('..');

sampleFreq = 4e6;
sampleTime = 1/sampleFreq;
fc = 1.1e6;
resampleFactor = 32;

data_source = [3 0 1 0 2 3 3 2 1 0 2 3 1 2 0 3 1 2 1 0 3 1 2 3 1 2 3 0 1 2 0 1];
modulated_data = pskmod(data_source, 4);

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
% filtSignal = fft(origSignal);
% filtSignal(100:800) = 0;
% filtSignal = 2*ifft(filtSignal);

filtSignal = origSignal;

filtSignal = resample(filtSignal,1,resampleFactor);
fftSignal = fft(filtSignal);
demodSignal = pskdemod(fftSignal,4);

figure;
stem(data_source,'-o'); hold on;
stem(demodSignal,'r--x'); hold off;


