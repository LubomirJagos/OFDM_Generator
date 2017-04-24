close all;
clear all;

addpath('..');

lw410 = LW410Interface();
lw410.wave_data

data_source = [1 5 3 2 8 9 2 4 5 6 3 9 0 10 3 11 15 2 0 2 4 3 12 13 15 3 2 7 10 2 5 7 8 1 10 11 2 15 13];
modulated_data = qammod(data_source, 16);

% data_source = [2 3 7 1 2 5 4 7 5 6 1 0 1 2 5 0 4 3 0 5 5 3 2 6 5 4 0 5 6 3 2 5 1 0 6 7 1 7 3 7 5 5 1 0 2 0 6 1 2 7 3 4 1];
% modulated_data = qammod(data_source, 8);

ofdm_signal = ifft(modulated_data);
ofdm_signal = resample(ofdm_signal,32,1);   %resample to have more samples per symbol

fc = 0.7e6;
tc = [0:lw410.sampletime:(length(ofdm_signal)-1)*lw410.sampletime];
mixedSignal = ofdm_signal .* exp(i*2*pi*fc*tc);
realMixedSignal = real(mixedSignal);
imagMixedSignal = imag(mixedSignal);

figure;
plot(realMixedSignal); hold on;
plot(imagMixedSignal, 'r');
plot(realMixedSignal+imagMixedSignal, 'g'); hold off;

xAxis = linspace(0, 1/lw410.sampletime, 1024);
figure;
subplot(211); plot(ofdm_signal);
subplot(212); plot(xAxis, abs(fft(ofdm_signal, 1024))); title('OFDM signal spectrum'); xlabel('f[kHz]');
figure;
subplot(211); plot(mixedSignal);
subplot(212); plot(xAxis, abs(fft(mixedSignal, 1024))); title('OFDM signal spectrum'); xlabel('f[kHz]');

lw410.wave_data(realMixedSignal,1);
%lw410.wave_data(realMixedSignal+imagMixedSignal,1);




