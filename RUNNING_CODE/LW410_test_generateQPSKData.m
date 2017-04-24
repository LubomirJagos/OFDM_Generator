close all;
clear all;

addpath('..');

lw410 = LW410Interface();
lw410.wave_data

data_source = [2 2 1 2 0 3 1 3 3 3 2 0 1 2 0 1 3 1 0 2 3 1 0 3 1 2 0 3 1 0 2 3 0 1 2 0 3 1 2 1 0 1 2 0 3 2 1 0 2 3 0 1 0 2 2 3 1 1 1 2 1 2 3 1];
modulated_data = pskmod(data_source, 4);

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




