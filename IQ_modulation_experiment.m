% OFDM generation and demodulation experiment 1

close all;
clear all;

addpath('LJ_SignalFunctions');          %auxiliary functions

fs = 4e6;
ts = 1/fs;
t = [0:ts:250e-6];
fftPoints = 4096;

specX = linspace(0, fs, fftPoints);
signal = 1.2*exp(i*2*pi*80e3*t) + 1.2*exp(i*2*pi*400e3*t);

% figure;
% plot(real(signal), 'b-x'); hold on;
% plot(imag(signal), 'r-x'); hold off;
% title('Original signal');
% legend('I','Q')

fc = 1.4e6;
carrier = exp(i*2*pi*fc*t);

% figure;
% plot(real(carrier), 'b-x'); hold on;
% plot(imag(carrier), 'r-x'); hold off;
% title('Original signal');
% legend('I','Q')

modSignal = signal .* carrier;

% figure;
% plot(real(modSignal), 'b-x'); hold on;
% plot(imag(modSignal), 'r-x'); hold off;
% title('IQ signal');

figure;
subplot(311); stem(specX, abs(fft(signal, fftPoints)));
subplot(312); stem(specX, abs(fft(carrier, fftPoints)));
subplot(313); stem(specX, abs(fft(modSignal, fftPoints)));
title('Signals spectrum');

modSignal2 = real(modSignal) + imag(modSignal);
figure;
plot(specX, abs(fft(modSignal2, fftPoints)), '-x');
title('IQ UPMIXED signal spectrum');

%this is downmixing, there has to be low-pass filter, because
%there are 4 images of signal because previous mixing
modSignal2 = modSignal2 .* (cos(2*pi*fc*t) - i*sin(2*pi*fc*t));
modSignal2 = real(modSignal2) + imag(modSignal2);

figure;
plot(specX, abs(fft(modSignal2, fftPoints)), '-x');
title('IQ DOWNMIXED signal spectrum');


lw410 = LW410Interface();
resampleFactor = ceil(ts / lw410.sampletime);
lwFs = 1 / lw410.sampletime;
modSignal2 = resample(modSignal2, resampleFactor, 1);

figure;
subplot(411); plot(modSignal2, '-x');
title('Modulated signal in time');
subplot(412); plot(linspace(0,lwFs,fftPoints), abs(fft(modSignal2, fftPoints)));
title('Signal spectrum');

% filtN = 20;                    % Filter order
% filtF = [0 0.125 0.25 1];         % Frequency band edges
% filtA = [1  1  0 0];           % Desired amplitudes
% filtB = firpm(filtN,filtF,filtA);
% filtSignal = filter(filtB, filtB, modSignal2);

filtFc = 500e3;
[b,a] = butter(6,filtFc/(fs/2));
filtSignal = filter(b,a,modSignal2);

subplot(413); freqz(b,a);
title('Butterworth filter frequency response');

subplot(414); plot(linspace(0,lwFs,fftPoints), abs(fft(filtSignal, fftPoints)));
title('Filtered signal spectrum');




%lw410.wave_data(modSignal2, 1);



