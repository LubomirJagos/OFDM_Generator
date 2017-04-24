close all;
clear all;

lw410 = LW410Interface();
%data_source = [3 3 0 2 1 3 3 2 0 1 2 3 0 1 2 0 3 2 1 0 2 3 0 1 0 1 2 0 1 2 3 1 2 2 1 0 1 2 3 1 1 0 2 2 3 1 1 0 2 0 3 3 2 0 1 3 1 2 3 0 1 2 1 3];
data_source = [3 2 1 3 3 2 0 1 2 1 2 0 3 2 1 0 2 1 0 1 2 0 1 2 3 2 1 0 1 2 3 1 1 0 2 2 3 1 1 0 2 0 3 3 2 0 1 3 1 2 3 0 1 2 1 3];

modulated_data = pskmod(data_source, 4);
ofdm_signal = ifft(modulated_data);
ofdm_signal = resample(ofdm_signal,32,1);   %resample to have more samples per symbol

fc = 0.7e6;
tc = [0:lw410.sampletime:(length(ofdm_signal)-1)*lw410.sampletime];
%mixedSignal = real(ofdm_signal) .* cos(2*pi*fc*tc) + imag(ofdm_signal)*sin(2*pi*fc*tc));
mixedSignal = real(ofdm_signal .* exp(i*2*pi*fc*tc));

figure;
xAxis = linspace(0, 1/lw410.sampletime, 1024);
subplot(211); plot(mixedSignal);
subplot(212); plot(xAxis, abs(fft(mixedSignal, 1024))); title('OFDM signal spectrum'); xlabel('f[kHz]');







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DEMODULATION
%   QPSK
%       - 128 carriers
%       - 13 bits for cyclic prefix
%       - 4 times upsampled

resampleFactor = 32;
hantekSampleFreq = 4e6;
hantekSampleTime = 1/hantekSampleFreq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ERROR is somewhere here in IQ demodulation
%

fc = 0.7e6;
tc = [0:hantekSampleTime:(length(ofdm_signal)-1)*hantekSampleTime];
origSignal = mixedSignal.*exp(i*2*pi*fc*tc);
%origSignal = mixedSignal .* (cos(2*pi*fc*tc) + i*sin(2*pi*fc*tc));

% Filter is not running good, probably is slightly changing phase
% characteristic or something else.
% filtFc = 200e3;
% [b,a] = butter(20,filtFc/(hantekSampleFreq/2));
% filtSignal = filter(b,a,real(origSignal)) + i*filter(b,a,imag(origSignal));

% Zeroing fft spectrum looks good.
filtSignal = fft(origSignal);
filtSignal(500:1000) = 0;
figure;
stem(filtSignal);
filtSignal = 2*ifft(filtSignal);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DEBUGGING%%%%%%%%%%%%%%%%%%%%%%%%%%
% origSignal = ofdm_signal;
% filtSignal = ofdm_signal;
filtSignal = conj(filtSignal);
% resampleFactor = 32;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
plot(abs(filtSignal), '-'); hold on;
plot(abs(ofdm_signal),'g-'); hold off;
figure;
plot(angle(filtSignal), '-'); hold on;
plot(angle(ofdm_signal),'g-'); hold off;



figure;
subplot(311);
plot(real(origSignal), 'b-x'); hold on;
plot(imag(origSignal), 'r-x'); hold off;
title('I and Q part of ofdm signal');
subplot(312); plot(abs(fft(origSignal, 1024)), '-x');
title('IQ signal before filtering spectrum');
subplot(313); plot(abs(fft(filtSignal, 1024)), '-x');
title('Filtered signal spectrum');




figure;
plot(real(ofdm_signal), '-'); hold on;
plot(real(filtSignal), 'r-'); hold off;











filtSignal = resample(filtSignal,1,resampleFactor);
fftSignal = fft(filtSignal);
demodSignal = pskdemod(fftSignal,4);

refSignal = [3 3 0 2 1 3 3 2 0 1 2 3 0 1 2 0 3 2 1 0 2 3 0 1 0 1 2 0 1 2 3 1 2 2 1 0 1 2 3 1 1 0 2 2 3 1 1 0 2 0 3 3 2 0 1 3 1 2 3 0 1 2 1 3];
figure;
stem(data_source,'-x'); hold on;
stem(demodSignal,'r--x'); hold off;


