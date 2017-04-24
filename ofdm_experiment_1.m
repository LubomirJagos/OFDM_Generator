% OFDM generation and demodulation experiment 1

close all;
clear all;

addpath('LJ_SignalFunctions');          %auxiliary functions

M = 4
nLength = 64
alphabet = [0:M-1]

%data = randsrc(1, nLength, alphabet);
%data = [3 3 0 3 3 2 0 2 2 0 3 1 0 0 1 0 1 2 1 3 3 0 3 1 0 3 2 0 3 3 2 2 3 3 0 3 0 3 0 0 3 0 1 2 2 2 0 0 0 0 2 3 3 3 1 3 1 2 0 1 0 0 2 2 2 1 0 1 0 3 3 1 3 2 2 2 3 3 1 1 1 1 3 0 2 2 1 3 3 3 3 2 0 0 1 0];

% pskmod just only map every bit to constellation point
% NO bit grouping into SYMBOL!
%   - if alphabet is only [0 1] it's using points 1+0i and 0+i
%signal = pskmod(data, M, pi/4);

fs = 40e3;
ts = 1/fs;
t = [0:ts:255*ts];
signal = 4.2*exp(i*2*pi*3e3*t) + 6.2*exp(i*2*pi*8e3*t);
signal = conj(signal);              %conjugte signal to get the right bottom spectral part

figure;
plot(real(signal), 'b-x'); hold on;
plot(imag(signal), 'r-x'); hold off;
title('Original signal');
legend('I','Q')

figure;
ofdmSignal256 = ifft(signal,256);           %normal IFFT
ofdmSignal512 = ifft(signal,512);           %added zeros = resampled output signal 2times
ofdmSignal768 = ifft(signal,768);           %3times resampled
plot(abs(ofdmSignal256), 'b-o'); hold on;
plot(abs(ofdmSignal512), 'r-x');
plot(abs(ofdmSignal768), 'g-x');
ofdmSignalResampled2 = resample(ofdmSignal256, 2, 1);    %resampled output signal
ofdmSignalResampled3 = resample(ofdmSignal256, 3, 1);

%resampled signal spectrum after multiplication (to preserve signal energy) MUST BE SAME AS ofdmSignal512
plot(abs(ofdmSignalResampled2)/2, 'black-o');
plot(abs(ofdmSignalResampled3)/3, 'yellow-o'); hold off;
title('OFDM signal'); xlabel('f[bin]'); ylabel('A[-]');

