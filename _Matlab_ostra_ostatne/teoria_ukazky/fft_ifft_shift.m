%  IFFT upsampling
%
%  by LuboJ.
%
close all; clear all;

N = 100;
fs = 10e6;
dt = 1/fs;

f = 0.2e6;
t = 0:dt:(N-1)*dt;
y = sin(2*pi*1e6*t) + 2.4*sin(2*pi*2e6*t) + 3.6*sin(2*pi*4.5e6*t);

specX = linspace(0,fs,N);
yfft = fft(y) ./ N;
yfft2 = fft(fftshift(y)) ./ N;
yfft3 = fftshift(fft(y)) ./ N;
yfft4 = fft(ifftshift(y)) ./ N;
yfft5 = ifftshift(fft(y)) ./ N;

%% Normal signal

figure;
subplot(321);
stem(specX,abs(yfft));
title('just fft');

subplot(323);
stem(specX,abs(yfft2));
title('fft(fftshift(y))');

subplot(325);
stem(specX,abs(yfft3));
title('fftshift(fft(y))');

subplot(324);
stem(specX,abs(yfft4));
title('fft(ifftshift(y))');

subplot(326);
stem(specX,abs(yfft5));
title('ifftshift(fft(y))');
