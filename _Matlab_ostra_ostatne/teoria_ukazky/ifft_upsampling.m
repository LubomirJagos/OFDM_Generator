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
% y = sin(2*pi*f*t);
%y = sawtooth(2*pi*f*t);
y = sin(2*pi*0.2e6*t) + 2.4*sin(2*pi*0.6e6*t) + 3.6*sin(2*pi*1.7e6*t) + 0.5*sin(2*pi*2.2e6*t);

specX = linspace(0,fs,N);
yfft = fft(y) ./ N;

upFactor = 23;
specX2 = linspace(0,upFactor*fs,upFactor*N);
yfft2 = [yfft(1:end/2) repmat(zeros(1,N),1,upFactor-1) yfft(end/2+1:end)];
y2 = ifft(yfft2) .* N * upFactor;   %to conserve energy!
t2 = 0:dt/upFactor:(upFactor*N-1)*dt/upFactor;

%% Normal signal

figure;
subplot(211);
plot(t,y);
title('Signal in time');
subplot(212);
plot(specX,abs(yfft));
title('Signal spectrum');

%% Upsampled signal

figure;
subplot(211);
plot(t2,y2);
title('Upsampled Signal in time');
subplot(212);
plot(specX2,abs(yfft2));
title('Upsampled Signal spectrum');

%% Comparison original vs. upsampled signal

figure;
plot(t,y,'-o');
title('Upsampled vs. normal signal');
grid on;
hold on;
plot(t2,y2,'-xr');
hold off;

figure;
plot(t2,resample(y,upFactor,1),'-o');
title('Upsampled vs. normal signal');
grid on;
hold on;
plot(t2,y2,'-xr');
hold off;
