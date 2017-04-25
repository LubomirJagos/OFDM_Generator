%  IFFT decimation
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
% y = sawtooth(2*pi*f*t);
%y = square(2*pi*f*t);
y = sin(2*pi*0.2e6*t) + 2.4*sin(2*pi*0.6e6*t) + 3.6*sin(2*pi*1.7e6*t) + 0.5*sin(2*pi*2.2e6*t);

specX = linspace(0,fs,N);
yfft = fft(y)/N;

y2 = [];
j = 1;
downFactor = 2.0;
for i = 1:downFactor:length(y)-1
    y2(j) = y(i);
    j = j + 1;
end
specX2 = linspace(0,fs/downFactor,N/downFactor);
yfft2 = fft(y2)/N;
t2 = 0:dt*downFactor:(N/downFactor-1)*dt*downFactor;

%% Normal signal

figure;
subplot(211);
plot(t,y);
grid on; title('Signal in time');
subplot(212);
plot(specX,abs(yfft));
grid on; title('Signal spectrum');

%% Upsampled signal

figure;
subplot(211);
plot(t2,y2);
grid on; title('Decimated Signal in time');
subplot(212);
plot(specX2,abs(yfft2));
grid on; title('Decimated Signal spectrum');

%% Comparison original vs. upsampled signal

figure;
plot(t,y,'-o');
grid on; title('Upsampled vs. normal signal');
grid on;
hold on;
plot(t2,y2,'-xr');
hold off;

