close all;
clear all;

spec = zeros(1,1024);
%with only this items, result signal is complex signal (imaginary part is pi/2 shifted)
spec(824) = 1;
spec(700) = 1;

%with this spectral items, imaginary signal is zero
% spec(326) = 1;
% spec(202) = 1;

figure;
plot(spec);

y = ifft(spec);
figure;
imag(y)
plot(real(y), 'b-x'); hold on;
plot(imag(y), 'r-x'); hold off;

t = [0:10e-6:1000e-6]
%sign before cosine decides if we get up spectral part or low spectral part
y2a = 4.2*sin(2*pi*4.2e3*t) + i*4.2*cos(2*pi*4.2e3*t);
y2b = 4.2*sin(2*pi*4.2e3*t) - i*4.2*cos(2*pi*4.2e3*t);
spec2a = abs(fft(y2a, 1024));
spec2b = abs(fft(y2b, 1024));
spec2X = linspace(0,1/10e-6, 1024);
plot(spec2X, spec2a); hold on;
plot(spec2X, spec2b, '-r'); hold off;
legend('sign +', 'sign -');



