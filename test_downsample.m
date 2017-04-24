%% Test downsample() and myDownsample() implementation behaviour.
%
% Compare Matlab and my function how it is to see. I wanted to just know if
% it's the same when I take every nth sample from original signal.

close all;
clear all;

fs = 47e3;
ts = 1/fs;
nLen = 420;
downsampleFactor = 7;

% Input signal generation
f = 2.3e3;
t = [0:ts:(nLen-1)*ts];
y = sin(2*pi*f*t);

figure;
subplot(311); stem(y);

% Downsample using Matlab function
%yDownsampled = downsample(y,downsampleFactor);
yDownsampled = resample(y,1,downsampleFactor);
subplot(312); stem(yDownsampled, 'rx');         % Matlab downsample

% My implementation.
yDownsampled2 = [];
for l = 1:downsampleFactor:length(y)
    yDownsampled2 = [yDownsampled2 y(l)];
end

subplot(313); stem(yDownsampled2, 'kx');        % My downsample

% Spectrum graph to see downsample effect, frequency is going higher
% ntimes
xAxis = linspace(0,fs,256);
figure;
plot(xAxis, abs(fft(y, 256))); hold on;
xlabel('f[Hz]');
plot(xAxis, abs(fft(yDownsampled2, 256)), 'r'); hold off;





