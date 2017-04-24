close all;
clear all;

fs = 100e3;
ts = 1/fs;
f = 42e3;
nLen = 10;

t = [0:ts:(nLen-1)*ts];
y = exp(j*2*pi*f*t);

figure;
subplot(211);plot(real(y),'-x');
subplot(212);plot(imag(y),'r-x');

y = resample(y,2,1);

figure;
subplot(211);plot(real(y),'-x');
subplot(212);plot(imag(y),'r-x');
