%test LW410Interface, sending data directly into device.

close all;
clear all;

%gpib_addr = 0
lw410 = LW410Interface();
lw410.sampletime = 250e-9;
nLength = 500;
f = 0.3e6;
t = [0:lw410.sampletime:lw410.sampletime*(nLength-1)];

%LW410 level set is RMS value not amplitude!   WRONG!
%When increase bandpass filter (max 10MHz) amplitude is 2times higher!
y = 2.33*sin(2*pi*f.*t)/2;

lw410.wave_data(y, 1);              %waveData, outputEnabled?

figure
plot(t,y);

