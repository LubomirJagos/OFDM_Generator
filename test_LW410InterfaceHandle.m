%test LW410Interface, sending data directly into device.

close all;
clear all;

%gpib_addr = 0
lw410 = LW410Interface();
lw410.sampletime = 2.5e-9;
aaa = lw410;
aaa.sampletime = 25e-6;

lw410.sampletime
aaa.sampletime



lw410.gpib_addr = 2;
nLength = 400;

disp(lw410.idn());

f = 0.4e3;
t = [0:lw410.sampletime:lw410.sampletime*nLength];
y = sin(2*pi*f.*t);
lw410.wave_data(y, 1);              %waveData, outputEnabled?

figure
plot(t, y, '-x');


