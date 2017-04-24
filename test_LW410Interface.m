%test LW410Interface, sending data directly into device.

close all;
clear all;

%gpib_addr = 0
lw410 = LW410Interface();
lw410.sampletime = 25e-6;
nLength = 400;

disp(lw410.idn());
lw410.sampletime

f = 0.4e3;
t = [0:lw410.sampletime:lw410.sampletime*nLength];

%input signal has different waves in himself to see how LW410 is setting output
%from this test is obvious, that LW410 uses only real part
y = sawtooth(2*pi*f.*t) + i*square(2*pi*f.*t);

lw410.wave_data(y, 1);              %waveData, outputEnabled?

figure
plot(t, real(y), 'b-x', t, imag(y), 'r-x', t, abs(y), 'g-x');
legend('real', 'imag', 'abs');

