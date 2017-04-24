%   THIS EXAMPLE GENERATE SIGNAL AT FM Radio frequencies, if wire is
%   connected to the output it's received by usual radio station and you
%   can listen pitch (there is too less samples do something more,
%   at 400MHz and 250k samples we can modify just 250e3/400e6 appr. 0.5us)

clear all;
close all;

lecroySampleFreq = 400e6;
lecroySampleTime = 1/lecroySampleFreq;

fs = 400e6;
ts = 1/fs;
nLength = 250e3;

fm1 = 9e3;
fm2 = 4.2e3;
fm3 = 6.2e3;

fc = 103e6;
mi = 5;

t=0:ts:ts*(nLength-1);

% m = sin(2*pi*fm1*t) + ...
%     2*sin(2*pi*fm2*t) + ...
%     3*sin(2*pi*fm3*t);
m = sin(2*pi*fm1*t);

subplot(3,1,1);
plot(t,m);
xlabel('Time');
ylabel('Amplitude');
title('Message Signal');
grid on;

c=sin(2*pi*fc*t);
subplot(3,1,2);
plot(t,c);
xlabel('Time');
ylabel('Amplitude');
title('Carrier Signal');
grid on;

y=sin(2*pi*fc*t+(mi.*m));%Frequency changing w.r.t Message
subplot(3,1,3);
plot(t,y);
xlabel('Time');
ylabel('Amplitude');
title('FM Signal');
grid on;

%txSig = resample(y, lecroySampleFreq, fs);          %signal on generator output
txSig = y;
lw410 = LW410Interface(1, lecroySampleTime);           %gpibAddr, sampletime
lw410.wave_data(txSig, 1);
