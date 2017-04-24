%LuboJ.
close all;
clear all;

hantekSampleFreqInt = 13;
rxLength = 4000;
hantekSampleFreq = Hantek6022BE_GetSampleFreq(hantekSampleFreqInt);
hantekSampleTime = 1/hantekSampleFreq;
    % 0 - 10  = 48MHz
    % 11      = 16MHz
    % 12      = 8MHz
    % 13      = 4MHz
    % 14 - 24 = 1MHz
    % 25 = 500kHz
    % 26 = 200kHz
    % 27 = 100kHz

hantekCh1VoltDiv = 6;
hantekCh2VoltDiv = 6;
    % 0 = 20mV
    % 1 = 50mV
    % 2 = 100mV
    % 3 = 200mV
    % 4 = 500mV
    % 5 = 1V
    % 6 = 2V
    % 7 = 5V

[ch1Data, ch2Data]= Hantek6022BE_ReadingData(0,0,rxLength,hantekCh1VoltDiv, hantekCh2VoltDiv,hantekSampleFreqInt,0,0,0,0,0,0,0,0);
Modulated = double(ch1Data);

figure;
plot(Modulated, '-x');
title('RX data');

figure;
specX = linspace(0,hantekSampleFreq,4096);
plot(specX, abs(fft(Modulated, 4096)));
title('RX signal spectrum');