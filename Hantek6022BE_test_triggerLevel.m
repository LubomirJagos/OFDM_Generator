clear all;
close all;

%Nb is the number of bits to be transmitted
T=1e-3;                 %'value' rate is assumed to be 1 value/ms;
samplesPerValue = 500;
txData=[0 5 0 1 2 1 0 1 2 3 2 1 0]/2;

ts = T/samplesPerValue;
fs = 1/ts;

auxTxData = [];
for i = 1:length(txData)
    auxTxData = [auxTxData txData(i)*ones(1,samplesPerValue)];
end
txData = auxTxData;
clear auxData;

lecroySampleFreq = 400e3;
lecroySampleTime = 1/lecroySampleFreq;

hantekSampleFreqInt = 14;
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

txData = resample(txData, lecroySampleFreq, fs);          %signal on generator output
rxLength = floor(1.5*(hantekSampleFreq/lecroySampleFreq)*length(txData));                        %scope received signal, more samples to be sure that I have whole signal

figure(1);
subplot(211);
plot(txData, '-x'); title('Generator Modulated data resampled');
    
lw410 = LW410Interface(1, lecroySampleTime);           %gpibAddr, sampletime
lw410.wave_data(txData, 1);

disp('Waiting for generator');
pause(2);
disp('Done.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Reading data, here is where actually trigger level test is happening

triggerSource = 1;
triggerSweep = 2;
    % 0 = AUTO
    % 1 = Normal
    % 2 = SINGLE
triggerLevel = 30;          % <---- should be more than zero to run correctly!
triggerSlope = 0;
    % 0 = rise
    % 1 = fall
    
nHTrigPosition = 0;         % 0 - 100

a = 0;
[ch1Data, ch2Data, triggerPointIndex]= Hantek6022BE_ReadingData(0,0,rxLength,hantekCh1VoltDiv, hantekCh2VoltDiv,hantekSampleFreqInt,triggerSweep,triggerSource,triggerLevel,triggerSlope,nHTrigPosition,0,0,0);
rxData = double(ch1Data);

% rxData = resample(double(ch1Data), fs, hantekSampleFreq);
% tc = [0:ts:(length(rxData)-1)*ts];

figure(1);
subplot(212);
plot(rxData, 'r-x');
hold on; plot(double(ch2Data), 'g-x'); hold off;
title('RX data');

disp(horzcat('Triggered by point x = ', num2str(triggerPointIndex)));


