clear all;
close all;

leCroySampleFreq = 4e6; %HARDCODED HERE
data_source = [1 5 3 2 8 9 2 4 5 6 3 9 0 10 3 11 15 2 0 2 4 3 12 13 15 3 2 7 10 2 5 7 8 1 10 11 2 15 13];
%data_source = [2 3 7 1 2 5 4 7 5 6 1 0 1 2 5 0 4 3 0 5 5 3 2 6 5 4 0 5 6 3 2 5 1 0 6 7 1 7 3 7 5 5 1 0 2 0 6 1 2 7 3 4 1];

hantekSampleFreqInt = 13;
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

rxLength = 130e3;

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

[ch1Data, ch2Data, triggerPointIndex]= Hantek6022BE_ReadingData(0,0,rxLength,hantekCh1VoltDiv, hantekCh2VoltDiv,hantekSampleFreqInt,triggerSweep,triggerSource,triggerLevel,triggerSlope,nHTrigPosition,0,0,0);
rxData = double(ch1Data);

dlmwrite('RXQPSK128carriers.dat',rxData,',');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DEMODULATION
%   QPSK
%       - 128 carriers
%       - 13 bits for cyclic prefix
%       - 4 times upsampled

resampleFactor = 32;

figure(1);
subplot(311);
plot(rxData, 'r-x');
hold on; plot(double(ch2Data), 'g-x'); hold off;
title('RX data');

disp(horzcat('Triggered by point x = ', num2str(triggerPointIndex)));

subplot(312);
triggerPointIndex = triggerPointIndex;
ofdm_signal = rxData(triggerPointIndex:triggerPointIndex+length(data_source)*resampleFactor-1);
plot(ofdm_signal, '-x');
title('Cutted OFDM signal');

subplot(313);
plot(abs(fft(ofdm_signal, 1024)), '-x');
title('RX ofdm signal');

fc = 0.7e6;
tc = [0:hantekSampleTime:(length(ofdm_signal)-1)*hantekSampleTime];
demixedSignal = ofdm_signal .* exp(-i*2*pi*fc*tc);

filtSignal = fft(demixedSignal);
filtSignal(50:1200) = 0;
filtSignal = ifft(filtSignal);

figure;
plot(abs(fft(filtSignal)), 'g-x');


xAxis = linspace(0,4e6,1024);
figure;
subplot(311);
plot(real(demixedSignal), 'b-x'); hold on;
plot(imag(demixedSignal), 'r-x'); hold off;
title('I and Q part of ofdm signal');
subplot(312); plot(xAxis, abs(fft(demixedSignal, 1024)), '-x');
title('RX ofdm signal spectrum');
subplot(313); plot(abs(fft(filtSignal, 1024)), '-x');
title('Filtered signal spectrum');

demixedSignal = resample(filtSignal,1,resampleFactor);
fftSignal = fft(demixedSignal);
demodSignal = qamdemod(fftSignal,16,pi/2);
%demodSignal = qamdemod(fftSignal,8,pi/2);

figure;
stem(data_source); hold on;
stem(demodSignal, 'r-x'); hold off;
legend('TX symbols', 'RX symbols');

