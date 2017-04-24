clear all;
close all;

leCroySampleFreq = 4e6; %HARDCODED HERE
resampleFactor = 32;
rxLength = 8600;

%THIS:
%data_source = [0 1 2 3 0 1 2 0 1 0 3 1 3 1 0 2 1 3 1 0 1 2 1 2 1 0 3 1 0 3 1 1 2 2 1 0 1 1 0 2 1 3 1 2 1 3 0 1 2 0 1 2 1 3 3 3 1 3 3 3 2 2 0 1];
%OR THIS:
data_source = dlmread('sendData.dat',',');
data_source = pskdemod(data_source, 4, 0);  %READ WITHOUT PHASE OFFSET

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

fc = 0.7e6;
    %IQ demodulator frequency

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

dlmwrite('receivedData.dat',rxData,',');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DEMODULATION
%   QPSK
%       - 128 carriers
%       - 13 bits for cyclic prefix
%       - 4 times upsampled

figure(1);
subplot(311);
plot(rxData, 'r-x');
hold on; plot(double(ch2Data), 'g-x'); hold off;
title('RX data');

disp(horzcat('Received data length = ', num2str(length(rxData))));
disp(horzcat('Triggered by point x = ', num2str(triggerPointIndex)));

subplot(312);
triggerPointIndex = triggerPointIndex;
ofdm_signal = rxData(triggerPointIndex:triggerPointIndex+length(data_source)*resampleFactor-1);
plot(ofdm_signal, '-x');
title('Cutted OFDM signal');

subplot(313);
plot(abs(fft(ofdm_signal, 1024)), '-x');
title('RX ofdm signal');

tc = [0:hantekSampleTime:(length(ofdm_signal)-1)*hantekSampleTime];
origSignal = ofdm_signal .* exp(-i*2*pi*fc*tc);

filtSignal = origSignal;
filtSignal = fft(origSignal);
filtSignal(300:2000) = 0;
filtSignal = ifft(filtSignal);

xAxis = linspace(0,4e6,1024);
figure;
subplot(311);
plot(real(origSignal), 'b-x'); hold on;
plot(imag(origSignal), 'r-x'); hold off;
title('I and Q part of ofdm signal');
subplot(312); plot(xAxis, abs(fft(origSignal, 1024)), '-x');
title('RX ofdm signal spectrum');
subplot(313); plot(xAxis, abs(fft(filtSignal, 1024)), '-x');
title('Filtered signal spectrum');

filtSignal = resample(filtSignal,1,resampleFactor);
fftSignal = fft(filtSignal);
demodSignal = pskdemod(fftSignal,4, pi/2);

figure;
stem(data_source); hold on;
stem(demodSignal, 'r-x'); hold off;
legend('original','decoded');

