%clear all;
% close all;

leCroySampleFreq = 4e6; %HARDCODED HERE

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

figure;
subplot(311);
plot(rxData, 'r-x');
hold on; plot(double(ch2Data), 'g-x'); hold off;
title('RX data');

disp(horzcat('Triggered by point x = ', num2str(triggerPointIndex)));

subplot(312);
%ofdm_signal = resample(rxData(triggerPointIndex:triggerPointIndex+72192),1,4);
%ofdm_signal = rxData(triggerPointIndex:triggerPointIndex+72192);

triggerPointIndex = triggerPointIndex;
%ofdm_signal = rxData(triggerPointIndex:triggerPointIndex+resampleFactor*(32));

% Idea about taking more chunks and made average, but probably not good one
% for now.
%
% ofdm_signal = zeros(1,resampleFactor*32);
% for i = 1:3
%     ofdm_signal = ofdm_signal + rxData(triggerPointIndex+i*resampleFactor*32:triggerPointIndex+(i+1)*resampleFactor*(32)-1);
% end
% ofdm_signal = ofdm_signal / 3;


plot(ofdm_signal, '-x');
title('Cutted OFDM signal');

subplot(313);
plot(abs(fft(ofdm_signal, 1024)), '-x');
title('RX ofdm signal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ERROR is somewhere here in IQ demodulation
%

fc = 0.7e6;
tc = [0:hantekSampleTime:(length(ofdm_signal)-1)*hantekSampleTime];
%origSignal = ofdm_signal.*exp(i*2*pi*fc*tc);
origSignal = ofdm_signal .* (cos(2*pi*fc*tc) + i*sin(2*pi*fc*tc));

filtFc = 500e3;
[b,a] = butter(6,filtFc/(hantekSampleFreq/2));
filtSignal = filter(b,a,origSignal);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DEBUGGING%%%%%%%%%%%%%%%%%%%%%%%%%%
% origSignal = ofdm_signal;
% filtSignal = ofdm_signal;
% resampleFactor = 32;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
subplot(311);
plot(real(origSignal), 'b-x'); hold on;
plot(imag(origSignal), 'r-x'); hold off;
title('I and Q part of ofdm signal');
subplot(312); plot(abs(fft(origSignal, 1024)), '-x');
title('IQ signal before filtering spectrum');
subplot(313); plot(abs(fft(filtSignal, 1024)), '-x');
title('Filtered signal spectrum');

filtSignal = resample(filtSignal,1,resampleFactor);
fftSignal = fft(filtSignal);
demodSignal = pskdemod(fftSignal,4);

refSignal = [3 3 0 2 1 3 3 2 0 1 2 3 0 1 2 0 3 2 1 0 2 3 0 1 0 1 2 0 1 2 3 1 2 2 1 0 1 2 3 1 1 0 2 2 3 1 1 0 2 0 3 3 2 0 1 3 1 2 3 0 1 2 1 3];
figure;
stem(refSignal,'-x'); hold on;
stem(demodSignal,'r--x'); hold off;





