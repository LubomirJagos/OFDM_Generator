clear all;
close all;

leCroySampleFreq = 4e6; %HARDCODED HERE
rxLength = 50000;
fc = 700e3;
    %IQ demodulator frequency
numCarriers = 17;
numSymbols = 30;
resampleFactor = 23;

hDemod = comm.QPSKDemodulator('PhaseOffset', pi/4+pi/2,'SymbolMapping','Binary');  %demodulate data from file
hDemod2 = comm.QPSKDemodulator('PhaseOffset', pi/4,'SymbolMapping','Binary'); %demodulate rxData

figure;
hDemod2.constellation;

%THIS:
% data_source = [0 1 2 3 0 1 2 0 1 0 3 1 3 1 0 2 1 3 1 0 1 2 1 2 1 0 3 1 0 3 1 1 2 2 1 0 1 1 0 2 1 3 1 2 1 3 0 1 2 0 1 2 1 3 3 3 1 3 3 3 2 2 0 1];
%OR THIS:
sendData = dlmread('sendData.dat',',');
data_source = step(hDemod, sendData).'  %READ WITHOUT PHASE OFFSET

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

%dlmwrite('receivedData.dat',rxData,',');

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
triggerPointIndex = triggerPointIndex+3;    %offset empirical 3 
ofdmSignal = rxData(triggerPointIndex:triggerPointIndex+length(data_source)*resampleFactor-1);
%ofdmSignal = dlmread('sendOFDMSeq.dat',',');     %signal as from generator
ofdmSignal = dlmread('sendOFDMSignal.dat',',');  %complex number signal

% % % % % % % filtFc = 100e3;
% % % % % % % [b,a] = butter(6,filtFc/(hantekSampleFreq/2));
% % % % % % % ofdmSignal = filter(b,a,ofdmSignal);

subplot(312);
plot(abs(ofdmSignal));
title('Cutted OFDM signal');

subplot(313);
plot(abs(fft(ofdmSignal, 1024)), '-x');
title('RX ofdm signal');

tc = [0:hantekSampleTime:(length(ofdmSignal)-1)*hantekSampleTime];
demixedSignal = ofdmSignal .* exp(-j*2*pi*fc*tc);
%origSignal = ofdmSignal .* (cos(2*pi*fc*tc) - j*sin(2*pi*fc*tc));

filtSignal = fft(demixedSignal);
filtSignal(end/4:3*end/4) = 0;
filtSignal = ifft(filtSignal);

xAxis = linspace(0,4e6,1024);
figure;
subplot(311);
plot(real(demixedSignal), 'b-x'); hold on;
plot(imag(demixedSignal), 'r-x'); hold off;
title('I and Q part of ofdm signal');
subplot(312); plot(xAxis, abs(fft(demixedSignal, 1024)), '-x');
title('RX ofdm signal spectrum');
subplot(313); plot(xAxis, abs(fft(filtSignal, 1024)), '-x');
title('Filtered signal spectrum');


%filtSignal = resample(filtSignal,1,resampleFactor);
filtSignal = downsample(filtSignal,resampleFactor);
filtSignal = reshape(filtSignal, numCarriers, numSymbols);
fftSignal = fft(filtSignal);
fftSignal = reshape(fftSignal, 1, numCarriers*numSymbols);
demodSignal = step(hDemod2, fftSignal.').';


figure;
stem(data_source); hold on;
stem(demodSignal, 'r-x'); hold off;
legend('original','decoded');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ADDITIONAL GRAPH
%constPoints = 1/sqrt(2)*[1+j -1+j -1-j 1-j];
constPoints = hDemod2.constellation*max(abs(fftSignal));

%normFftSignal = fftSignal/max(fftSignal);
wrongPoints = [];
for i = 1:length(demodSignal)
    if (data_source(i) ~= demodSignal(i))
%         wrongPoints = [wrongPoints normFftSignal(i)];
         wrongPoints = [wrongPoints fftSignal(i)];
    end
end


figure;
plot(real(fftSignal),imag(fftSignal),'b.', ...
    real(constPoints),imag(constPoints),'gx', ...
    real(wrongPoints),imag(wrongPoints),'r.');
grid on;


dlmwrite('rxDataFftSignal.dat',',');
figure;
subplot(221); plot(real(sendData), 'b');
title('real(sendData)');
subplot(223); plot(imag(sendData), 'r');
title('imag(sendData)');
subplot(222); plot(real(fftSignal), 'b');
title('real(rxModData)');
subplot(224); plot(imag(fftSignal), 'r');
title('imag(rxModData)');









