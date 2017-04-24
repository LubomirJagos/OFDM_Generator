close all;
clear all;

%% Start conditions

lecroySampleFreq = 400e3;
lecroySampleTime = 1/lecroySampleFreq;
dcOffset = 0;
ampl = 1.5;
f = 50e3;
f2 = 57e3;
f3 = 112e3;
f4 = 223e3;
f5 = 40e3;
txLength = 600;
rxLength = 1024;
rxFFTPoints = 1024;

hantekSampleFreqInt = 14;
hantekSampleFreq = Hantek6022BE_GetSampleFreq(hantekSampleFreqInt);
hantekSampleTime = 1/hantekSampleFreq;
    % 0 - 10  = 48MHz
    % 11      = 16MHz
    % 12      = 8MHz
    % 13      = 4MHz
    % 14 - 24 = 1MHz

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
    
%PSK settings
M = 4;

%% Generate and upload wave

lw410 = LW410Interface(1, lecroySampleTime);           %gpibAddr, sampletime

%sinewave, sampletime is from lecroy generator
t = [0:lecroySampleTime:(txLength-1)*lecroySampleTime];
% y = dcOffset +              ...
%     ampl*sin(2*pi*f.*t) +   ...
%     ampl*sin(2*pi*f2.*t) +  ...
%     ampl*sin(2*pi*f3.*t) +  ...
%     ampl*sin(2*pi*f4.*t) +  ...
%     ampl*sin(2*pi*f5.*t)    ...
% ;

y = dcOffset +              ...
    ampl*sin(2*pi*f.*t)     ...
;

% Ukazka resampling, pretoze mam svoj signal ktory bol samplovany nejakym
% kmitoctom a chcem ho vyslat cez LW410, tak musim spravit resample na
% jeho kmitocet.
%
% Apropo pri prijme signalu na Hantek ho osciloskop vidi tak ako je v reale
% lebo uz bol prisposobeny na vyslanie signal, takze ho staci len prijat a
% vsetky veliciny su ok. (len NEMUSI PLATIT ze fs LW410 == Hantek !)
%
% f = 47.42e3;
% mySampleFreq = 42.42e3;
% mySampleTime = 1/mySampleFreq;
% myT = [0:mySampleTime:(txLength-1)*mySampleTime];
% y = dcOffset +              ...
%     ampl*sin(2*pi*f.*myT)     ...
% ;
% y = resample(y, lecroySampleFreq, mySampleFreq);
% t = linspace(0,length(y)*lecroySampleTime, length(y));



figure
plot(t, y, '-x'); grid on;
title('Generated signal samples'); 

%should be PSK
% data = randsrc(1, nLength, [0:M-1]);
% dataMod = pskmod(data, M);
% figure
% hold on;
% title('PSK generated data');
% plot(real(dataMod), 'b-x');
% plot(imag(dataMod), 'r-x')
% hold off;

lw410.wave_data(y, 1);

disp('Waiting for generator');
pause(2);
disp('Done.');

%% Reading set data

[ch1Data, ch2Data]= Hantek6022BE_ReadingData(0,0,rxLength,hantekCh1VoltDiv, hantekCh2VoltDiv,hantekSampleFreqInt,0,0,0,0,0,0,0,0);
ch1Data = double(ch1Data);
ch2Data = double(ch2Data);
trx = 0:hantekSampleTime:(rxLength-1)*hantekSampleTime;

figure
hold on; grid on;
plot(trx, ch1Data, 'b-x');
plot(trx, ch2Data, 'r-x');
hold off;

specCh1Data = abs(fft(ch1Data, rxFFTPoints));
specCh2Data = abs(fft(ch2Data, rxFFTPoints));
fAxisRx = linspace(0,hantekSampleFreq, rxFFTPoints);

figure
hold on; grid on;
plot(fAxisRx, specCh1Data, 'b-x');
plot(fAxisRx, specCh2Data, 'r-x');
hold off;
