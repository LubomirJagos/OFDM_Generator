clear all;
close all;

addpath('..');

%Nb is the number of bits to be transmitted
T=1e-3;                 %Bit rate is assumed to be 1 bit/s;
samplesPerBit = 1000;
%bits to be transmitted
b=[1 0 1 0 1 0 0 0 1 0 1 0 0 1 1 0 0]
ts = T/samplesPerBit;
fs = 1/ts;
%Rb is the bit rate in bits/second
 
 
NRZ_out=[];
RZ_out=[];
Manchester_out=[];
  
%Vp is the peak voltage +v of the NRZ waveform
Vp=1;
%Here we encode input bitstream as Bipolar NRZ-L waveform
for index=1:size(b,2)
 if b(index)==1
 NRZ_out=[NRZ_out ones(1,samplesPerBit)*Vp];
 elseif b(index)==0
 NRZ_out=[NRZ_out ones(1,samplesPerBit)*(-Vp)];
 end
end


% Generated bit stream impulses
figure(1);
stem(b);
xlabel('Time (seconds)-->')
ylabel('Amplitude (volts)-->')
title('Impulses of bits to be transmitted');
% NOT NEEDED TO SEE NRZ data
% figure(2);
% plot(NRZ_out);
% xlabel('Time (seconds)-->');
% ylabel('Amplitude (volts)-->');
% title('Generated NRZ signal');
 
tc=0:ts:T*length(b)-ts;
%Frequency of the carrier
fc=1e3;                     %fc << fs
%Here we generate the modulated signal by multiplying it with 
%carrier (basis function)

%Modulated=NRZ_out.*(sqrt(2/T)*cos(2*pi*fc*tc));
Modulated=NRZ_out.* cos(2*pi*fc*tc);

% not needed to see now
% figure; hold on;
% plot(tc, Modulated);
% xlabel('Time (seconds)-->');
% ylabel('Amplitude (volts)-->');
% title('BPSK Modulated signal');
% hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LuboJ. part

lecroySampleFreq = 400e3;
lecroySampleTime = 1/lecroySampleFreq;

hantekSampleFreqInt = 25;
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

Modulated = resample(Modulated, lecroySampleFreq, fs);          %signal on generator output
rxLength = floor(1.5*length(Modulated));                        %scope received signal, more samples to be sure that I have whole signal

figure
plot(Modulated, '-x'); title('Generator Modulated data resampled');
    
lw410 = LW410Interface(1, lecroySampleTime);           %gpibAddr, sampletime
lw410.wave_data(Modulated, 1);

disp('Waiting for generator');
pause(2);
disp('Done.');


[ch1Data, ch2Data]= Hantek6022BE_ReadingData(0,0,rxLength,hantekCh1VoltDiv, hantekCh2VoltDiv,hantekSampleFreqInt,0,0,0,0,0,0,0,0);
Modulated = double(ch1Data);
%Modulated = Modulated/max(Modulated);      %NEEDED?? normalization???

Modulated = resample(double(ch1Data), fs, hantekSampleFreq);
tc = [0:ts:(length(Modulated)-1)*ts];

figure;
plot(Modulated, 'r-x');
title('RX data');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


y=[];
%We begin demodulation by multiplying the received signal again with 
%the carrier (basis function)
%demodulated=Modulated.*(sqrt(2/T)*cos(2*pi*fc*tc));
demodulated=Modulated .* cos(2*pi*fc*tc);

%Here we perform the integration over time period T using trapz 
%Integrator is an important part of correlator receiver used here

%LuboJ.
%there is no synchronization, so I could be on any sample for beginning
%no iteration over whole samples there is same samples left in the end,
%to be sure not exceed matrix dimension

for i=1:samplesPerBit:length(demodulated)-samplesPerBit
 y=[y trapz(tc(i:i+samplesPerBit-1),demodulated(i:i+samplesPerBit-1))];
end
received=y>0;
figure;
stem(received)
title('Impulses of Received bits');
xlabel('Time (seconds)-->');
ylabel('Amplitude (volts)')