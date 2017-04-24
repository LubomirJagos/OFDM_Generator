clear all;
close all;

addpath('..');

%Nb is the number of bits to be transmitted
T=1e-3;                 %Bit rate is assumed to be 1 bit/s;
samplesPerBit = 500;
%bits to be transmitted
b=[1 0 1 0 1 0 1 1 0 0 1 0]
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
subplot(211);
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
Modulated_memory = Modulated;                   %just to remember for analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%
%   FAKE SYNCHRONIZATION

Modulated = [zeros(1,samplesPerBit) ones(1,samplesPerBit) zeros(1,samplesPerBit) Modulated];

%%%%%%%%%%%%%%%%%%%%%%%%%%

% not needed to see now
% figure; hold on;
% plot(tc, Modulated);
% xlabel('Time (seconds)-->');
% ylabel('Amplitude (volts)-->');
% title('BPSK Modulated signal');
% hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LuboJ. part

% lecroySampleFreq = 400e3;
% lecroySampleTime = 1/lecroySampleFreq;
% 
% hantekSampleFreqInt = 24;
% hantekSampleFreq = Hantek6022BE_GetSampleFreq(hantekSampleFreqInt);
% hantekSampleTime = 1/hantekSampleFreq;
%     % 0 - 10  = 48MHz
%     % 11      = 16MHz
%     % 12      = 8MHz
%     % 13      = 4MHz
%     % 14 - 24 = 1MHz
%     % 25 = 500kHz
%     % 26 = 200kHz
%     % 27 = 100kHz
% 
% hantekCh1VoltDiv = 6;
% hantekCh2VoltDiv = 6;
%     % 0 = 20mV
%     % 1 = 50mV
%     % 2 = 100mV
%     % 3 = 200mV
%     % 4 = 500mV
%     % 5 = 1V
%     % 6 = 2V
%     % 7 = 5V
% 
% Modulated = resample(Modulated, lecroySampleFreq, fs);          %signal on generator output
% rxLength = floor(1.5*(hantekSampleFreq/lecroySampleFreq)*length(Modulated));                        %scope received signal, more samples to be sure that I have whole signal
% 
% subplot(212);
% plot(Modulated, '-x'); title('Generator Modulated data resampled');
%     
% lw410 = LW410Interface(1, lecroySampleTime);           %gpibAddr, sampletime
% lw410.wave_data(Modulated, 1);
% 
% disp('Waiting for generator');
% pause(2);
% disp('Done.');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
