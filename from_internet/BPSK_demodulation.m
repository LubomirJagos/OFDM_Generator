% [ch1Data, ch2Data]= Hantek6022BE_ReadingData(0,0,rxLength,hantekCh1VoltDiv, hantekCh2VoltDiv,hantekSampleFreqInt,0,0,0,0,0,0,0,0);
% ModulatedSens = double(ch1Data);
% ModulatedSens = Modulated/max(Modulated);
% 
% ModulatedSens = resample(double(ch1Data), fs, hantekSampleFreq);
% tc = [0:ts:(length(ModulatedSens)-1)*ts];
% 
% figure;
% plot(ModulatedSens, 'r-x');
% title('RX data before cut');
    
%                   ^^^^
%                   ||||
%  to sense data uncomment top part 
%  --------------------------------------                                  
%  to evaluate data uncomment bottom part
%                   ||||
%                   vvvv

%%%LuboJ.
Modulated = ModulatedSens(1:end);
figure;
subplot(211);
plot(Modulated, 'r-x');
title('RX data after cut');

tc=0:ts:ts*(length(Modulated)-1);
fc=1e3;
%%%%%%%%%

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

subplot(212);
stem(received)
title('Impulses of Received bits');
xlabel('Time (seconds)-->');
ylabel('Amplitude (volts)')

