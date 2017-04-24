close all;
%% BPSK setting
%
% * bits, b=[1 0 1 0 1 0 0 1 1 0 0 1 1 0 0]
% * output coding NRZ

figure
stem(NRZ_out);
title('NRZ output signal');


%% Device setup:
%
% LeCroy LW410
%
% * lecroySampleFreq = 400e3;
% * leCroySampleTime = 2.5 us
%
% Hantek 6022BE
%
% * hantekSampleFreq = 500kHz
% * hantekSampleTime = 2 us
% * num. of samples (rxLength) = 9000
%
% Missing:
%
% * *forget to save input signal into generator*

figure;
set(gcf,'units','points','position',[10,10,800,520]);
subplot(211);
stem(NRZ_out);
title('NRZ intput signal');

subplot(212);
plot(demodulated, '-x');
title('Demodulated signal');
xlabel('#sample[-]');
ylabel('scope output [ ??? in progress ??? ]');

%% Received and transmitted bits
%
% BPSK demodulation:
%
% * <http://drmoazzam.com/matlab-code-bpsk-modulation-and-demodulation-with-explanation/ BPSK modulation and demodulation>
%
% Problems:
%
% * made correction are maybe not right (at mixing with carrier removed coefficient sqrt(1/T) )
% * synchronization not right, input and output signal shifted

figure;
subplot(121); stem(b, '-O');
title('Transmitted bits');
xlabel('bit[-]');
ylabel('bit value');
subplot(122); stem(received, 'r-O');
title('Received bits');
xlabel('bit[-]');
ylabel('bit value');

