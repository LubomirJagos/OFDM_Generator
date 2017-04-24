close all
clear all

oversampling = 32;
numCarriers = 18;
numSymbols = 4;

fs = 4e6;
ts = 1/fs;
data_source = [1 2 0 2 1 0 3 1 2 3 1 2 0 0 0 0 0 0 0 0 0 0 3 1 0 1 2 2 1 0 2 1 2 2 2 2 2 2 2 2 2 2 1 0 0 1 2 3 3 2 1 0 0 1 2 3 3 2 1 0 0 1 2 3 0 1 2 3 0 1 2 3];
%data_source = randsrc(1,72,[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]);
%data_source = randsrc(1,72,[0 1 2 3 4 5 6 7]);

figure;
plot(data_source);
title('Original');

hModulator =  comm.QPSKModulator('PhaseOffset', pi/4,'SymbolMapping','Binary');
%hModulator =  comm.RectangularQAMModulator(16);

hDemodulator =  comm.QPSKDemodulator('PhaseOffset', pi/4+pi,'SymbolMapping','Binary');
%hDemodulator =  comm.RectangularQAMDemodulator(16);
modData = step(hModulator, data_source.').';
data_matrix = reshape(modData, numCarriers, numSymbols);

ifft_data = ifft(data_matrix,numCarriers);
[rows_ifft_data cols_ifft_data]=size(ifft_data);
len_ofdm_data = rows_ifft_data*cols_ifft_data;
ofdmSignal = reshape(ifft_data, 1, len_ofdm_data);

figure;
plot(real(ofdmSignal), 'b-x'); hold on;
plot(imag(ofdmSignal), 'r-x'); hold off;


% IQ modulation
ofdmSignal = resample(ofdmSignal,oversampling,1);   %resample to have more samples per symbol
fc = 700e3;
tc = [0:ts:(length(ofdmSignal)-1)*ts];
%mixedSignal = ofdmSignal .* exp(j*2*pi*fc*tc);
mixedSignal = ofdmSignal .* (cos(2*pi*fc*tc)+j*sin(2*pi*fc*tc));



mixedSignal = 2*real(mixedSignal);
%mixedSignal = awgn(mixedSignal,33,0);   %added noise


xAxis = linspace(0,fs,1024);
figure;
plot(xAxis, abs(fft(ofdmSignal, 1024))); hold on;
plot(xAxis, abs(fft(mixedSignal, 1024)), 'r'); hold off;


%For QAM(16) no carrier ERROR!
%fc = 699.88e3;
deMixedSignal = mixedSignal .* exp(-j*2*pi*fc*tc);
%deMixedSignal = mixedSignal .* (cos(2*pi*fc*tc) + j*sin(2*pi*fc*tc));

deMixedSignal = fft(deMixedSignal);
deMixedSignal(end/4:3*end/4) = 0;
deMixedSignal = ifft(deMixedSignal);

%deMixedSignalResampled = resample(deMixedSignal,1,oversampling);
 deMixedSignalResampled = downsample(deMixedSignal(1:end),oversampling);

figure;
plot(abs(fft(deMixedSignal))); hold on;
title('deMixedSignal');

deMixedSignalResampled = reshape(deMixedSignalResampled, numCarriers, numSymbols);
fftSignal = fft(deMixedSignalResampled);
fftSignal = reshape(fftSignal, 1, numCarriers*numSymbols);
demodSignal = step(hDemodulator, fftSignal.').';

figure;
plot(demodSignal);
title('Demodulated');

figure;
plot(real(fftSignal), imag(fftSignal),'b.');
title('Received data');

figure;
subplot(221); plot(real(modData));
subplot(223); plot(imag(modData), 'r');
subplot(222); plot(real(fftSignal));
subplot(224); plot(imag(fftSignal), 'r');
title('Modulated signal, fftSignal comparison');

