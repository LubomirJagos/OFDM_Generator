close all
clear all

oversampling = 32;
numCarriers = 18;
numSymbols = 4;

fs = 4e6;
ts = 1/fs;
data_source = randsrc(1,72,[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]);

figure;
plot(data_source);
title('Original');

hModulator =  comm.RectangularQAMModulator(16);
hDemodulator =  comm.RectangularQAMDemodulator(16);

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



%fc = 699.9e3;  %carrier frequency error
mixedSignal = real(mixedSignal);
mixedSignal = awgn(mixedSignal,27,0);   %SNR = 27dB, low noisy, WORST CASE! border to see it by eyes
% mixedSignal = awgn(mixedSignal,33,0);   %SNR = 33dB, added noise
% mixedSignal = awgn(mixedSignal,39,0);   %SNR = 39dB, added noise, bad
                                            % GOOD case


xAxis = linspace(0,fs,1024);
figure;
plot(xAxis, abs(fft(ofdmSignal, 1024))); hold on;
plot(xAxis, abs(fft(mixedSignal, 1024)), 'r'); hold off;

deMixedSignal = mixedSignal .* exp(-j*2*pi*fc*tc);
%deMixedSignal = mixedSignal .* (cos(2*pi*fc*tc) + j*sin(2*pi*fc*tc));

%%%%%FILTERING without this it's not running! filtering is neededing,
%%%%%otherwis high frequencies are creating noise
%%%Filter by zeroing FFT spectrum

%%%Filter by frequency spectrum zeroing FFT(y) --> zeroing --> IFFT(y)
%%%Better result than butterworth filter
deMixedSignal = fft(deMixedSignal);
deMixedSignal(end/4:3*end/4) = 0;
deMixedSignal = ifft(deMixedSignal);

%%%Filter by butterworth filer
% filtFc = 500e3;
% [b,a] = butter(2,filtFc/(fs/2));
% deMixedSignal = filter(b,a,deMixedSignal);
%%%%%END FILTERING

%deMixedSignalResampled = resample(deMixedSignal,1,oversampling);
 deMixedSignalResampled = downsample(deMixedSignal,oversampling);

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

