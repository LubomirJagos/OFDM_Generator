close all;
clear all;

M = 2;
nLength = 17;
sampleFreq = 100e3;
sampleTime = 1/sampleFreq;
fftPoints = 4096;

oversample = 21;            %oversampling for modulator
fc = 350e3;                 %modulator oscillator frequency

%% Generating random data signal
data = randsrc(1, nLength, [0:M-1]);
modData = pskmod(data, M);

figure
stem(data);

%signal analysis
spec = fft(modData, fftPoints);
xFreqAxis = linspace(0,sampleFreq, fftPoints);
figure
plot(xFreqAxis, log10(abs(spec)/max(abs(spec)))); title('pskmod spectrum'); xlabel('f[Hz]'); ylabel('A[-]');
hold on;

%% Upsampling data signal
sampleFreq = oversample*sampleFreq;
sampleTime = 1/sampleFreq;
nLength = oversample*nLength;

% generating carrier oscillator wave
tc = [sampleTime:sampleTime:nLength*sampleTime];
yI = cos(2*pi*fc.*tc);
yQ = sin(2*pi*fc.*tc);

xFreqAxis = linspace(0,sampleFreq, fftPoints);

dfs = sampleFreq/fftPoints;
carrierBin = ceil(fc/dfs);              %carrier frequency in spectrum bin order

%signal now gain energy because there is more samples => more energy
modData = resample(modData,oversample,1);

%% Simulatin IQ modulator
mixedSignalI = modData .* yI;
mixedSpecI = fft(mixedSignalI,fftPoints);

% cutting top part of mixed spectrum
mixedSpecI(1:carrierBin) = 0.1;
mixedSpecI(end/2:end) = 0.1;
mixedSignalI = ifft(mixedSpecI, fftPoints);

plot(xFreqAxis,log10(abs(mixedSpecI)/max(abs(mixedSpecI))), 'r'); title('Normalized mixed signal spectrum I');
hold off;

figure
plot(xFreqAxis,abs(fft(modData,fftPoints))); title('UPSAMPLED Modulation data spectrum');

% figure;
% plot(abs(fft(modData))./10);

% cutting top part for Q
mixedSignalQ = modData .* yQ;
mixedSpecQ = fft(mixedSignalQ,fftPoints);

figure
hold on;
%before cutting spectrum (high frequencies)
%plot(xFreqAxis,log10(abs(mixedSpecQ)/max(abs(mixedSpecQ))), 'b'); title('Normalized mixed signal spectrum Q');
mixedSpecQ(1:carrierBin) = 0.1;
mixedSpecQ(end/2:end) = 0.1;
mixedSignalQ = ifft(mixedSpecQ, fftPoints);
plot(xFreqAxis,log10(abs(mixedSpecQ)/max(abs(mixedSpecQ))), 'r'); title('Normalized mixed signal spectrum Q');
hold off;

%RESULTING MIXED SIGNAL IQ
mixedSignal = abs(mixedSignalI) + abs(mixedSignalQ);
figure
plot([sampleTime:sampleTime:sampleTime*length(mixedSignal)], mixedSignal); title('Result waveform IQ signal');
figure
hold on;
plot([sampleTime:sampleTime:sampleTime*length(mixedSignal)], abs(mixedSignalI)); title('Tx signal I,Q');
plot([sampleTime:sampleTime:sampleTime*length(mixedSignal)], abs(mixedSignalQ), 'r');
hold off;



% -------------------------------------------------------------------
%                   DEMODULATION
% -------------------------------------------------------------------
tc = [sampleTime:sampleTime:fftPoints*sampleTime];
yI = cos(2*pi*fc*tc);
yQ = sin(2*pi*fc*tc);

rxI = mixedSignal .* yI;
rxSpecI = fft(rxI,fftPoints);
rxQ = mixedSignal .* yQ;
rxSpecQ = fft(rxQ,fftPoints);

figure;
hold on;
plot(rxI, 'b-x')
plot(rxQ, 'r-x')
hold off;

carrierBin = 100e3;
dfs = sampleFreq/fftPoints;
carrierBin = ceil(fc/dfs);              %carrier frequency in spectrum bin order
rxSpecI(1:carrierBin) = 0.1;
rxSpecI(end/2:end) = 0.1;
rxSpecQ(1:carrierBin) = 0.1;
rxSpecQ(end/2:end) = 0.1;
rxI = ifft(rxSpecI, fftPoints);
rxQ = ifft(rxSpecQ, fftPoints);
rxSig = abs(rxI) + abs(rxQ);

figure
plot(rxSig); title('Received signal IQ abs'); hold on;

% red sinewave in graph
a = [sampleTime:sampleTime:fftPoints*sampleTime];
b = 100e3;
c = sin(2*pi*b.*a);
plot(c,'r'); hold off;

% rxI = resample(rxI,1,oversample);
% rxQ = resample(rxQ,1,oversample);

rxIresamp = [];
rxQresamp = [];
j = 1;
for i = 1:oversample:length(rxI)
    rxIresamp(j) = rxI(i);
    rxQresamp(j) = rxQ(i);
    j = j + 1;
end

scatterplot(horzcat(abs(rxIresamp)', abs(rxQresamp)'));





