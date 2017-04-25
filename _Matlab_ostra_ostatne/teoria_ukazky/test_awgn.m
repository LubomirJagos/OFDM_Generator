% Bc. Lubomir Jagos, 5.4.2017
%
close all;

% usrp fs = 100MHz
fftLen = 4096;
interpolation = 1000;

%baseband signal oversampling, baseband = fs_usrp/interpolation
%signal BW = baseband / oversample
oversample = 2;

fs = 100e6/interpolation;
ds = 1/fs;

t = [0:fftLen-1] .* ds;
fc1 = 20e3;
fc2 = 30e3;
fc3 = 40e3;
fc4 = 50e3;
fc5 = 4e3;
sine1 = exp(2*pi*fc1*1j .* t)
sine2 = exp(2*pi*fc2*1j .* t)
sine3 = exp(2*pi*fc3*1j .* t)
sine4 = exp(2*pi*fc4*1j .* t)
sine5 = exp(2*pi*fc5*1j .* t)
%sine5 = 0.5*hilbert(sawtooth(2*pi*fc5*t));          %sawtooth aby som videl co to zrobi
                                                    %POZOR v spektre su
                                                    %harmonicke cize
                                                    %4,8,12,16,20,24,...
                                                    %kHz
sine = sine1 + sine2 + sine3 + sine4 + sine5;

%priprava na vypocet celkoveho vykonu signalu, najprv musim umocnit
%jednotlive nasnimane hodnoty (je jedno ci v casovej alebo frekvencnej
%oblasti) a potom kazdu vzorku vydelit dlzkou FFT, aby som dostal realne
%hodnoty
sigPwr = (fft(sine).^2 ) ./ fftLen;

% PRE SPRAVNE AWGN() JE NUTNE ZADAT SPRAVNU HODNOTU ENERGIE SIGNALU!

sineNoise = awgn(sine,50,10*log10(sum(abs(sine.^2))));
sigPwrNoise = (fft(sineNoise).^2)./fftLen;

specX = linspace(0,fs,fftLen);

figure;
plot(real(sineNoise),'r'); hold on;
plot(real(sine)); hold off;
figure;
plot(specX, 20*log10(abs(fft(sineNoise))/fftLen),'r'); hold on;
plot(specX, 20*log10(abs(fft(sine))/fftLen),'--'); hold off;


%vykon signalu je merany v signal processingu nad 1Ohm takze staci scitat a
%zlogaritmovat mocniny zosnimanych vzorkov. VYSLEDKY SA MUSIA ROVNAT

disp('Signal power from fft spec dBW');
10*log10(sum(abs(sigPwr)))
disp('Signal power of noisy sinewave, freq domain dBW');
10*log10(sum(abs(sigPwrNoise)))
disp('Signal power from time domain dBW');
10*log10(sum(abs(sine.^2)))
disp('Signal power of noisy sinewave, time domain dBW');
10*log10(sum(abs(sineNoise.^2)))
disp('Noise power dBW');
noisePwr = 10*log10(sum(abs(sineNoise.^2)) - sum(abs(sine.^2)))


%toto je len moja uvaha, a asi to ani nie je dobre, ale doriesit to!
disp('Average energy of noise component')
10^(noisePwr/10)/fftLen    %divided two times, first to get real power in real world second time to get for one spectral component
                                  %ak ma nejaka nosna ampl. 1 a su styri,
                                  %tak dokopy je to vykon 4W ak je SNR =
                                  %10dB, tak na kazdej nosnej bude mat sum
                                  %desatinu celkoveho vykonu, tj. zhruba
                                  %4*0.1 = 0.4V na kazdej nosnej, to je len
                                  %tak priblizne vysvetlene





