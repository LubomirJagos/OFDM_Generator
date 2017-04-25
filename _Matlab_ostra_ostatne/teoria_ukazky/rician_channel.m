close all; clear all;

%interp = 256;
%fs = 100e6/interp;
fs = 20e6;
ts = 1/fs;

%%  Generovanie signalu
dt = ts;

    % 420k vzoriek za 1780ms
    % 42k vzoriek za 220ms
    % 4.2k vzoriek za 85ms
    N = 4200;
    
f1 = 4.0e6;
f2 = 4.1e6;
f3 = 4.2e6;
f4 = 4.3e6;
f5 = 4.4e6;

f6 = 4.05e6;
f7 = 4.15e6;
f8 = 4.25e6;
t = [0:dt:dt*(N-1)];
sigin = exp(2*1j*pi*f1*t) + ...
        exp(2*1j*pi*f2*t) + ...
        exp(2*1j*pi*f3*t) + ...
        exp(2*1j*pi*f4*t) + ...
        exp(2*1j*pi*f5*t) + ...
        exp(2*1j*pi*f6*t) + ...
        exp(2*1j*pi*f7*t) + ...
        exp(2*1j*pi*f8*t);

figure;
plot(real(sigin));
hold on;
plot(imag(sigin),'-r');
hold off;
xlabel('N [-]');
ylabel('A [V]');
title('Vstupny signal');

%%  Vytvorenie modelu kanalu
%     CHAN = ricianchan(TS, FD, K) constructs a frequency-flat ("single
%     TS is the sample time of the input signal, in seconds.
%     FD is the maximum Doppler shift, in Hertz.
%     K is the Rician K-factor in linear scale.
%     CHAN = ricianchan(TS, FD, K, TAU, PDB) constructs a frequency-selective
%     
    % Doplerov posun urcuje sirku pasma ktora bude zarusena. Napr. ak
    % je signalom sinusovka na f = 4MHz a maximalny doplerovsky posun
    % 0.3MHz, tak vystupny signal zabera pasmo 3.7-4.3MHz
    maxDoplerShift = 100;


% simulacia viacerych ciest, prve pole tau je delay jednotlivych ciest,
% druhe pole je ich zisk
%chanModel = rayleighchan(ts, maxDoplerShift, [12e-6,5e-6,79e-6],[3,10,-5]);
% len jedna cesta
chanModel = rayleighchan(ts, maxDoplerShift);
set(chanModel, 'PathDelays', 10e-6);
%set(chanModel, 'AvgPathGaindB', 0);   %toto sa zatial nijak neprejavilo

%%  Signal po prechode kanalom
tic;
sigout = filter(chanModel, sigin);
sigout = awgn(sigout,10,'measured');
benchmark = toc;
disp(horzcat('Channel filter time: ', num2str(benchmark)));
      
xspec = linspace(0,fs,N);

figure;
plot(real(sigout));
xlabel('N [-]');
ylabel('A [V]');
title('Vystupny signal po prechode rayleighovym kanalom');
% hold on;
% plot(imag(sigout), '-r')
% hold off;

figure;
plot(xspec, 20*log10(abs(fft(sigout))/N));
xlabel('f [Hz]');
ylabel('A [dBV]');
title('Spektrum vstupneho a signalu po skresleni kanalom');
grid on;
hold on;
plot(xspec, 20*log10(abs(fft(sigin)/N)),'--r');
hold off;

'Energia vstupneho signalu (J, 1Ohm):'
sum(abs(sigin).^2/N)
'Energia signalu po prechode kanalom (J, 1Ohm)'
sum(abs(sigout).^2/N)




