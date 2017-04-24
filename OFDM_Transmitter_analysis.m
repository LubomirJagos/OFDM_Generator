clear all;
f = 4.242e3;                %carrier frequency
f2 = 1.786e3;
T = 1/f;                    %period
dt = T/24;
N = 1e3;                  %num. of samples
t = 0:dt:N*dt;               %generating time axis
A = 2.7e3;
A2 = 1.3e3;
fftPoints = 2^10;

wave = A*sin(2 * pi * f .* t);
%wave2 = A2*sin(2 * pi * f2 .* t);
%wave = wave + wave2;


figure(1);
subplot(311); plot(t, wave);

waveSpec = abs(fft(wave, fftPoints));
wavePhase = angle(fft(wave, fftPoints))
freqAxis = linspace(0, 1/dt, fftPoints);

subplot(312); stem(freqAxis, waveSpec);
subplot(313); stem(freqAxis, log10(waveSpec),'r');

figure(2);
stem(freqAxis, wavePhase);

%stem(real(fft));
