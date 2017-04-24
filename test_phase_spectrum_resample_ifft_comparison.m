f = 4e6;
sampletime = 1/f;

nLength = 64;
resampleFactor = 7;

t = 0:sampletime:nLength*sampletime;

A = 42.7;
f1 = 3e6;
y = A*sin(2*pi*f1*t);
yResampled = resample(y,resampleFactor,1);

yIfft = ifft(y, nLength*resampleFactor);
yIfftResampled = resample(ifft(y),resampleFactor,1);

figure;
subplot(311); plot(y, '-x');
title('Signal in time');
subplot(312); plot(abs(fft(y, 256)), '-x');
title('Signal Spectrum');
subplot(313); plot(angle(fft(y, 256)), '-x');
title('Signal phase');
figure;
subplot(311); plot(yResampled, '-x');
title('Signal in time');
subplot(312); plot(abs(fft(yResampled, 256)), '-x');
title('Signal Spectrum');
subplot(313); plot(angle(fft(yResampled, 256)), '-x');
title('Signal phase');

figure;
subplot(311); plot(downsample(yResampled,resampleFactor), '-x');
title('Signal in time');
subplot(312); plot(abs(fft(downsample(yResampled,resampleFactor), 256)), '-x');
title('Signal Spectrum');
subplot(313); plot(angle(fft(downsample(yResampled,resampleFactor), 256)), '-x');
title('Signal phase');


% For IFFT

fftPoints = 512;

figure;
subplot(311); plot(abs(yIfft), '-x');
title('Signal in time');
subplot(312); plot(abs(fft(yIfft, fftPoints)), '-x');
title('Signal Spectrum');
subplot(313); plot(angle(fft(yIfft, fftPoints)), '-x');
title('Signal phase');
figure;
subplot(311); plot(abs(yIfftResampled), '-x');
title('Signal in time');
subplot(312); plot(abs(fft(yIfftResampled, fftPoints)), '-x');
title('Signal Spectrum');
subplot(313); plot(angle(fft(yIfftResampled, fftPoints)), '-x');
title('Signal phase');





