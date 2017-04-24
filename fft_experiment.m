clear all;

sampletime = 10e-6;
t = [0:sampletime:1e-3];
f1 = 10.0e3;
f2 = 15.0e3;
y = sin(2*pi*f1.*t) + 0.4*sin(2*pi*f2.*t);
plot(y);

fftPoints = 256;
spec = fft(y,fftPoints);
spec2 = fftshift(fft(ifftshift(y), fftPoints));
specX = linspace(0,1/sampletime,fftPoints);

plot(specX, abs(spec));
hold on;
plot(specX, abs(spec2), 'r');
hold off;

figure
newY = ifft(spec, fftPoints);
plot(y);
figure
plot(newY, 'r');

