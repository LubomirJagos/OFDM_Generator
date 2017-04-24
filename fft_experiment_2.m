clear all;

%%
% example of decimation, adding zeros to the end of original signal and
% say to fft that sampletime is few times smaller to move spectrum to
% the higher frequencies and whatch what would happen with its bandwidth

sampletime = 10e-6;
t = [0:sampletime:1e-3];
f1 = 10.0e3;
f2 = 15.0e3;
y = sin(2*pi*f1.*t) + 0.4*sin(2*pi*f2.*t);

% oversampling original signal by adding zero on the end

overSample = 3;
y2 = [y (overSample-1)*zeros(1, length(y))];
%yHalf = floor(length(y2)/2)
%y2 = [y2(1:yHalf) zeros(1, length(y2)) y2(yHalf+1:end)];

figure
plot(y2, 'r');
hold on;
plot(y)

fftPoints = 256;
spec = fft(y,fftPoints);
spec2 = fft(y2,fftPoints);
specX = linspace(0,1/sampletime,fftPoints);
specX2 = linspace(0,1/(sampletime/overSample),fftPoints);

figure
plot(specX, abs(spec));
hold on;
plot(specX2, abs(spec2),'r');

figure
newY = ifft(spec, fftPoints);
plot(newY, 'r');
hold on;
plot(y);

