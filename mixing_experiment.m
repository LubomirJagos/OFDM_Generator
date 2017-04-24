clear all;

sampleFreq = 400e3;
sampletime = 1/sampleFreq;         % LeCroy sampling time
fftPoints = 256;

% example of mixing two singals in Matlab
f1 = 12.42e3;
f2 = 15.78e3;
f3 = 18.10e3;
fc = 42.42e3;
t = [0:sampletime:10/f2];
yc = sin(2*pi*fc.*t);
y1 = 1*sin(2*pi*f1.*t);
y2 = 0.5*sin(2*pi*f2.*t);
y3 = 0.3*sin(2*pi*f3.*t);

y = y1 + y2 + y3;

spec = fft(y,fftPoints);
spec2 = fft(y .* yc,fftPoints);
xFreqAxis = linspace(0,1/sampletime,fftPoints);


figure
stem(xFreqAxis, abs(spec), 'x');
hold on;
stem(xFreqAxis, abs(spec2),'rx');

