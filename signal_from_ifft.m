clear all;

nLength = 1024;
ifftPoints = 1024;
sampleFreq = 4e6;
sampleTime = 1/sampleFreq;
spec = zeros(1, nLength);
%spec(73) = 1024;
spec(nLength - 73) = 1024;
figure
stem(spec);

y = ifft(spec, ifftPoints);
xTimeAxis = sampleTime:sampleTime:length(y)*sampleTime;
figure
plot(xTimeAxis, y); title('Signal in time');