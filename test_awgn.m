%% Test functio awgn(), its behaving

fs = 47e3;
ts = 1/fs;
nLen = 300;

f = 3.72e3;
t = [0:ts:(nLen-1)*ts];
y = sin(2*pi*f*t);
tNoise = awgn(y,3,0);       %add noise with SNR = dBW, signal is 0dBW

figure;
plot(y); hold on;
plot(yNoise, 'r'); hold off;

