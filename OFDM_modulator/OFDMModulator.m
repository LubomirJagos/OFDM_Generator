% Lubomir Jagos
%
%   Porovnanie vlastnej implementacie a Matlab OFDM bloku
%

close all;

hMod = comm.OFDMModulator( ...
    'FFTLength', 256, ...
    'NumGuardBandCarriers', [0;0], ...
    'CyclicPrefixLength', 0);

hMod2 = comm.OFDMModulator( ...
    'FFTLength', 256, ...
    'NumGuardBandCarriers', [12;15], ...
    'CyclicPrefixLength', 32);

data = [ones(1,64) zeros(1,64) ones(1,64) zeros(1,64)]';
data = randi([0 1],1,256)';
data = repmat([1 2 3 4 5 6 7 8],1,32)';
data2 = horzcat(repmat([1 2 3 4 5 6 7 8],1,28),[1 2 3 4 5])';
out = step(hMod,data);
out2 = step(hMod2,data2);
figure
plot(abs(fft(out)))
hold on;
plot(data,'--r')
hold off;
figure
hold on
plot(real(out))
plot(real(ifft(data)),'--r')
hold off

%spravit!

data2ifft = ifft(data2);

figure
plot(real(out2));
hold on;
plot(real(horzcat(  zeros(1,12),data2ifft(end-32+1:end)',data2ifft(33:end)',zeros(1,15)   )),'--r');
hold off;