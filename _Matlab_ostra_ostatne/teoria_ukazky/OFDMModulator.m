% Lubomir Jagos
%
%   Porovnanie vlastnej implementacie a Matlab OFDM bloku
%

clear all; close all;

%%  Definicia modulatorov.
%   Jeden je len tak holy aby som videl ci robi ifft() a druhy ma nastavene
%   aj ostatne parametre.

hMod = comm.OFDMModulator( ...
    'FFTLength', 256, ...
    'NumGuardBandCarriers', [0;0], ...
    'CyclicPrefixLength', 0);

fftLen = 256;
leftGuard = 12;
rightGuard = 15;
prefixLen = 32;

hMod2 = comm.OFDMModulator( ...
    'FFTLength', fftLen, ...
    'NumGuardBandCarriers', [leftGuard;rightGuard], ...
    'CyclicPrefixLength', prefixLen);

%%   Porovnanie vystupu matlab modulatoru a normalneho prevedenia.
%       Bez ochranneho intervalu a cyklickeho prefixu je modulator len
%       ifft().

data = [ones(1,64) zeros(1,64) ones(1,64) zeros(1,64)]';
data = randi([0 1],1,256)';
data = repmat([1 2 3 4 5 6 7 8],1,32)';

%data2 = horzcat(repmat([1 2 3 4 5 6 7 8],1,28),[1 2 3 4 5])';
data2 =  horzcat(repmat([2 2 0 4 10 2 0 2],1,28),[1 2 3 4 5])';

out = step(hMod,data);
out2 = step(hMod2,data2);

% figure                      %porovnanie vystupu modulatoru a ifft()
% plot(abs(fft(out)))
% hold on;
% plot(data,'--r')
% hold off;
% figure                      %len tak na porovnanie vykreslenie realnej casti signalu
% hold on
% plot(real(out))
% plot(real(ifft(data)),'--r')
% hold off

%% Porovnanie OFDM s ochrannym intervalom a cyklickym prefixom
%   Takze aj systemovy modulator OFDM robi to iste co moj, otazka je ako
%   teraz spravne zostavit spektrum, pretoze existuju funkcie fftshift() a
%   ifftshift() :D

%data2ifft = ifft([ zeros(leftGuard,1); data2; zeros(rightGuard,1)]);
%data2out = [data2ifft(end-prefixLen+1:end); data2ifft];
data2ifft = ifft(ifftshift([ zeros(1,leftGuard) data2' zeros(1,rightGuard)]));
data2out = [data2ifft(end-prefixLen+1:end)'; data2ifft'];

figure;
plot(data2);

figure
plot(real(out2));
hold on;
plot(real(data2out),'--r');
hold off;

figure
plot(ifftshift(imag(out2)));
hold on;
plot(-ifftshift(imag(data2out)),'--r');
hold off;

figure;
plot(real(fftshift(fft(out2)'))); hold on;
plot(fliplr(real(fftshift(fft(data2out)'))), '--g'); hold off;
