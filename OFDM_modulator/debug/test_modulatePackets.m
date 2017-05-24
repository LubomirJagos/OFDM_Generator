%
%   LuboJ.
%
clear all; close all;

%
%   Create modulators.
%
headMod = comm.BPSKModulator('PhaseOffset', pi);
payloadMod = comm.QPSKModulator('PhaseOffset', 3/4*pi, 'BitInput', true);

fs = 200e3;
fftLen = 64;            %must be even!
packetLen = 96;
nProcessPackets = 1;

%
%   Generate data.
%
sync1 = [0., 0., 0., 0., 0., 0., 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 0., 0., 0., 0., 0.];
sync2 = [0, 0, 0, 0, 0, 0, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, 0, 1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 0, 0, 0, 0, 0];
occupiedCarriers = [(-26:-22) (-20:-8) (-6:-1) 1:6 8:20 22:26]+1;
pilotCarriers = [-21 -7 7 21]+1;
pilotSymbols = [1 1 1 -1];          %amplitudes

nSymbols = ceil(packetLen/length(occupiedCarriers));

for k = 1:length(occupiedCarriers)
    if (occupiedCarriers(k) < 1)
        occupiedCarriers(k) = occupiedCarriers(k) + fftLen;
    end
end
for k = 1:length(pilotCarriers)
    if (pilotCarriers(k) < 1)
        pilotCarriers(k) = pilotCarriers(k) + fftLen;
    end
end

% dataIn = randi(255,1,packetLen+4);
% dataIn = zeros(1,packetLen+4);

dataIn = [];
f = fopen('_test_mod_randseq_1.txt','r');

sigOut = [];
dPackets = [];
dFrames = [];
for j = 1:nProcessPackets
    dataIn = fread(f,packetLen,'uint8')';

    header = generateHeader(packetLen+4,j-1);       %plus 4 because there is CRC added in payload
    payload = de2bi(generatePayload(dataIn), 8)';
    payload = reshape(payload,1,numel(payload));

    modHeader = step(headMod, header')';
    modPayload = step(payloadMod, payload')';
    packet = [modHeader modPayload];
    dPackets = [dPackets packet];

    %
    %   Until here it looks good.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Allocator is NOT RUNNING CORRECTLY!
    %
    frameAux = zeros(1,fftLen);
    frame = [];
    for k = 1:nSymbols
        frameAux(occupiedCarriers) = modPayload(1:length(occupiedCarriers));
        frameAux(pilotCarriers) = pilotSymbols(1:length(pilotCarriers));
        frame = [frame frameAux];
    end
    sig = [sync1 sync2 modHeader frame];
    dFrames = [dFrames sig];

%     timeAxis = [1:1:length(sig)].*1/fs;
%     figure;
%     plot(timeAxis, real(sig));
%     title('OFDM real signal');
%     grid on;
%     hold on;
%     plot(timeAxis, imag(sig),'-rx');
%     hold off;

    %
    %   IFFT
    %
    for k = 1:fftLen:length(sig)-fftLen
    %     sigOut = [sigOut ifftshift(ifft(ifftshift(sig(k:k+fftLen-1))))];
        sigOut = [sigOut ifft(sig(k:k+fftLen-1))];
    end
end

fclose(f);
release(headMod);
release(payloadMod);

%
%   Upsampling
%   prechod od digitalu --> analog
%
%sigOut = resample(sigOut, 6,1);

timeAxis = [1:1:length(sigOut)].*1/fs;
figure;
plot(timeAxis, real(sigOut));
hold on;
plot(timeAxis, imag(sigOut),'-r');
hold off;

                    % df = fs/length(sigOut);
                    % specAxis = [-fs/2+df:df:fs/2];
                    % figure;
                    % plot(specAxis, 20*log10(abs(fft(sigOut))));
                    % grid on;

%
%   Debug checking diff btw. muxed stream in gnuradio and calculated in
%   matlab
%
f2 = fopen('_test_mod_muxed_1.txt','r');
gnuradioData = fread(f2,2*length(dPackets),'float32');
fclose(f2);
gnuradioData = reshape(gnuradioData,2,length(dPackets));
gnuradioData = gnuradioData(1,:) + 1i*gnuradioData(2,:);

figure;
plot(real(gnuradioData));
title('Output muxed data header and payload, GNURadio');
grid on;
hold on;
plot(imag(gnuradioData),'-r');
hold off;

figure;
plot(real(dPackets));
title('CALCULATED header and payload');
grid on;
hold on;
plot(imag(dPackets),'-r');
hold off;

figure;
plot(abs(dPackets-gnuradioData));
title('Comparison gnuradio, calculated from muxer');
ylim([-1 1]);

%
%   Debug check allocator.
%
f2 = fopen('_test_mod_allocated_1.txt','r');
gnuradioData = fread(f2,2*length(dFrames),'float32');
fclose(f2);
gnuradioData = reshape(gnuradioData,2,length(dFrames));
gnuradioData = gnuradioData(1,:) + 1i*gnuradioData(2,:);

figure;
plot(real(gnuradioData));
title('Output allocated carriers, GNURadio');
grid on;
hold on;
plot(imag(gnuradioData),'-r');
hold off;

figure;
plot(real(dFrames));
title('CALCULATED frame');
grid on;
hold on;
plot(imag(dFrames),'-r');
hold off;

figure;
plot(abs(dFrames-gnuradioData));
title('Comparison gnuradio calculated FRAMES');
ylim([-1 1]);


