%
%   LuboJ.
%
clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input parameters.
%
fs = 200e3;
fftLen = 64;
cpLen = 16;
packetLen = 96;
nProcessPackets = 100;
occupiedCarriers = [39:43 45:57 59:64 2:7 9:21 23:27];  % <-------- It has to be this way, don't change order, otherwise it's not giving right results!
pilotCarriers = [44 58 8 22];
pilotSymbols = [1 1 1 -1];
sync1 = [0., 0., 0., 0., 0., 0., 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 0., 0., 0., 0., 0.];
sync2 = [0, 0, 0, 0, 0, 0, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, 0, 1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 0, 0, 0, 0, 0];

nSymbols = ceil(packetLen/length(occupiedCarriers));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Create modulators.
%
headMod = comm.BPSKModulator('PhaseOffset', pi);
payloadMod = comm.QPSKModulator('PhaseOffset', 3/4*pi, 'BitInput', true);

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Generovanie packetov
%
f = fopen('block_tests_files/test_cp_In.txt','r');
dataIn = [];
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
end

fclose(f);
release(headMod);
release(payloadMod);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   POZOR MODULACIOU SA MENI DLZKA PAKETOV!
%
%   Gnuradio seka packety po novej dlzke packetov: 
%       (length(modHeader)+length(modPayload))
%
%   potom vlozi synchronizacnu postupnost
%
dFrames = allocateCarriers(         ...
    dPackets,                       ...
    fftLen,                         ...
    length(packet),                 ...
    nProcessPackets,                ...
    occupiedCarriers,               ...
    pilotCarriers,                  ...
    pilotSymbols,                   ...
    sync1,                          ...
    sync2);

benchmark = toc

ifftSig = []
for k = 1:fftLen:fftLen*nProcessPackets
    ifftChunk = ifft(ifftshift(dFrames(k:k+fftLen-1))).*fftLen;
    ifftSig = [ifftSig ifftChunk(end-cpLen+1:end) ifftChunk];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Results comparison with Gnuradio %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = fopen('block_tests_files/test_cp_Out.txt','r');
gnuradioSamp = arrayToComplex(fread(f,2*length(ifftSig),'float32')');    % <----- pozor normovanie!!
fclose(f);

figure;
plot(real(ifftSig));
title('Calculated IFFT signal with cyclic prefix');
hold on;
plot(imag(ifftSig));
hold off;

figure;
plot(real(gnuradioSamp));
title('Gnuradio CP block output');
hold on;
plot(imag(gnuradioSamp));
hold off;

figure;
plot(abs(ifftSig-gnuradioSamp));
title('Gnuradio vs calculated comparison');
% ylim([-1 1]);

