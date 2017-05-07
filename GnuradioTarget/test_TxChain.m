%
%   LuboJ.
%
clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input parameters.
%
fs = 200e3;
fftLen = 128;
cpLen = 16;
packetLen = 96;    
nProcessPackets = 2;

%(range(-26, -21) + range(-20,-7) + range(-6, 0) + range(1,7) + range(8, 21) + range(22, 27),)
%occupiedCarriers = [39:43 45:57 59:64 2:7 9:21 23:27];  % <-------- It has to be this way, don't change order, otherwise it's not giving right results!

%(range(-26, -21) + range(-6, 0) + range(8, 21) + range(22, 27),)
%occupiedCarriers = [39:43 59:64 9:21 23:27];

%(range(-26, -21) + range(8, 21) + range(22, 27),)
% occupiedCarriers = [39:43 9:21 23:27];

%(range(-30, -25) + range(8, 21) + range(24, 30),)
occupiedCarriers = [33:43 45:57 59:64 1:7 9:21 23:100];


pilotCarriers = [44 58 8 22];
pilotSymbols = [1 1 1 -1];
sync1 = [0., 0., 0., 0., 0., 0., 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 0., 0., 0., 0., 0.];
sync2 = [0, 0, 0, 0, 0, 0, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, 0, 1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 0, 0, 0, 0, 0];
% sync1 = [sync1 sync1];
% sync2 = [sync2 sync2]

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
f = fopen('block_tests_files/ofdm_sig_in.txt','r');
dataIn = [];
dPackets = [];
dFrames = [];
for j = 1:nProcessPackets
    dataIn = fread(f,packetLen,'uint8')';
%     dataIn = randi(255,1,packetLen)

    header = generateHeader(packetLen+4,j-1);       %plus 4 because there is CRC added in payload
    payload = de2bi(generatePayload(dataIn), 8)';
    payload = reshape(payload,1,numel(payload));

    modHeader = step(headMod, header')';
    
        %debug LuboJ
        modHeader = [modHeader(1:48) zeros(1,80)];
        %%%%nezabudnut! zmazat
    
    modPayload = step(payloadMod, payload')';
    packet = [modHeader modPayload];
    dPackets = [dPackets packet];
end

mySig = dPackets;

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

f = fopen('block_tests_files/ofdm_sig_out.txt','r');
gnuradioSamp = arrayToComplex(fread(f,2*length(mySig),'float32')');    % <----- pozor normovanie!!
fclose(f);

figure;
plot(real(mySig));
title('Calculated signal');
hold on;
plot(imag(mySig));
hold off;

figure;
plot(real(gnuradioSamp));
title('Gnuradio output signals');
hold on;
plot(imag(gnuradioSamp));
hold off;

figure;
plot(abs(mySig-gnuradioSamp));
title('Gnuradio vs calculated comparison');
% ylim([-1 1]);

