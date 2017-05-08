%
%   LuboJ.
%
clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input parameters.
%

fftLen = 64;
packetLen = 96;
cpLen = fftLen/4;
nProcessPackets = 10;
sync1 = [0., 0., 0., 0., 0., 0., 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 0., 0., 0., 0., 0.];
sync2 = [0, 0, 0, 0, 0, 0, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, 0, 1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 0, 0, 0, 0, 0];
% sync1 = [sync1 sync1];
% sync2 = [sync2 sync2];

%
%   It has to be this way, same order as in python don't change it,
%   otherwise it's not giving right results!
%
%(range(-26,-21) + range(-20,-7) + range(-6, 0) + range(1,7) + range(8,21) + range(22,27),)
%[39:43 45:57 59:64 2:7 9:21 23:27]
occupiedCarriers = [ ...
    (fftLen+1-26):(fftLen-21) ...
    (fftLen+1-20):(fftLen-7) ...
    (fftLen+1-6):(fftLen-0) ...
    1+1:7 ...
    1+8:21 ...
    1+22:27 ...
];
% pilotCarriers = (-51,-18,-15,16,19);
pilotCarriers = [ ...
    (fftLen+1-21) ...
    (fftLen+1-7) ...
    1+7          ...
    1+21          ...
];
pilotSymbols = [1 1 1 -1];

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
    headerLen = 8*floor(length(occupiedCarriers)/8);
    header = [header zeros(1,headerLen-length(header))];
    header = header(1:headerLen);
    
    payload = de2bi(generatePayload(dataIn), 8)';
    payload = reshape(payload,1,numel(payload));

    modHeader = step(headMod, header')';    
    modPayload = step(payloadMod, payload')';
    packet = [modHeader modPayload];
    dPackets = [dPackets packet];
end
packetLen = length(packet);

% mySig = dPackets;

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
packetLen = length(dFrames)/nProcessPackets;

benchmark = toc

% mySig = dFrames;

ifftSig = []
for k = 1:fftLen:packetLen*nProcessPackets
    ifftChunk = ifft(ifftshift(dFrames(k:k+fftLen-1))).*fftLen;
    ifftSig = [ifftSig ifftChunk(end-cpLen+1:end) ifftChunk];
end

mySig = ifftSig;

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

