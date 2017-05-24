%
%   LuboJ.
%
clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input parameters.
%

ifftLen = 256;
packetLen = 420;
cpLen = ifftLen/4;
nProcessPackets = 1;
sync1 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,1.4142,0,1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,1.4142,0,1.4142,0,1.4142,0,1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,1.4142,0,-1.4142,0,1.4142,0,1.4142,0,1.4142,0,1.4142,0,0,0,1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,1.4142,0,1.4142,0,1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,0,0,-1.4142,0,-1.4142,0,1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,1.4142,0,1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
sync2 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,0,-1,0,-1,0,1,0,1,0,1,0,-1,0,-1,0,-1,0,1,0,-1,0,1,0,1,0,-1,0,-1,0,-1,0,-1,0,-1,0,1,0,-1,0,-1,0,1,0,-1,0,-1,0,1,0,-1,0,1,0,1,0,1,0,-1,0,1,0,-1,0,0,0,1,0,1,0,1,0,-1,0,1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,1,0,-1,0,1,0,1,0,1,0,-1,0,1,0,1,0,1,0,-1,0,1,0,1,0,1,0,1,0,1,0,-1,0,-1,0,1,0,-1,0,-1,0,-1,0,-1,0,1,0,1,0,0,0,-1,0,-1,0,1,0,-1,0,1,0,1,0,-1,0,-1,0,-1,0,1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,1,0,1,0,-1,0,1,0,-1,0,-1,0,1,0,1,0,-1,0,1,0,-1,0,1,0,1,0,1,0,1,0,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
occupiedCarriers = [-102:-76 -74:-40 -38:-12 -10:-1 1:10 12:38 40:74 76:102];
pilotCarriers = [-103 -75 -39 -11 11 39 75 103];
pilotSymbols = [1,1,1,1,1,1,1,-1];

for k = 1:length(occupiedCarriers)
    if (occupiedCarriers(k) < 0)
        occupiedCarriers(k) = occupiedCarriers(k) + ifftLen + 1;
    else
        occupiedCarriers(k) = occupiedCarriers(k) + 1;
    end
end
for k = 1:length(pilotCarriers)
    if (pilotCarriers(k) < 0)
        pilotCarriers(k) = pilotCarriers(k) + ifftLen + 1;
    else
        pilotCarriers(k) = pilotCarriers(k) + 1;
    end
end

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
if (f == -1)
    f = fopen('ofdm_sig_in.txt','r');
end

dataIn = [];
dPackets = [];
dFrames = [];
for j = 1:nProcessPackets
    dataIn = fread(f,packetLen,'uint8')';

    header = generateHeader(packetLen+4,j-1);       %plus 4 because there is CRC added in payload    
   headerLen = 8*floor(length(occupiedCarriers)/8);   %for 128 and 64
%    it's floor, otherwise ceil, but more than 128 is not running anyway
%     headerLen = 8*ceil(length(occupiedCarriers)/8);
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
    ifftLen,                         ...
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
for k = 1:ifftLen:packetLen*nProcessPackets
    ifftChunk = ifft(ifftshift(dFrames(k:k+ifftLen-1))).*ifftLen;
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

