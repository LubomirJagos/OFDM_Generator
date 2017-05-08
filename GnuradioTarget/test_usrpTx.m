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
nProcessPackets = 42;

%
%   It has to be this way, same order as in python don't change it,
%   otherwise it's not giving right results!
%
%(range(-60,-52) + range(-50, -20) + range(-10,-9) + range(24, 60),)
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

sync1 = [0., 0., 0., 0., 0., 0., 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 0., 0., 0., 0., 0.];
sync2 = [0, 0, 0, 0, 0, 0, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, 0, 1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 0, 0, 0, 0, 0];

nSymbols = ceil(packetLen/length(occupiedCarriers));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Create modulators.
%
headMod = comm.BPSKModulator('PhaseOffset', pi);
payloadMod = comm.QPSKModulator('PhaseOffset', 3/4*pi, 'BitInput', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Generovanie packetov
%
f = fopen('volacka.wav','r');
dataIn = [];
dPackets = [];
dFrames = [];

tic;    %benchmark
for j = 1:nProcessPackets
%       dataIn = fread(f,packetLen,'uint8')';
%       dataIn = double(dataIn);
      
%    dataIn = randi(255,1,packetLen);

    txStr = [uint8(' Tak si so mnou opakuj, mat je viac nez, Marika Gombitova, strasne pekna piesen. Pustit.')];
    txStr = [txStr 10]; %newline
    txStr = [txStr randi([65 90],1,packetLen-length(txStr))];
    dataIn = double(txStr);
    dataIn = dataIn(1:packetLen);
    
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
packetGenBenchmark = toc;

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
    ifftChunk = ifft(ifftshift(dFrames(k:k+fftLen-1))) .* fftLen;
    ifftSig = [ifftSig ifftChunk(end-cpLen+1:end) ifftChunk];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Transmit generated signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f = fopen('_ofdm_signal_tx.txt');
% ifftSig = arrayToComplex(fread(f,20e3,'float32')');
% fclose(f);

ifftSig = 1i.*real(ifftSig) + imag(ifftSig);
ifftSig = ifftSig .* 0.05;

Tx = comm.SDRuTransmitter(...
  'Platform','N200/N210/USRP2',...
  'IPAddress', '192.168.10.2', ...
  'CenterFrequency', 2.45e9, ...
  'InterpolationFactor', 250,  ...
  'Gain', 30    ...
);

for counter = 1:2000
  disp('Transmitting seq. num. ');
  counter
  step(Tx, ifftSig');
end

release(Tx);

