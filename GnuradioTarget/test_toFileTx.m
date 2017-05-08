%
%   LuboJ.
%
clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input parameters.
%
fftLen = 64;
cpLen = fftLen/4;
packetLen = 96;
nProcessPackets = 200;

%
%   It has to be this way, same order as in python don't change it,
%   otherwise it's not giving right results!
%
%(range(-26,-21) + range(-20,-7) + range(-6,0) + range(1,7) + range(8,21) + range(22,27),)
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

sync1 = [0., 0., 0., 0., 0., 0., 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 0., 0., 0., 0., 0.];
sync2 = [0, 0, 0, 0, 0, 0, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, 0, 1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 0, 0, 0, 0, 0];

% sync1 = [sync1 sync1];
% sync2 = [sync2 sync2];

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

    if (rem(j,10) == 1)
        txStr = uint8(' Aby sa mi neustale neopakovali spravy a mal to take pestrejsie vysielam 10 sprav.');
    elseif (rem(j,10) == 2)
        txStr = uint8(' Toto je len 1. sprava aby som videl ci to funguje :D Uz som to skusil vypnut a zapnut :D');
    elseif (rem(j,10) == 3)
        txStr = uint8(' Dalsia sprava. Som rad ze som ju prijal. Sprava cislo dva sa ozyva :D Podme na jedno.');
    elseif (rem(j,10) == 4)
        txStr = uint8(' Dufam ze za toto a spravene statnice dostanem titul. Som rad ze to aspon funguje.');
    elseif (rem(j,10) == 5)
        txStr = uint8(' Choose life, fucking television. Trainspotting 2 is really great movie. Watch it!');
    elseif (rem(j,10) == 6)
        txStr = uint8(' Slabucke vysielanie. Skusme nejaku presmycku. Zlutoucky kun upenlive prosikal v noci.');
    elseif (rem(j,10) == 7)
        txStr = uint8(' Ach jaj tyoto pisanie je zabavnejsie nez to kodenie, ale malo by to fungovat. Dufam ze hej/.');
    elseif (rem(j,10) == 8)
        txStr = uint8(' I choose not to choose life. And reasons? There are no reasons when you gpot heroin.');
    elseif (rem(j,10) == 9)
        txStr = uint8(' Tak toto je zvysok 9. Zvysok 9 je posledny :D Potom to zacne cele cyklicky odznova :D');
    else
        txStr = uint8(' Tak zvysok 0. Toz Konecna. Viavt UREL 2017! Hurray, good bye by LuboJ. >D');
    end

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
  'Gain', 20    ...
);

for counter = 1:1000
  disp('Transmitting seq. num. ');
  counter
  step(Tx, ifftSig');
end

release(Tx);

