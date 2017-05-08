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
nProcessPackets = 10;
f = fopen('block_tests_files/ofdm_sig_in.txt','r');
if (f == -1)
    f = fopen('GnuradioTarget/block_tests_files/ofdm_sig_in.txt','r');
end

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
dataIn = [];
dPackets = [];
dFrames = [];

tic;    %benchmark
for j = 1:nProcessPackets
    dataIn = fread(f,packetLen,'uint8')';
    
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

figure;
plot(real(dPackets));
hold on;
plot(imag(dPackets),'-r');
hold off;

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
for k = 1:fftLen:length(dFrames-fftLen+1)
    ifftChunk = ifft(ifftshift(dFrames(k:k+fftLen-1))) .* fftLen;
    ifftSig = [ifftSig ifftChunk(end-cpLen+1:end) ifftChunk];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Write to file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ifftSig = 1i.*real(ifftSig) + imag(ifftSig);
ifftSig = ifftSig .* 0.05;

f = fopen('block_tests_files/ofdm_sig_outScript.txt','w');
if (f == -1)
    f = fopen('GnuradioTarget/block_tests_files/ofdm_sig_outScript.txt','w');
end

dSigReal = real(ifftSig)';
dSigImag = imag(ifftSig)';
fwrite(f, reshape([dSigReal dSigImag]',1,2*length(dSigReal)), 'float32');

fclose(f);





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

