%
%   LuboJ.
%
clear all; close all;

%
%   Create modulators.
%
headMod = comm.BPSKModulator('PhaseOffset', pi);
payloadMod = comm.QPSKModulator('PhaseOffset', 3/4*pi, 'BitInput', true);

%
%   Generate data.
%
sync1 = [0., 0., 0., 0., 0., 0., 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 0., 0., 0., 0., 0.];
sync2 = [0, 0, 0, 0, 0, 0, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, 0, 1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 0, 0, 0, 0, 0];
occupiedCarriers1 = [1:6 8:20 22:26]+1;
occupiedCarriers2 = [-26:-22 -20:-8 -6:-1]+1;
occupiedCarriers = [occupiedCarriers1 occupiedCarriers2]+32;
pilotCarriers = [7 21 43 57]+1;
pilotSymbols = [1 1 1 -1];          %amplitudes

fs = 200e3;
fftLen = 64;            %must be even!
packetLen = 96;
nSymbols = ceil(packetLen/length(occupiedCarriers));
nProcessPackets = 500;

f = fopen('_matlab_muxstream_for_gnuradio.txt','w');

for j = 1:nProcessPackets
    dataIn = randi(255,1,packetLen);

    header = generateHeader(packetLen+4,j-1);       %plus 4 because there is CRC added in payload
    payload = de2bi(generatePayload(dataIn), 8)';
    payload = reshape(payload,1,numel(payload));

    modHeader = step(headMod, header')';
    modPayload = step(payloadMod, payload')';
    packet = [uint8(modHeader) uint8(modPayload)];
    fwrite(f, packet);
end
fclose(f);
