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
fftLen = 64;
packetLen = 96;
nProcessPackets = 5;
occupiedCarriers = [39:43 45:57 59:64 2:7 9:21 23:27];  % <-------- It has to be this way, don't change order, otherwise it's not giving right results!
pilotCarriers = [44 58 8 22];
pilotSymbols = [1 1 1 -1];
sync1 = [0., 0., 0., 0., 0., 0., 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 0., 0., 0., 0., 0.];
sync2 = [0, 0, 0, 0, 0, 0, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, 0, 1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 0, 0, 0, 0, 0];

nSymbols = ceil(packetLen/length(occupiedCarriers));

dataIn = [];
f = fopen('_test_mod_randseq_1.txt','r');

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
%   Gnuradio seka packety po novej dlzke packetov: 
%       (length(modHeader)+length(modPayload))
%
%   potom vlozi synchronizacnu postupnost
%

dFrames = [];
for k = 1:length(packet):nProcessPackets*length(packet)
    frame = allocateCarriers(           ...
        dPackets(k:k+length(packet)-1), ...
        fftLen,                         ...
        packetLen,                      ...
        nProcessPackets,                ...
        occupiedCarriers,               ...
        pilotCarriers,                  ...
        pilotSymbols,                   ...
        sync1,                          ...
        sync2);
    dFrames = [dFrames frame];
end

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


