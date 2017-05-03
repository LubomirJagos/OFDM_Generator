%
%   Citanie dat z binarnych suborov z gnuradia.
%
%   Dlzka paketu 96B vzoriek = 64B data + 32B crc
%
close all; clear all;

%% Citanie hlavickovych bitov
nPacket = 96;
nPacketProcess = 200;
nRead = (nPacket+4)*nPacketProcess;

f = fopen('gnuradioCRCSeq.txt','r');
headerBytes = fread(f,nRead,'uint8')';
fclose(f);

calcPackets = [];

for i = 1:nPacket:nPacket*nPacketProcess
    crc = gnuradioCRC2(headerBytes(i:i+nPacket-1));
    crcAux = de2bi(crc);
    crc1 = bi2de(crcAux(1:8));
    crc2 = bi2de(crcAux(9:16));
    crc3 = bi2de(crcAux(17:24));
    crc4 = bi2de([crcAux(25:end) zeros(1,32-length(crcAux))]);
    calcPackets = [calcPackets headerBytes(i:i+nPacket-1) crc1 crc2 crc3 crc4];
end

f = fopen('gnuradioCRCseqOUT.txt','r');
header = fread(f,nRead,'uint8')';
fclose(f);

figure;
subplot(211); stem(header);
hold on;
stem(calcPackets,'-xr');
hold off;
title('Hlavicka paketov');
subplot(212);
plot(header - double(calcPackets));

