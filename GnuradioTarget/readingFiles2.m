%
%   Citanie dat z binarnych suborov z gnuradia.
%
%   Dlzka paketu 96B vzoriek = 64B data + 32B crc
%
close all; clear all;

%% Citanie hlavickovych bitov
nPacket = 96;
nPacketProcess = 80;
nRead = (nPacket+4)*nPacketProcess;

f = fopen('_crc_seqRand_crc_rand_out_1.txt','r');
headerBytes = fread(f,nRead,'uint8')';
fclose(f);

calcPackets = [];
for i = 1:(nPacket+4):(nPacket+4)*nPacketProcess
    crc = gnuradioCRC(headerBytes(i:i+nPacket-1));
    crcAux = de2bi(crc);
    crc1 = bi2de(crcAux(1:8));
    crc2 = bi2de(crcAux(9:16));
    crc3 = bi2de(crcAux(17:24));
    crc4 = bi2de([crcAux(25:end) zeros(1,32-length(crcAux))]);
    calcPackets = [calcPackets headerBytes(i:i+nPacket-1) crc1 crc2 crc3 crc4];
end

figure;
subplot(211); stem(headerBytes);
hold on;
stem(calcPackets,'-xr');
hold off;
title('Hlavicka paketov');
subplot(212);
plot(headerBytes - calcPackets);


%% Citanie vystupneho signalu
% 
% nSig = 100e3;
% 
% f2 = fopen('_ofdm_signal_tx.txt','r');
% fseek(f2,1e3,'bof');
% signalOut = fread(f2,nSig,'float32');       %komplexny signal 4B realna cast, 4B imaginarna cast
% signalOut2 = reshape(signalOut, 2, length(signalOut)/2);
% fclose(f2);
% 
% figure;
% plot(signalOut2(1,1:end));
% title('Output signal from GnuRadio');
% grid on;
% hold on;
% plot(signalOut2(2,1:end),'-r');
% hold off;
% 
%% Citanie vstupnej sekvencie
% 
% nSeq = nHeader;
% 
% f = fopen('_ofdm_rand_seq.txt','r');
% inputSeq = fread(f,nSeq,'uint8');
% fclose(f);
% 
% figure;
% stem(inputSeq);
% title('Vstupna sekvencia bitov');
