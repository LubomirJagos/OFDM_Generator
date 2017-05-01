%
%   Citanie dat z binarnych suborov z gnuradia.
%
%   Dlzka paketu 96B vzoriek = 64B data + 32B crc
%
close all; clear all;

%% Citanie hlavickovych bitov
nHeader = 100;
nPacket = 96;

f = fopen('_crc_seq_const_out_1.txt','r');
headerBits = fread(f,nHeader,'uint8')';
fclose(f);

crc = gnuradioCRC(headerBits(1:nPacket))
crcAux = de2bi(crc)
crc1 = bi2de(crcAux(1:8));
crc2 = bi2de(crcAux(9:16));
crc3 = bi2de(crcAux(17:24));
crc4 = bi2de(crcAux(25:end));

figure;
subplot(311); stem(headerBits);
title('Hlavicka paketov');
subplot(312); stem(headerBits);
title('1. paket');
subplot(313); stem([headerBits(1:nPacket) crc1 crc2 crc3 crc4]);
title('Spocitany paket');

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
