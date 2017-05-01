nPackets = 20;
nLen = nPackets*100;

% f = fopen('_crc_seqRand_crc_rand_in_1.txt','r');
f = fopen('_crc_const_seq_in_1.txt','r');
notag = fread(f, nLen, 'uint8')';
fclose(f);

% f2 = fopen('_crc_seqRand_crc_rand_out_1.txt','r');
f2 = fopen('_crc_seq_const_out_1.txt','r');
tag = uint32(fread(f2, nLen, 'uint8'))';
fclose(f2);

notag2 = [];
crcplot = [];
for i = 1:96:(nPackets)*96
    crc = gnuradioCRC2(notag(i:i+96-1));
    crcAux = de2bi(crc);
    crcAux = [crcAux zeros(1,32-length(crcAux))];
    crc1 = bi2de(crcAux(1:8));
    crc2 = bi2de(crcAux(9:16));
    crc3 = bi2de(crcAux(17:24));
    crc4 = bi2de(crcAux(25:32));
    crcplot = [crcplot zeros(1,96) crc1 crc2 crc3 crc4];
    notag2 = [notag2 notag(i:i+96-1) crc1 crc2 crc3 crc4];
end

figure;
plot(tag,'-o');
hold on;
plot(notag2,'-xr');
%plot(crcplot,'-xg');
hold off;

figure;
plot(tag-notag2)


