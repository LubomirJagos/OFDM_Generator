nPackets = 80;
nLen = (nPackets+1)*104;

f = fopen('_ofdm_rand_seq.txt','r');
noHeader = fread(f, nLen, 'uint8')';
fclose(f);

f2 = fopen('_ofdm_header_bits.txt','r');
header = uint32(fread(f2, nLen*8, 'uint8'))';
fclose(f2);

headerBytes = [];
for k = 1:104*8:nLen*8
    headerPacketLength = bi2de(header(k:k+11));
    headerNumber = bi2de(header(k+12:k+23));
    headerCRC = bi2de(header(k+24:k+32));
    headerData = [];
    for j = 33:8:104        
        headerData = bi2de(header(j:j+7));
    end
    headerBytes = [headerBytes headerPacketLength headerNumber headerCRC headerData];
end

out = [];
for k = 1:100:nPackets*100
    out = [out noHeader(k:k+101) 0 0 0];
end

out2 = [];
for j = 1:8:104*8*10        
    out2 = [out2 bi2de(header(j:j+7))];
end



figure;
plot(headerBytes, '-o');
hold on;
plot(out, '-xr');
hold off;

figure;
stem(out2);