function header = generateHeader(lenPacket, numPacket)
    %masking packet length 0x0FFF
    %dlzka paketu je ulozena v 12b cisle, tak toto je maskovanie horneho
    %bajtu
    nPacketForCRC = bitand(bitshift(lenPacket,8)+bitshift(lenPacket,-8),65535); %dlzka paketu
    %poradie paketu je tiez ulozene v 12b tak toto je maskovanie horneho
    %bajtu
    numPacketForCRC = bitand(bitshift(numPacket,8) + bitshift(numPacket,-8),65535); %poradie paketu

    %bitove usporiadanie v subore je LittleEndian, preto treba
    %preusporiadat jednotlive bajty dlzky paketu a jeho poradia
    crc8 = gnuradioCRC8([            ...
        bitshift(nPacketForCRC,-8)   ...
        bitand(nPacketForCRC,255)    ...
        bitshift(numPacketForCRC,-8) ...
        bitand(numPacketForCRC,255)  ...
    ]);

    nPacket12b = de2bi(lenPacket);
    if (length(nPacket12b) < 12)
        nPacket12b = [nPacket12b zeros(1,12-length(nPacket12b))];
    else
        nPacket12b = nPacket12b(1:12);
    end

    numPacket12b = de2bi(numPacket);
    if (length(numPacket12b) < 12)
        numPacket12b = [numPacket12b zeros(1,12-length(numPacket12b))];
    else
        numPacket12b = numPacket12b(1:12);
    end

%     header = [       ...
%         nPacket12b   ...
%         numPacket12b ...
%         de2bi(crc8)  ...
%     ];
%     header = [header zeros(1,48-length(header))];


    %zostavenie hlavicky:
    %        |dlzka paketu 12b| poradie paketu 12b|CRC8|
    %  tu je to ulozene uz predpripravene pre BPSK modulator, tj. rozsekane
    %  na jednicky a nuly 'Repack Bits'
    header = [ ...
        nPacket12b(1:8)    ...
        nPacket12b(9:12)   ...
        numPacket12b(1:4)  ...
        numPacket12b(5:12) ...
        de2bi(crc8,8)      ...
    ];
end

