function header = generateHeader(lenPacket, numPacket)
    nPacketForCRC = bitand(bitshift(lenPacket,8)+bitshift(lenPacket,-8),65535);
    numPacketForCRC = bitand(bitshift(numPacket,8) + bitshift(numPacket,-8),65535);

    crc8 = gnuradioCRC8([            ...
        bitshift(nPacketForCRC,-8)   ...
        bitand(nPacketForCRC,255)    ...
        bitshift(numPacketForCRC,-8) ...
        bitand(numPacketForCRC,255)  ...
    ]);

    %
    %   Creating header, full header is 48b long
    %
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

    header = [       ...
        nPacket12b   ...
        numPacket12b ...
        de2bi(crc8)  ...
        zeros(1,24)  ...
    ];
end