%
%   LuboJ. testing protocol formatter
%

%
%   Input values = values for processing
%

crcLen = 4;             %CHECK IF THERE IS 'Stream CRC32' block used before

nPacket = 96 + crcLen;
nPacket = bitand(uint16(nPacket),hex2dec('0FFF'));  %just 12 bits
nPacketProcess = 1;
nRead = nPacketProcess*nPacket;
nHeader = 48;                       %48Bytes (ie. [0 0 1 1 1 0])
numPacket = 575;
numPacket = bitand(uint16(numPacket),hex2dec('0FFF'));  %just 12 bits

%
%   Output values = target values
%

f2 = fopen('testProtocolFormatterOut.txt','r');
fseek(f2,numPacket*nHeader,'bof');                  %move in file for concrete packet position
headerBytes = fread(f2,nHeader,'uint8')';
fclose(f2);

disp('Protocol Formatter Output, HEX values');
packetLen = dec2hex(bi2de(headerBytes(1:12)))
packetNum = dec2hex(bi2de(headerBytes(13:24)))
packetCRC = dec2hex(bi2de(headerBytes(25:32)))
%rest 33:48 are just zeros

%
%   PROCESSING INPUT VALUES
%

disp('Processed values');
% numPacketForCRC = swapbytes(numPacket);
crcIn = [];
% if (nPacket > 255)
    nPacketForCRC = bitand(bitshift(nPacket,8)+bitshift(nPacket,-8),65535);
% else
%     nPacketForCRC = nPacket;
% end
% if (numPacket > 255)
    numPacketForCRC = bitand(bitshift(numPacket,8) + bitshift(numPacket,-8),65535);
% else
%     numPacketForCRC = numPacket;
% end

crc8 = gnuradioCRC8([            ...
    bitshift(nPacketForCRC,-8)   ...
    bitand(nPacketForCRC,255)    ...
    bitshift(numPacketForCRC,-8) ...
    bitand(numPacketForCRC,255)  ...
]);

    %
    %   Creating header, full header is 48b long
    %
    nPacket12b = de2bi(nPacket,12);
    if (length(nPacket12b) < 12)
        nPacket12b = [nPacket12b zeros(1,12-length(nPacket12b))];
    else
        nPacket12b = nPacket(1:12);
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
    dec2hex(bi2de(header(1:12)))
    dec2hex(bi2de(header(13:24)))
    dec2hex(bi2de(header(25:48)))


