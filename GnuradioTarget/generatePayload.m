function packet = generatePayload(payload)
    crc = gnuradioCRC32(payload);
    packet = [                                ...
        payload                               ...
        uint8(bitand(crc,255))                ...
        uint8(bitand(bitshift(crc,-8),255))   ...
        uint8(bitand(bitshift(crc,-16),255))  ...
        uint8(bitand(bitshift(crc,-24),255))  ...
    ];
end