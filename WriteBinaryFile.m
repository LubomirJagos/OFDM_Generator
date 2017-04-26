%
%   M = count of state
%       QPSK = 4
%       BPSK = 2
%
function y = WriteBinaryFile(M,headLen,dataLen)
    %
    %   Zapis barkrovej postupnosti do binarneho suboru.
    %
    gen = comm.BarkerCode( ...
        'SamplesPerFrame',headLen, ...
        'Length',13, ...
        'OutputDataType', ...
        'int8');
    seq = gen.step()';
    
    f = fopen('OFDMHead.bin','w');
    fwrite(f,seq,'uint8');
    fclose(f);

    %
    %   Write only ones [1 1 ... 1]
    %
    %seq2 = ones(1,dataLen);
    seq2 = randi([0 1],1,dataLen);    
    f2 = fopen('OFDMData.bin','w');
    fwrite(f2,seq2,'uint8');
    fclose(f2);
    
    y = 0;  % no error;
end
