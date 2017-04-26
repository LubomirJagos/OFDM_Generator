%
%   M = count of state
%       QPSK = 4
%       BPSK = 2
%
function y = WriteBinaryFile(M)
    %
    %   Zapis barkrovej postupnosti do binarneho suboru.
    %

    gen = comm.BarkerCode( ...
        'SamplesPerFrame',64, ...
        'Length',13, ...
        'OutputDataType', ...
        'int8');
    seq = gen.step();

    seq = [seq repmat([1 0],1,128 - 64/log2(M))];
    
    f = fopen('OFDMHead.bin','w');
    fwrite(f,seq,'uint8');
    fclose(f);

    seq2 = ones(1,88);
    f2 = fopen('OFDMData.bin','w');
    fwrite(f,seq,'uint8');
    fclose(f);
    
    y = 0;  % no error;
end
