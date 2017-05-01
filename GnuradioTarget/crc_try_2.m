function Output = crc_try_2(Input)
    GenPoly = de2bi(79764919); 
    %G(x)=x32+x26+x23+x22+x16+x12+ x11+ x10+ x8+ x7+ x5+ x4+ x2+x+1
    BufferInit = zeros(1,32);
    Input = [ Input zeros(1,32)];
        for i = 1:length(Input)
            temp1 = BufferInit(end);
            temp2 = temp1*GenPoly;
            for j = length(BufferInit):-1:2
                BufferInit(j) = xor(temp2(j), BufferInit(j-1));
            end
            BufferInit(1) = xor(Input(i), temp2(1));
        end
    Output = fliplr(BufferInit)
%     Output = BufferInit;
end