function y = reversePaddedDecNumBits(A, padTo)
    y = bi2de(fliplr(de2bi(A,padTo)));
end