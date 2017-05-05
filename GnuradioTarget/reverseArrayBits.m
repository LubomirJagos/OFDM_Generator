function y = reverseArrayBits(A, padTo)
    A = de2bi(A);
    if (length(A) < padTo)
        A = [A zeros(1,padTo-length(A))];
    end
    y = bi2de(fliplr(A));
end