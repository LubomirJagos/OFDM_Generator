function y = arrayToComplex(inArray)
    inArray = reshape(inArray,2,length(inArray)/2);
    y = inArray(1,:) + 1i*inArray(2,:);
end
