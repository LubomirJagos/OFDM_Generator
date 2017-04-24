%IQ modulation done by multiplying by complex exponential function
function [modData] = modulateIQexp(complexData, fc, t)
    modData = complexData(i)*exp(i*2*pi*fc);
end
