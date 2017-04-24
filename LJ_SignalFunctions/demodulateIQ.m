% IQ modualtion function
%   Takes complex data and perform IQ modulation with carrier fc, each
%   symbol has as many samples as length tc vector.
%   t - time vector for one bit information

function [demodData] = demodulateIQ(rxData, fc, t)
    y = [];
    for(i=1:1:length(rxData)/length(t))   %compute symbols count      
        dataI = rxData((i-1)*length(t)+1:i*length(t)).*cos(2*pi*fc*t); 
        dataQ = rxData((i-1)*length(t)+1:i*length(t)).*-sin(2*pi*fc*t);
        y = horzcat(y,[dataI; dataQ]);
    end
    
    demodData = y;                        %result is matrix with two rows (I and Q)
end
