% IQ modualtion function
%   Takes complex data and perform IQ modulation with carrier fc, each
%   symbol has as many samples as length t vector.
%   t - time vector for one bit information

function [modData] = modulateIQ(complexData, fc, t)
    y = [];
    y_in=[];
    y_qd=[];

    dataI = real(complexData);
    dataQ = imag(complexData);
    for(i=1:length(complexData))
        y1 = dataI(i) * cos(2*pi*fc*t);     % inphase component
        y2 = dataQ(i) * -sin(2*pi*fc*t) ;    % Quadrature component
        y_in = [y_in y1];                     % inphase signal vector
        y_qd = [y_qd y2];                     %quadrature signal vector

        y = [y y1-y2];                        % modulated signal vector (I - Q)
    end
    
    modData = y;
end
