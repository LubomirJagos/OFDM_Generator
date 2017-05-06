%
%   LuboJ.
%
%   clear start
%
clear all; close all;

%
%   Settings.
%
fs = 200e3;
fftLen = 64;
packetLen = 96;
nProcessPackets = 1;

sync1 = [0., 0., 0., 0., 0., 0., 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 0., 0., 0., 0., 0.];
sync2 = [0, 0, 0, 0, 0, 0, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, 0, 1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 0, 0, 0, 0, 0];

occupiedCarriers = [3:64];
% pilotCarriers = [1 2];
% pilotSymbols = [-10 -10];

%pilotCarriers = (range(-26, -21) + range(-20, -7) + range(-6, 0) + range(1, 7) + range(8, 21) + range(22, 27),)
occupiedCarriers = [39:43 45:57 59:64 2:7 9:21 23:27];  % <-------- It has to be this way, don't change order, otherwise it's not giving right results!
pilotCarriers = [44 58 8 22];
pilotSymbols = [10 10 10 -10];

nCarriers = length(occupiedCarriers);
nSymbols = ceil(packetLen/nCarriers);

dataIn = [];
myData = [];
% f = fopen('_meranie_5_5/allocator_in.txt','r');
f = fopen('block_tests_files/test_allocator_96pLen_In.txt','r');

    while (nProcessPackets > 0)
        myData = [myData sync1 sync2];
        needToRead = packetLen;
        while (needToRead > 0)
            if (needToRead > nCarriers)
                [dataIn count] = fread(f, 2*nCarriers, 'float32');
            else
                [dataIn count] = fread(f, 2*needToRead, 'float32');
            end
            dataIn = arrayToComplex(dataIn');
            count = count/2;
            if (count == nCarriers)
                needToRead = needToRead - count;
            else
                dataIn = [dataIn zeros(1,nCarriers-count)];
                needToRead = 0;
            end
            fftFrame = zeros(1,fftLen);
            fftFrame(occupiedCarriers) = dataIn;
            fftFrame(pilotCarriers) = pilotSymbols;
            fftFrame = circshift(fftFrame', floor(fftLen/2))';
            myData = [myData fftFrame];
        end        
        nProcessPackets = nProcessPackets-1;
    end

fclose(f);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Debug checking diff btw. muxed stream in gnuradio and calculated in
%   matlab
%
%
% f2 = fopen('_meranie_5_5/allocator_out.txt','r');
f2 = fopen('block_tests_files/test_allocator_96pLen_Out.txt','r');
gnuradioData = arrayToComplex(fread(f2,2*length(myData),'float32')');
fclose(f2);

figure;
stem(real(gnuradioData));
title('Output allocator, GNURadio');
grid on;
hold on;
stem(imag(gnuradioData),'-r');
hold off;

figure;
stem(real(myData));
title('CALCULATED allocation');
grid on;
hold on;
stem(imag(myData),'-r');
hold off;

figure;
stem(abs(myData-gnuradioData));
title('Comparison gnuradio and calculated');
% ylim([-3 10]);
