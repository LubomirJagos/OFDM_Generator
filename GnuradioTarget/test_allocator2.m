%
%   LuboJ.
%
%   clear start
%
clear all; close all;

fftLen = 64;
packetLen = 96;
nProcessPackets = 5;
occupiedCarriers = [39:43 45:57 59:64 2:7 9:21 23:27];  % <-------- It has to be this way, don't change order, otherwise it's not giving right results!
pilotCarriers = [44 58 8 22];
pilotSymbols = [10 10 10 -10];
sync1 = [0., 0., 0., 0., 0., 0., 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., -1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., -1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 1.41421356, 0., 0., 0., 0., 0., 0.];
sync2 = [0, 0, 0, 0, 0, 0, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, 0, 1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 0, 0, 0, 0, 0];

f = fopen('block_tests_files/test_allocator_96pLen_In.txt','r');
data = fread(f, 2*500, 'float32');
data = arrayToComplex(data');
fclose(f);

myData = allocateCarriers(          ...
    data,                           ...
    fftLen,                         ...
    packetLen,                      ...
    nProcessPackets,                ...
    occupiedCarriers,               ...
    pilotCarriers,                  ...
    pilotSymbols,                   ...
    sync1,                          ...
    sync2);


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
