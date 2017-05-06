close all; clear all;

dataIn = [];
packetLen = 124;

f = fopen('_meranie_5_5/allocator_in.txt','r');

dataIn = fread(f,2*packetLen,'float32')';
dataIn = reshape(dataIn,2,packetLen);
dataIn = dataIn(1,:) + 1i*dataIn(2,:);

fclose(f);

figure;
stem(real(dataIn));
hold on;
stem(imag(dataIn),'-r');
hold off;
