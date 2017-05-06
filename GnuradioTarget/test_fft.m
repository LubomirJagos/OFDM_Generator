%
%   LuboJ.
%

%
%   Processing
%
% dataIn = [0:63];
% dataIn = zeros(1,64); dataIn(16) = 64;
dataIn = [0,8-7j,0,0,1+5j,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12-3j,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
dataOut = ifft(ifftshift(dataIn));

figure;
plot(real(dataOut));
title('Calculated FFT');
hold on;
plot(imag(dataOut));
hold off;

%
%   Check
%
f = fopen('block_tests_files/test_fft.txt','r');
samples = arrayToComplex(fread(f,inf,'float32')').*1/64;    % <----- pozor normovanie!!
fclose(f);

figure;
plot(real(samples));
title('Gnuradio FFT64 output, 0:63');
hold on;
plot(imag(samples));
hold off;
