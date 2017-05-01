% LuboJ.

f = fopen('_test_ofdm_taggedStream.txt','r');
sig = fread(f,300,'float32');
fclose(f);
sig = reshape(sig,2,length(sig)/2);
sig = sig(1,:) + 1i*sig(2,:);

figure;
plot(real(sig));
title('OFDM znackovany stream dat - gnuradio');
hold on;
plot(imag(sig),'-r');
hold off;
