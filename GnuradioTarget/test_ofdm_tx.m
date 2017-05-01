% LuboJ.

f = fopen('_test_ofdm_tx.txt','r');
sig = fread(f,20e3,'float32');
fclose(f);
sig = reshape(sig,2,length(sig)/2);
sig = sig(1,:) + 1i*sig(2,:);

figure;
plot(real(sig));
title('OFDM signal z gnuradia');
hold on;
plot(imag(sig),'-r');
hold off;

figure;
plot(10*log10(abs(ifftshift(fft(ifftshift(sig(1:255)))))));
title('OFDM signal spektrum');
for k = 1:128:100*128
    plot(10*log10(abs(ifftshift(fft(ifftshift(sig(k:k+255)))))));
    pause(0.1);
end





