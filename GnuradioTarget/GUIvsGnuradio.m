%
%   LuboJ.
%
% close all
% clear all

nRead = 5000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f2 = fopen('GnuradioTarget/block_tests_files/ofdm_sig_outGUI.txt','r');
if (f2 == -1)
    f2 = fopen('block_tests_files/ofdm_sig_outGUI.txt','r');
end
mySig = arrayToComplex(fread(f2,2*nRead,'float32')');
fclose(f2);

% mySig = 1i*real(mySig) + imag(mySig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nRead > length(mySig))
    nRead = length(mySig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = fopen('GnuradioTarget/block_tests_files/ofdm_sig_out.txt','r');
if (f1 == -1)
    f1 = fopen('block_tests_files/ofdm_sig_out.txt','r');
end
gnuradioSamp = arrayToComplex(fread(f1,2*nRead,'float32')');
fclose(f1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% index = 10;
% mySig = mySig*gnuradioSamp(index)/mySig(index);

figure;
plot(real(mySig));
title('Calculated signal');
hold on;
plot(imag(mySig));
hold off;

figure;
plot(real(gnuradioSamp));
title('Gnuradio output signals');
hold on;
plot(imag(gnuradioSamp));
hold off;

figure;
subplot(211);
plot(real(mySig) - real(gnuradioSamp));
title('Gnuradio vs calculated real comparison');
subplot(212);
plot(imag(mySig) - imag(gnuradioSamp));
title('Gnuradio vs calculated imag comparison');
% ylim([-1 1]);

% frames = ifft(ifftshift(gnuradioSamp));
% figure;
% plot(real(frames));
% hold on;
% plot(imag(frames), '-r');
% hold off;