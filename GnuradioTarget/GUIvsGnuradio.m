%
%   LuboJ.
%

f2 = fopen('block_tests_files/ofdm_sig_outGUI.txt','r');
mySig = arrayToComplex(fread(f2,420,'float32')');    % <----- pozor normovanie!!
fclose(f2);

f = fopen('block_tests_files/ofdm_sig_out.txt','r');
gnuradioSamp = arrayToComplex(fread(f,2*length(mySig),'float32')');    % <----- pozor normovanie!!
fclose(f);

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
plot(abs(mySig-gnuradioSamp));
title('Gnuradio vs calculated comparison');
% ylim([-1 1]);

