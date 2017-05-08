%
%   LuboJ.
%
close all
clear all

nRead = 5000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f2 = fopen('GnuradioTarget/block_tests_files/ofdm_sig_outGUI.txt','r');
if (f2 == -1)
    f2 = fopen('block_tests_files/ofdm_sig_outGUI.txt','r');
end
mySig = arrayToComplex(fread(f2,2*nRead,'float32')');
fclose(f2);

%mySig = 1i*real(mySig) + imag(mySig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f2 = fopen('GnuradioTarget/block_tests_files/ofdm_sig_outScript.txt','r');
% if (f2 == -1)
%     f2 = fopen('block_tests_files/ofdm_sig_outScript.txt','r');
% end
% mySig = arrayToComplex(fread(f2,nRead,'float32')');
% fclose(f2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nRead > length(mySig))
    nRead = length(mySig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = fopen('GnuradioTarget/block_tests_files/ofdm_sig_outScript.txt','r');
if (f == -1)
    f = fopen('block_tests_files/ofdm_sig_outScript.txt','r');
end
gnuradioSamp = arrayToComplex(fread(f,2*nRead,'float32')');
fclose(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f = fopen('block_tests_files/ofdm_sig_out.txt','r');
% if (f == -1)
%     f = fopen('GnuradioTarget/block_tests_files/ofdm_sig_out.txt','r');
% end
% gnuradioSamp = arrayToComplex(fread(f,2*nRead,'float32')');    % <----- pozor normovanie!!
% fclose(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
plot(abs(mySig-gnuradioSamp));
title('Gnuradio vs calculated comparison');
% ylim([-1 1]);

