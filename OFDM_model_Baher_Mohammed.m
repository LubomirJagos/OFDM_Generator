%==========================================================================
% The mfile investigates the generation, transmission and reception of
% the OFDM signal without channel noise or HPA effect
% Written By     : Baher Mohamed
%==========================================================================
clear all
clc
close  
%   ---------------
%   A: Setting Parameters
%   ---------------
M = 4;                          %   QPSK signal constellation
no_of_data_points = 64;        %   have 64 data points
block_size = 8;                 %   size of each ofdm block
cp_len = ceil(0.1*block_size);  %   length of cyclic prefix
no_of_ifft_points = block_size;           %   8 points for the FFT/IFFT
no_of_fft_points = block_size;
%   ---------------------------------------------
%   B:  %   +++++   TRANSMITTER    +++++
%   ---------------------------------------------
%   1.  Generate 1 x 64 vector of data points phase representations
data_source = randsrc(1, no_of_data_points, 0:M-1);
figure(1)
stem(data_source); grid on; xlabel('data points'); ylabel('transmitted data phase representation')
title('Transmitted Data "O"')
%   2.  Perform QPSK modulation
qpsk_modulated_data = pskmod(data_source, M);
figure(2)
scatterplot(qpsk_modulated_data);title('qpsk modulated transmitted data')
%   3.  Do IFFT on each block
%   Make the serial stream a matrix where each column represents a pre-OFDM
%   block (w/o cyclic prefixing)
%   First: Find out the number of colums that will exist after reshaping
num_cols=length(qpsk_modulated_data)/block_size;
data_matrix = reshape(qpsk_modulated_data, block_size, num_cols);

%   Second: Create empty matix to put the IFFT'd data
cp_start = block_size-cp_len;
cp_end = block_size;

%   Third: Operate columnwise & do CP
for i=1:num_cols,
    ifft_data_matrix(:,i) = ifft((data_matrix(:,i)),no_of_ifft_points);
    %   Compute and append Cyclic Prefix
    for j=1:cp_len,
       actual_cp(j,i) = ifft_data_matrix(j+cp_start,i);
    end
    %   Append the CP to the existing block to create the actual OFDM block
    ifft_data(:,i) = vertcat(actual_cp(:,i),ifft_data_matrix(:,i));
end

%   4.  Convert to serial stream for transmission
[rows_ifft_data cols_ifft_data]=size(ifft_data);
len_ofdm_data = rows_ifft_data*cols_ifft_data;

%   Actual OFDM signal to be transmitted
ofdm_signal = reshape(ifft_data, 1, len_ofdm_data);
figure(3)
plot(real(ofdm_signal)); xlabel('Time'); ylabel('Amplitude');
title('OFDM Signal');grid on;
%   ------------------------------------------
%   E:  %   +++++   RECEIVER    +++++
%   ------------------------------------------

%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   LubomirJ. analysis OFDM signal
%   ofdm_signal

fftN = 2^10;
spec = abs(fft(abs(ofdm_signal), fftN));
spec2 = abs(fft(real(ofdm_signal), fftN));
figure(4)
stem(spec); grid on; title('OFDM serial signal spectrum');
hold on;
stem(spec2, 'r');

%
%   END analysis
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   1.  Pass the ofdm signal through the channel
recvd_signal = ofdm_signal;

%   4.  Convert Data back to "parallel" form to perform FFT
recvd_signal_matrix = reshape(recvd_signal,rows_ifft_data, cols_ifft_data);

%   5.  Remove CP
recvd_signal_matrix(1:cp_len,:)=[];

%   6.  Perform FFT
for i=1:cols_ifft_data,
    %   FFT
    fft_data_matrix(:,i) = fft(recvd_signal_matrix(:,i),no_of_fft_points);
end

%   7.  Convert to serial stream
recvd_serial_data = reshape(fft_data_matrix, 1,(block_size*num_cols));

%figure(11);
figure(5)
plot(real(recvd_serial_data)); grid on; title('Real RX Data');
figure(6)
plot(imag(recvd_serial_data)); grid on; title('Imag RX Data');
figure(7)
scatterplot(recvd_serial_data); title('RX Data');
%title('My Rcvd Serial Data');

%   8.  Demodulate the data
qpsk_demodulated_data = pskdemod(recvd_serial_data,M);
figure(8)
scatterplot(qpsk_modulated_data);title('qpsk modulated received data')

figure(9)
stem(qpsk_demodulated_data,'rx');
grid on;xlabel('data points');ylabel('received data phase representation');title('Received Data "X"')   




