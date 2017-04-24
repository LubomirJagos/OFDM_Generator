% % This would be QPSK modulation to compare directly from Matlab
% %
% % % Create a 16-PSK modulator System object with bits as inputs and Gray-coded signal constellation
% % hModulator = comm.PSKModulator(4,'BitInput',true);
% % 
% % % Change the phase offset to pi/16 
% % hModulator.PhaseOffset = 0;%pi/16; 
% % 
% % % Modulate and plot the data
% % modData = step(hModulator, txData); 
% % constellation(hModulator);
% % 
% % %constellation points
% % figure(1);
% % scatterplot(modData);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % Create binary data for 24, 4-bit symbols 
% clear all;
% close all;
% 
% %
% %   SW generating samples, setting variables
% %
% 
% %txData = randi([0 1],96,1);
% % Data for QPSK have to be transposed to create column
% txData = [0 0 1 1 0 0 1 1 1 0 0 1 0 1 0 0]';
% 
% % XXXXXXXXXXXXXXXXXXXXXXX QPSK modulatio  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% data_NZR=2*txData-1; % Data Represented at NZR form for QPSK modulation
% s_p_data=reshape(data_NZR,2,length(txData)/2);  % S/P convertion of data
% 
% 
% br=4.2e3; %Let us transmission bit rate  1000000
% f=br; % minimum carrier frequency
% T=1/br; % bit duration
% 
% % compute sample frequency
% samplesPerBit = 99;         % how many samples is for one bit, because bit duration > one sample
% fs = samplesPerBit/T;
% ts = 1/fs;
% 
% t=ts:ts:T; % Time vector for one bit information
% 
% y=[];
% y_in=[];
% y_qd=[];
% for(i=1:length(txData)/2)
%     y1=s_p_data(1,i)*cos(2*pi*f*t); % inphase component
%     y2=s_p_data(2,i)*sin(2*pi*f*t) ;% Quadrature component
%     y_in=[y_in y1]; % inphase signal vector
%     y_qd=[y_qd y2]; %quadrature signal vector
%     y=[y y1+y2]; % modulated signal vector
% end
% 
% Tx_sig=y; % transmitting signal after modulation
% tt=ts:ts:(T*length(txData))/2;
% 
% figure(1);
% 
% subplot(3,1,1);
% plot(tt,y_in,'-x'), grid on;
% title(' wave form for inphase component in QPSK modulation ');
% xlabel('time(sec)');
% ylabel(' amplitude(volt0');
% 
% subplot(3,1,2);
% plot(tt,y_qd,'-x'), grid on;
% title(' wave form for Quadrature component in QPSK modulation ');
% xlabel('time(sec)');
% ylabel(' amplitude(volt0');
% 
% subplot(3,1,3);
% plot(tt,Tx_sig,'r-x'), grid on;
% title('QPSK modulated signal (sum of inphase and Quadrature phase signal)');
% xlabel('time(sec)');
% ylabel(' amplitude(volt0');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %
% %   HW Control
% %
% 
% lecroySampleFreq = 400e3;
% lecroySampleTime = 1/lecroySampleFreq;
% 
% hantekSampleFreqInt = 14;
% hantekSampleFreq = Hantek6022BE_GetSampleFreq(hantekSampleFreqInt);
% hantekSampleTime = 1/hantekSampleFreq;
%     % 0 - 10  = 48MHz
%     % 11      = 16MHz
%     % 12      = 8MHz
%     % 13      = 4MHz
%     % 14 - 24 = 1MHz
%     % 25 = 500kHz
%     % 26 = 200kHz
%     % 27 = 100kHz
% 
% hantekCh1VoltDiv = 6;
% hantekCh2VoltDiv = 6;
%     % 0 = 20mV
%     % 1 = 50mV
%     % 2 = 100mV
%     % 3 = 200mV
%     % 4 = 500mV
%     % 5 = 1V
%     % 6 = 2V
%     % 7 = 5V
% 
% 
% Rx_sig = resample(Tx_sig, lecroySampleFreq, fs);          %signal on generator output
% rxLength = floor(2*(hantekSampleFreq/lecroySampleFreq)*length(Rx_sig));                        %scope received signal, more samples to be sure that I have whole signal
% 
% figure(2);
% subplot(211);
% plot(Rx_sig, '-x'); title('Generator Modulated data resampled');
%     
% lw410 = LW410Interface(1, lecroySampleTime);           %gpibAddr, sampletime
% lw410.wave_data(Rx_sig, 1);
% 
% disp('Waiting for generator');
% pause(2);
% disp('Done.');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   HANTEK reading data
% 
% [ch1Data, ch2Data, triggerPointIndex]= Hantek6022BE_ReadingData(0,0,rxLength,hantekCh1VoltDiv, hantekCh2VoltDiv,hantekSampleFreqInt,0,0,0,0,0,0,0,0);
% rxData = double(ch1Data);
% 
% rxData = resample(double(ch1Data), fs, hantekSampleFreq);
% tc = [0:ts:(length(rxData)-1)*ts];
% 
% figure(2);
% subplot(212);
% plot(rxData, 'r-x');
% title('RX data');
% 
% disp(horzcat('Triggered by point x = ', num2str(triggerPointIndex)));
% 
% 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%   Demodulate data
%

%manual cutting 24.11.2016
%rxData = rxData(21:end);


Rx_data = [];
Rx_sig = rxData; % Received signal
for(i=1:1:(length(txData)+10)/2)         % THIS IN CASE TO SEE
                                           % EXTENDED SYMBOLS
%for(i=1:1:(length(txData))/2)              % txData here as count of transmitted bits! LuboJ.
    
    %%XXXXXX inphase coherent dector XXXXXXX
    Z_in=Rx_sig((i-1)*length(t)+1:i*length(t)).*cos(2*pi*f*t); 
    % above line indicat multiplication of received & inphase carred signal
    
    Z_in_intg=(trapz(t,Z_in))*(2/T);% integration using trapizodial rull
    if(Z_in_intg>0) % Decession Maker
        Rx_in_data=1;
    else
       Rx_in_data=0; 
    end
    
    %%XXXXXX Quadrature coherent dector XXXXXX
    Z_qd=Rx_sig((i-1)*length(t)+1:i*length(t)).*sin(2*pi*f*t);
    %above line indicat multiplication ofreceived & Quadphase carred signal
    
    Z_qd_intg=(trapz(t,Z_qd))*(2/T);%integration using trapizodial rull
        if (Z_qd_intg>0)% Decession Maker
        Rx_qd_data=1;
        else
       Rx_qd_data=0; 
        end
        
        Rx_data=[Rx_data  Rx_in_data  Rx_qd_data]; % Received Data vector
end

figure(3)
subplot(211);
%plot(rxData(21:end), '-x');                     %cutted!
plot(rxData, '-x');
title('Received signal sampling set back to original');

%DRAWING VERTICAL LINE TO DISTINGUISH SYMBOLS
samplesPerBit2 = samplesPerBit;
xVerPoints = [0:samplesPerBit2:(length(txData)-1)*samplesPerBit2];
xVerPoints = [xVerPoints; xVerPoints];
xVerPoints = xVerPoints(:)';
yVerPoints = ylim;
yVerPoints = repmat(yVerPoints,1,length(txData));
%draw lines separately, if it would be one command there would be diagonal
%lines connecting vertical ones!
for i = 1:2:length(xVerPoints)
    hLineObj = line(xVerPoints(i:i+1), yVerPoints(i:i+1));
    set(hLineObj, 'color', 'red');    
end

subplot(212);
stem(Rx_data,'linewidth',3) 
title('Information after Receiveing ');
grid on;

figure(4);
subplot(211);
plot(Tx_sig, '-x');
title('Signal at generator output');

%DRAWING VERTICAL LINE TO DISTINGUISH SYMBOLS
samplesPerBit2 = floor(lecroySampleFreq/fs*samplesPerBit);
xVerPoints = [0:samplesPerBit2:(length(txData)-1)*samplesPerBit2];
xVerPoints = [xVerPoints; xVerPoints];
xVerPoints = xVerPoints(:)';
yVerPoints = ylim;
yVerPoints = repmat(yVerPoints,1,length(txData));
%draw lines separately, if it would be one command there would be diagonal
%lines connecting vertical ones!
for i = 1:2:length(xVerPoints)
    hLineObj = line(xVerPoints(i:i+1), yVerPoints(i:i+1));
    set(hLineObj, 'color', 'red');    
end

subplot(212);
stem(txData);
title('Original bits');

% figure
% scatterplot([Rx_data(1:end/2)' Rx_data(end/2+1:end)']);






% JUST FOR FUN, MOVING VERTICAL LINE OVER PLOT
% figure;
% plot(rxData, 'b-x');
% lineObj = line([0 0], ylim);
% set(lineObj, 'LineWidth', 3, 'Color', 'green');
% for i = 0:4:length(rxData)
%     set(lineObj, 'XData', [i i], 'YData', ylim);
%     pause(0.2);
% end
