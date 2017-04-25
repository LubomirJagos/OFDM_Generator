% Bc. Lubomir Jagos, 5.4.2017
%
% Streamovanie dat do USRP.
%

% c = parcluster
% j = createJob(c);

%Rx = comm.SDRuReceiver
%info(Rx)

% Tx = comm.SDRuTransmitter('Platform','N200','IPAddress','192.168.20.2')
% info(Tx)

% OTAZKA Z MathWorks fora
%
% I'm trying to stream from the input to output of an ettus X310 in
% Simulink. The block diagram has just 3 blocks, SDRu receiver,
% buffer, & SDRu Transmitter. The decimation (&interpolation) value
% is 200 points so the sampling rate is only 1 MHz. I maxed the
% frame length in the receiver to 375000. The buffer size
% is 2*375000. I maxed the FastSendDatagramThreshold value
% to 1500 in the windows registry. It's a 1 gbit Ethernet connection
% and a fast PC. Even at 1 MHz I'm getting a lot of underflows (I
% assume from the transmitter) and loss of data. I have dreams of
% going all the way to 200 MHz (with some simple signal processing
% in between.) That may be a pipe dream. Suggestions?? How do I get
% the throughput up?

% hMod = comm.DPSKModulator('BitInput',true);
% hMod = comm.GeneralQAMModulator;   
hMod = comm.OFDMModulator;
set(hMod,'FFTLength',4096);

Rx = comm.SDRuReceiver(...
              'Platform','N200/N210/USRP2', ...
              'IPAddress', '192.168.10.2',  ...
              'CenterFrequency', 3e6,       ...
              'DecimationFactor', 256,      ...
              'FrameLength', 2048)


figure;
for i=1:100 
    disp('Rx sequence: ');
    i
    rxData = step(Rx);
    plot(abs(fft((double(rxData)))));
    pause(0.3);
    %plot(imag(rxData),'r');
end
          
release(Rx);



