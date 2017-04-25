% Bc. Lubomir Jagos, 5.4.2017
%
% Streamovanie dat do USRP.
%

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

% % c = parcluster
% % j = createJob(c);

% hMod = comm.DPSKModulator('BitInput',true);
% hMod = comm.GeneralQAMModulator;   
disp('Creating OFDM modulator');
hMod = comm.OFDMModulator;
%set(hMod,'FFTLength',128);
%set(hMod,'FFTLength',1024);
set(hMod,'FFTLength',2048);
%set(hMod,'FFTLength',4096);
%set(hMod,'FFTLength',8192);
%set(hMod,'FFTLength',16384);

disp('Creating USRP TX');
interpolation = 200;
Tx = comm.SDRuTransmitter(...
              'Platform','N200/N210/USRP2', ...
              'IPAddress', '192.168.10.2',  ...
              'CenterFrequency', 8e6,       ...
              'InterpolationFactor', interpolation)
%set(Tx,'CenterFrequency',8e6);
set(Tx,'EnableBurstMode',false);
%set(Tx,'UnderrunOutputPort',true);

%pre OFDM modulator
set(hMod,'NumGuardBandCarriers',[0;0]);
set(hMod,'CyclicPrefixLength',0);

guardBandLen = get(hMod,'NumGuardBandCarriers');
fftLen = get(hMod, 'FFTLength') - guardBandLen(1) - guardBandLen(2);

oversample = 2;
fs = 100e6/interpolation;
filtFc = 300e3;
[b,a] = butter(3,filtFc/fs);
filtSignal = [];

disp('Start generate sequence.');
k=[];
sig = [];
counter = 0;
    data = reshape([ones(1,256);zeros(1,256);zeros(1,256);zeros(1,256);zeros(1,256);zeros(1,256);zeros(1,256);zeros(1,256);],2048,1);
while counter < 5000
    display(horzcat('Generate seq ', num2str(counter)));
    data = randi([0 1], fftLen, 1);
    modSignal = step(hMod, data);
%     filtSignal = filter(b,a,modSignal);
%     filtSignal = resample(step(hMod, data),oversample,1);
    %modSignal = upsample(step(hMod, data),oversample); %nepouzivat, lebo
                                                        %kopiurje cyklicky signal
%     modSignal = awgn(modSignal,20,0);
    %modSignal = awgn(step(hMod, data),20,0);

%     step(Tx, filtSignal);
    step(Tx, modSignal);    
    counter = counter + 1;
end

release(Tx);



