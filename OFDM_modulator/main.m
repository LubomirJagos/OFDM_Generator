function varargout = main(varargin)
% main MATLAB code for main.fig
%      main, by itself, creates a new main or raises the existing
%      singleton*.
%
%      H = main returns the handle to a new main or the handle to
%      the existing singleton*.
%
%      main('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in main.M with the given input arguments.
%
%      main('Property','Value',...) creates a new main or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main

% Last Modified by GUIDE v2.5 01-Jun-2017 08:58:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_OpeningFcn, ...
                   'gui_OutputFcn',  @main_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% LuboJ., add gnuradio functions
addpath('utils');
feature('DefaultCharacterSet','UTF-8');

% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main (see VARARGIN)

% Choose default command line output for main
handles.output = hObject;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   LuboJ. 6.11.2016
%set panels invisible, only default one OFDM.
handles.data.activePanel = 'ofdmWavePanel';
set(handles.ofdmWavePanel,'visible','on');
set(handles.deviceSettingsPanel, 'visible','off');
set(handles.signalOutputPanel, 'visible','off');

handles.lw410 = LW410Interface();
handles.ofdm.dataInputMethod = 'randsrc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%26.4.2017
handles.ofdm.headData = 0;

% Update handles structure
handles.langENG = 1;

% reading images for display output device
handles.usrpDeviceImg = imread('img/usrp.png');
handles.lw410Img = imread('img/lw410.png');

guidata(hObject, handles);

% UIWAIT makes main wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in butSetOutput.
function butSetOutput_Callback(hObject, eventdata, handles)
% hObject    handle to butSetOutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%disp(horzcat('Current Panel: ', handles.data.activePanel));
disp(horzcat('Sequence length #samples: ', num2str(length(handles.data.seq))));

handles.lw410.wave_data(handles.data.seq, 1);

set(handles.statusText, 'String', horzcat('Uploading sequence. ', num2str(length(handles.data.seq)), ' samples'))
guidata(hObject, handles);
    
% --- Executes on button press in ofdmBut.
function ofdmBut_Callback(hObject, eventdata, handles)
% hObject    handle to ofdmBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data.activePanel = 'ofdmWavePanel';
set(handles.ofdmWavePanel,'visible','on')
set(handles.deviceSettingsPanel, 'visible','off');
set(handles.signalOutputPanel, 'visible','off');

guidata(hObject, handles);  %how to obtain data


% --- Executes on button press in deviceSetingsBut.
function deviceSetingsBut_Callback(hObject, eventdata, handles)
% hObject    handle to deviceSetingsBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data.activePanel = 'deviceSettingsPanel';   %set semaphor
set(handles.deviceSettingsPanel, 'visible','on');
set(handles.ofdmWavePanel,'visible','off');
set(handles.signalOutputPanel,'visible','off');

guidata(hObject, handles);  %how to obtain data


% --- Executes on button press in generateOFDMWaveBut.
function generateOFDMWaveBut_Callback(hObject, eventdata, handles)
% hObject    handle to generateOFDMWaveBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data.activePanel = 'signalOutputPanel';   %set semaphor
set(handles.deviceSettingsPanel, 'visible','off');
set(handles.ofdmWavePanel,'visible','off');
set(handles.signalOutputPanel,'visible','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OFDM generating
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% Setting sample rate (default 400kHz)
%--------------------------------------------------------------------------
deviceType = 0;
fs = str2double(get(handles.noDeviceSampleRate, 'string'))*1000;

if (get(handles.usrpDevice, 'Value') == 1)
    fs = str2double(get(handles.USRPSampleRate, 'string'))*1000;
    deviceType = 1;
end;
if (get(handles.lw410Device, 'Value') == 1)
    fs = str2double(get(handles.lw410SampleRate, 'string'))*1000;
    
    contents = cellstr(get(handles.gpibAddrPopup,'String'));
    gpib_addr = str2num(contents{get(handles.gpibAddrPopup,'Value')});

    contents = cellstr(get(handles.sampleFreqPopup,'String'));
    lw410fs = contents{get(handles.sampleFreqPopup,'Value')};
    if (strcmp(lw410fs, '400 MHz'))
        lw410fs = 400e6;
    elseif (strcmp(lw410fs, '40 MHz'))
        lw410fs = 40e6;
    elseif (strcmp(lw410fs, '4 MHz'))
        lw410fs = 4e6;
    elseif (strcmp(lw410fs, '40 kHz'))
        lw410fs = 40e3;
    else
        lw410fs = 40e6;           % DEFAULT VALUE
    end
    deviceType = 2;
end;
fs              %output to command window

%--------------------------------------------------------------------------
% READING DATA FOR GENERATING WAVES FROM GUI by LuboJ.
%--------------------------------------------------------------------------
signalAmplitude = str2double(get(handles.editSignalAmplitude, 'string'))
ifftLen = str2double(get(handles.ifftLen, 'string'))
nSymbols = str2double(get(handles.numOfSymbolsInput, 'string'))
oversampleFactor = str2double(get(handles.oversamplingInput, 'string'))
noDeviceOversampleFactor = str2double(get(handles.noDeviceOversamplingInput, 'string'))
fc = 1e3;
                                                             
insertPrefixInterval= get(handles.intervalRadio, 'Value');
insertGuardInterval= get(handles.guardRadio, 'Value');

guardLength = eval(get(handles.guardLengthInput, 'string'))  %This is array[rows cols]
if (isempty(guardLength))
    guardLeft = [];
    guardRight = [];
else
    guardLeft = zeros(guardLength(1),1);
    guardRight = zeros(guardLength(2),1);
end
cpLen = str2double(get(handles.prefixLengthInput, 'string'))
txGain = str2double(get(handles.txGain, 'string'))

packetLen = str2double(get(handles.packetLen, 'string'))

%
%   Typ modulacie
%
modTypeList = get(handles.modulationTypePopup,'String');
payloadModType = modTypeList{get(handles.modulationTypePopup, 'Value')}
if (strcmp(payloadModType, 'QPSK'))
    M = 4
elseif (strcmp(payloadModType, 'BPSK'))
    M = 2
elseif (strcmp(payloadModType, 'PSK-8'))
    M = 8
end

handles.modM = M;  %setting for GUI

%
%   Nastavenie pilotnych tonov.
%
occupiedCarriers = eval(get(handles.occupiedCarriers, 'string'));
nCarriers = length(occupiedCarriers);       %alias pre lahsi pristup
for k = 1:length(occupiedCarriers)
    if (occupiedCarriers(k) < 0)
        occupiedCarriers(k) = ifftLen + 1 + occupiedCarriers(k);
    else
        occupiedCarriers(k) = 1 + occupiedCarriers(k);
    end
end

pilotCarriers = eval(get(handles.pilotTonesPositions, 'string'));
for k = 1:length(pilotCarriers)
    if (pilotCarriers(k) < 0)
        pilotCarriers(k) = ifftLen + 1 + pilotCarriers(k);
    else
        pilotCarriers(k) = 1 + pilotCarriers(k);
    end
end
pilotSymbols = eval(get(handles.pilotSymbols, 'string'));

numPilots = length(pilotCarriers);
usePilots = 0;
if (numPilots > 0)
    usePilots = 1;
end

enablePackets = get(handles.enablePackets,'value');
if (enablePackets)
    packetHeadModTypeList = get(handles.packetHeadModType,'String');
    packetHeadModType = packetHeadModTypeList(get(handles.packetHeadModType, 'Value'),:);
    if (strcmp(packetHeadModType, 'BPSK'))
        packetHeadMod = comm.BPSKModulator('PhaseOffset', pi);
        packetHeadModTypeFlag = 0;
    elseif (strcmp(packetHeadModType, 'QPSK'))
        packetHeadMod = comm.QPSKModulator('PhaseOffset', 3/4*pi, 'BitInput',true);
        packetHeadModTypeFlag = 1;
    end
    sync1 = eval(get(handles.sync1, 'string'));
    sync2 = eval(get(handles.sync2, 'string'));
    packetLen = str2double(get(handles.packetLen,'string'));
else
    packetLen = nCarriers;
end
packetLenOriginal = packetLen;


%
%   Data bits length count.
%
if (enablePackets)
    dataLen = packetLen*nSymbols;   %in this case nSymbols is packet count
else
    dataLen = nCarriers*nSymbols*log2(M);
    sync1 = [];
    sync2 = [];
end

%
%   Source data
%   Here is just justification of type, generation is inside loop.
%       0 = randsrc
%       1 = user defined
%       2 = file
%
dataSourceFile = 0;
if (strcmp(handles.ofdm.dataInputMethod,'randsrc'))
    dataSourceType = 0;
elseif (strcmp(handles.ofdm.dataInputMethod,'user'))
    dataSourceType = 1;
    data_source = eval(get(handles.seqInput, 'string'))';
    if (enablePackets)
        data_source = cellstr(data_source');
        data_source = data_source{1};
        data_source = [data_source randi([32 126],1,packetLen-length(data_source)-1) 10]';
        data_source = double(data_source);
        data_source = repmat(data_source(1:packetLen),nSymbols,1);
    end
elseif (strcmp(handles.ofdm.dataInputMethod,'file'))
    dataSourceType = 2;
    dataSourceFileInfo = dir(handles.ofdm.file)
    dataSourceFile = fopen(handles.ofdm.file, 'r')
end

%
%   Vytvorenie modulatoru
%       handles.payloadMod
%       handles.packetHeadMod
%
clear modulated_data;
if (strcmp(payloadModType, 'QPSK'))
    payloadMod = comm.QPSKModulator('PhaseOffset', 3/4*pi, 'BitInput', true);
    payloadModType = 0;
elseif (strcmp(payloadModType, 'BPSK'))
    payloadMod = comm.BPSKModulator('PhaseOffset', pi);
    payloadModType = 1;
elseif (strcmp(payloadModType, 'PSK-8'))
    payloadMod = comm.PSKModulator(...
        8, ...
        'BitInput',true, ...
        'SymbolMapping','Custom', ...
        'CustomSymbolMapping', [0 4 5 1 3 7 6 2]);  %gnuradio psk8 constellation
    payloadModType = 2;
%     payloadMod.constellation;
end
handles.payloadMod = payloadMod;

%   Nastavenie parametrov podla GUI.
if (get(handles.noChannelRadio, 'Value') == 1)
    useChannelModel = 0;
elseif (get(handles.awgnRadio, 'Value') == 1)
    useChannelModel = 1;
    awgnLevel = str2double(get(handles.awgnLevel, 'string'));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Values checking!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isError = 0;    %clear error flag
if (sum(ismember(pilotCarriers, occupiedCarriers)))
    msgbox('There is pilot tone at subcarrier position!','Error!');
    isError = 1;
end
if ((max(occupiedCarriers) > ifftLen) || (min(occupiedCarriers) < 1))
    msgbox('Allocated carriers values should be inside <-ifftLen/2, ifftLen/2-1) interval.','Error!');
    isError = 1;
end
if (enablePackets && length(sync1) ~= ifftLen)
    msgbox('1st sync word should have same length as IFFT.','Error!');
    isError = 1;
end
if (enablePackets && length(sync2) ~= ifftLen)
    msgbox('2nd sync word should have same length as IFFT.','Error!');
    isError = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Jadro, OFDM modulator
%       - LW410 inicializacia
%       - USRP inicializacia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (deviceType == 0)
    nfor = 1;
    oversampleFactor = noDeviceOversampleFactor;
elseif (deviceType == 1)
    fc = str2double(get(handles.carrierFreqInput, 'string'))*1e3
    nfor = str2num(get(handles.editTxCycleCount,'String'));
    usrpIp = get(handles.USRPIPaddr, 'string');
    try
        Tx = comm.SDRuTransmitter(...
          'IPAddress', usrpIp, ...
          'CenterFrequency', fc, ...
          'InterpolationFactor', floor(100e6/fs), ...
          'Gain', txGain);
    catch ME
        msgbox('There is some problem with USRP radio!','Error!');
        isError = 1;
    end
elseif (deviceType == 2)
    nfor = 1;
    oversampleFactor = 1;
    fc = str2double(get(handles.lw410carrier, 'string'))*1e3
    handles.lw410.fs = lw410fs;
    
    try
        handles.lw410.idn();
    catch ME
        msgbox('There is some problem with LW410 generator!','Error!');
        isError = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   In ERROR case, clear graphs!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isError)
    cla(handles.ofdmWaveAxes);
    cla(handles.ofdmWaveAxes3);
    cla(handles.ofdmWaveAxes4);
    return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cyclicGeneration = 1;
loopCounter = 0;

while (cyclicGeneration == 1)
    loopCounter = loopCounter+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Jadro, OFDM modulator
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %
        %   debug
        df1 = fopen('debug/GnuradioTarget/ofdm_sig_outGUI.txt','w');
        %

        packetLen = packetLenOriginal;
        
        %
        %   DATA SOURCE, READING DATA
        %       0 = randi
        %       1 = user input
        %               - if packet enabled using as string message
        %       2 = file
        %
        if (dataSourceType == 0)
            if (enablePackets)
                data_source = randi(255,dataLen,1);
            else
                data_source = randi([0 1],dataLen,1);
            end
        elseif (dataSourceType == 1)
%             DO NOTHING, MESSAGE ALREADY SET BEFORE
%             if (enablePackets == 1)
%             else
%             end
        elseif (dataSourceType == 2)
            %reading samples from file, packets means read bytes, otherwise read raw bits
            if (enablePackets)
                data_source = fread(dataSourceFile, dataLen, 'uint8');
            else
                data_source = fread(dataSourceFile, dataLen, 'ubit1');
            end
            
            if (feof(dataSourceFile))            %check if start from beginning of file
                fseek(dataSourceFile,0,'bof');
            end
            
            %check if enough bits was read, if no, zeros are added in the end
            if (length(data_source) < dataLen)      
%                 data_source = [data_source; zeros(dataLen-length(data_source),1)];
                data_source = [data_source; randi([32 92],dataLen-length(data_source),1)];
            end            
        end
        
        if (enablePackets)
            modulated_data = [];
            o = 0;
            for k = 1:packetLen:dataLen%-packetLen+1
                dataIn = double(data_source(k:k+packetLen-1));

                header = generateHeader(packetLen+4,o);       %plus 4 because there is CRC added in payload    
                o = o + 1;
                
                %for 64
                %for 128
                %headerLen = 8*floor(nCarriers/8);             %not documented detail that header is increasing due to carriers
                
                %THIS WOULD BE LITTLE EMPIRICAL :)
                %for 512 is ceil()-1
                %for 1024 is ceil()-5  nCarriers = 275 pilots = 4
                %for 1024 is ceil()-2  nCarriers = 334 pilots = 4
                %for 1024 is ceil()-1  nCarriers = 351 pilots = 4
                %for 1024 is ceil()    nCarriers = 368 pilots = 4
                %headerLen = 8*ceil(nCarriers/8);
                
                %hypothesis LOOKS OK! (for now all pass)
                %for 1024 nCarriers = 418 pilots = 4
                headerLen = 8*floor(nCarriers/8)+rem(nCarriers,8);

                header = [header zeros(1,headerLen-length(header))];
                header = header(1:headerLen)';
                
                if (packetHeadModTypeFlag == 0)
                    modHeader = step(packetHeadMod, header);                    %<----- BPSK modulated header
                elseif(packetHeadModTypeFlag == 1)
                    modHeader = conj(step(packetHeadMod,[header;zeros(length(header),1)]));  %QPSK
                end
                
                if (payloadModType == 0)    %QPSK
                    payload = de2bi(generatePayload(dataIn'), 8)';
                    payload = reshape(payload,numel(payload),1);        %there is conversion from bytes --> bits inside! that's why reshape
                    modPayload = conj(step(payloadMod, payload));
                elseif(payloadModType == 1) %BPSK
                    payload = de2bi(generatePayload(dataIn'), 8)';
                    payload = reshape(payload,numel(payload),1);        %there is conversion from bytes --> bits inside! that's why reshape
                    modPayload = step(payloadMod, payload);
                elseif(payloadModType == 2) %PSK8
                    %for gnuradio groups of 3 bits are each reversed so we
                    %can't send just so to the modulator
                    %here is some monster reshaping which does the same job
                    a = de2bi(generatePayload(dataIn'),8)';
                    b = reshape(a,1,numel(a));
                    b = [b zeros(1,3-rem(numel(b),3))];
                    payload = reshape(fliplr(reshape(b,3,numel(b)/3)')',1,numel(b))';
                    modPayload = (step(payloadMod, payload));
                end

                modulated_data = [modulated_data; modHeader; modPayload];
            end
            packetLen = length(modulated_data)/nSymbols;           %recount length of packet (there is now payload with CRC and header with CRC)
        else
            modulated_data = step(payloadMod, data_source);
            packetLen = nCarriers;
        end

        %
        %   Carrier allocation.
        %        
        frames = allocateCarriers(         ...
            modulated_data',                ...
            ifftLen,                        ...
            packetLen,                      ...
            nSymbols,                       ...
            occupiedCarriers,               ...
            pilotCarriers,                  ...
            pilotSymbols,                   ...
            sync1,                          ...
            sync2)';
        packetLen = length(frames)/nSymbols;    %renew packetLen
        
%         frames = ifft(ifftshift(reshape(frames,packetLen,nSymbols)));
%         if (insertPrefixInterval == 1)
%             if (cpLen > 0)
%                 ofdm_signal = [guardLeft; frames(end-cpLen+1:end,:); frames; guardRight];
%             else
%                 ofdm_signal = [guardLeft; frames; guardRight];
%             end
%             ofdm_signal = reshape(ofdm_signal,numel(ofdm_signal),1);
%         end
%         ofdm_signal = ofdm_signal .* ifftLen .* signalAmplitude;

        ofdm_signal = [];
        for k = 1:ifftLen:packetLen*nSymbols-ifftLen+1
            ifftChunk = ifft(ifftshift(frames(k:k+ifftLen-1))) .* ifftLen .* signalAmplitude;
            %ifftChunk = ifft(ifftshift(frames(k:k+ifftLen-1))) .* signalAmplitude;

            %
            % Prefix/Guard interval insertion
            %
%             if (insertPrefixInterval == 1)
%                 if (cpLen > 0)
%                     ofdm_signal = [ofdm_signal; guardLeft; ifftChunk(end-cpLen+1:end); ifftChunk; guardRight];
%                 else
%                     ofdm_signal = [ofdm_signal; guardLeft; ifftChunk; guardRight];
%                 end
%             end
            if (insertPrefixInterval == 1)
                ofdm_signal = [ofdm_signal; ifftChunk(end-cpLen+1:end); ifftChunk];
            elseif (insertGuardInterval == 1)
                ofdm_signal = [ofdm_signal; guardLeft; ifftChunk; guardRight];
            else
                ofdm_signal = [ofdm_signal; ifftChunk];
            end
                
        end
                       
        %
        %   USRP Tx
        %

        if (useChannelModel == 1)
            ofdm_signal_noisy = awgn(ofdm_signal, awgnLevel, 20*log10(sum(abs(ofdm_signal))));
        else
            % useChannel == 0, default, no noise
            % no change of signal, ideal case
            ofdm_signal_noisy = ofdm_signal;
        end
        
        %
        %   Transmission
        %

        %Signal already conjugate somewhere on the way! NOT USED!
        %
        %ofdm_signal_noisy = 1i.*real(ofdm_signal)+imag(ofdm_signal);
        %psd(ofdm_signal_noisy);

        if (oversampleFactor > 1)
            ofdm_signal = resample(ofdm_signal',oversampleFactor,1)' *oversampleFactor;
            ofdm_signal_noisy = resample(ofdm_signal_noisy',oversampleFactor,1)' *oversampleFactor;
        end
%        
        if (deviceType == 1)            
            for cycleCnt = 1:1:nfor
                disp(horzcat('generate samples, step ',num2str(cycleCnt)));
                step(Tx, ofdm_signal_noisy);
            end
        elseif (deviceType == 2)
            
            %LW410 data output - we have to mixed baseband with carrier in
            %software
            
            ofdm_signal_noisy = resample(ofdm_signal_noisy,floor(lw410fs/fs),1);
            ds = 1/lw410fs;
            t = [0:ds:(length(ofdm_signal_noisy)-1)*ds]';
            disp(horzcat('Sequence length #samples: ', num2str(length(ofdm_signal_noisy))));
            set(handles.statusText, 'String', horzcat('LW410 uploading sequence. ', num2str(length(ofdm_signal_noisy)), ' samples'))
            ofdm_signal_noisy = ofdm_signal_noisy .* exp(1j*2*pi*fc*t);
            ofdm_signal_noisy = real(ofdm_signal_noisy) + imag(ofdm_signal_noisy);

            handles.lw410.wave_data(ofdm_signal_noisy', 1);
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Debuging part only
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
%         figure;
%         plot(real(frames));
%         hold on;
%         plot(imag(frames));
%         hold off;
        dSigReal = real(ofdm_signal);
        dSigImag = imag(ofdm_signal);
        fwrite(df1, reshape([dSigReal dSigImag]',1,2*length(dSigReal)), 'float32');
        fclose(df1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Interactive part
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pause(0.2);
    cyclicGeneration = get(handles.cyclicGeneration,'value');

    signalAmplitude = str2double(get(handles.editSignalAmplitude, 'string'));
    awgnLevel = str2double(get(handles.awgnLevel, 'string'));
    if (get(handles.noChannelRadio, 'Value') == 1)
        useChannelModel = 0;
    elseif (get(handles.awgnRadio, 'Value') == 1)
        useChannelModel = 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Showing OFDM frames
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %clip frames length to not overload graph with data
    if (length(frames) > 10e3) 
        frames = frames(1:10e3);
        frames = frames(1:10e3);
    end
    
    axes(handles.ofdmWaveAxes3);
    stem(real(frames),'marker','none');
    set(gca,'Xtick',1:ifftLen/2:length(frames)+1,'XTickLabel',{['+' num2str(ifftLen) '|-' num2str(ifftLen)] 0},'fontsize',8);
    grid on;
    title('OFDM Frames in frequency spectrum')
    ylabel('U[V]');
    hold on;
    stem(imag(frames),'-r','marker','none')
    hold off;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Signal in time plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %clip signal to show, to not to display too much, calculations next are
    %based on this signal!
    if (length(ofdm_signal) > 10e3) 
        ofdm_signal = ofdm_signal(1:10e3);
        ofdm_signal_noisy = ofdm_signal_noisy(1:10e3);
    end
    
    axes(handles.ofdmWaveAxes);
    plot(zeros(1,ifftLen));
    title('OFDM signal');
    xlabel('sample[-]');
    ylabel('A[V]');
%     plot(linspace(0,fs,length(ofdm_signal)), real(ofdm_signal)+imag(ofdm_signal));
%     plot(linspace(0,fs,length(ofdm_signal)), real(ofdm_signal));
%     hold on; plot(linspace(0,fs,length(ofdm_signal)), imag(ofdm_signal),'r'); hold off;
    
    timeAxis_ms = [1:length(ofdm_signal)]*1/fs*1000;
    plot(timeAxis_ms, real(ofdm_signal)+imag(ofdm_signal));
    xlim([timeAxis_ms(1) timeAxis_ms(end)]);
    title('OFDM signal in time domain');
%     hold on; plot(timeAxis_ms, imag(ofdm_signal),'r'); hold off;
    grid on;
    xlabel('t[ms]'); ylabel('A[V]');

    disp('Energia OFDM signalu (dBm)');
    10*log10(sum(abs(ofdm_signal).^2))

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Signal after channel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sigSpec = 20*log10(abs(fft(ofdm_signal))./length(ofdm_signal));
    dfs = fs/length(ofdm_signal);
    minX = -fs/2;
    maxX = abs(minX)-dfs;
    xAxis = linspace(minX, maxX, length(ofdm_signal)).*1e-3;

    axes(handles.ofdmWaveAxes4);
    plot(xAxis, fftshift(sigSpec));

    ylim([floor(min(sigSpec)/10)*10 ceil(max(sigSpec)/10)*10]);
    xlim([xAxis(1) xAxis(end)]);

    if (useChannelModel == 0)
        title('Bandwidth signal spectrum'); xlabel('f[kHz]'); ylabel('A[dBm]');
    elseif (useChannelModel == 1)
        title('Signal after AWGN channel'); xlabel('f[kHz]'); ylabel('A[dBm]');
    else
        axes(handles.ofdmWaveAxes4);
        plot([]);
    end
    grid on;    

    % NEZABUDNUT ZE VYPOCET ENERGIE SIGNALU A ZASUMENEHO SIGNALU SA
    % ROBI Z PRESEMPLOVANEHO SIGNALU! <-- overit mozno je to dobre

    ofdm_signal_energy_noisy = 10*log10(sum(abs(ofdm_signal_noisy).^2));
    ofdm_signal_energy = 10*log10(sum(abs(ofdm_signal).^2));
    if (useChannelModel == 0)       % idealny kanal, vypis sily signalu
        set(handles.statusText, 'String', ...
            horzcat(    ...
                'Generating samples. no channel. ', ...
                'Signal energy ', ...
                num2str(ofdm_signal_energy),'dBm.', ...
                'Total channel energy ', ...
                num2str(ofdm_signal_energy_noisy),'dBm.' ...
            ));
    elseif (useChannelModel == 1)   % AWGN kanal parametre vypis
        set(handles.statusText, 'String', ...
            horzcat(    ...
                'Generating samples. Using AWGN channel. ', ...
                'Signal energy ', ...
                num2str(ofdm_signal_energy),'dBm.', ...
                'Total channel energy ', ...
                num2str(ofdm_signal_energy_noisy),'dBm.' ...
            ));
    else
    end
end     %koniec modulovania

if (dataSourceType == 2)
    fclose(dataSourceFile);
end

if (deviceType == 1)
    release(Tx);
end
if (enablePackets == 1)
    release(packetHeadMod);
end
release(payloadMod);

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in gpibAddrPopup.
function cyclicGeneration_Callback(hObject, eventdata, handles)

% --- Executes on selection change in gpibAddrPopup.
function gpibAddrPopup_Callback(hObject, eventdata, handles)
% hObject    handle to gpibAddrPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns gpibAddrPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from gpibAddrPopup


% --- Executes during object creation, after setting all properties.
function gpibAddrPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gpibAddrPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in detectDeviceBut.
function detectDeviceBut_Callback(hObject, eventdata, handles)
% hObject    handle to detectDeviceBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%TEST if is not empty device
% if (empty(handles.lw410))
%     handles.lw410 = LW410Interface();
% end
    contents = cellstr(get(handles.gpibAddrPopup,'String'));
    gpib_addr = str2num(contents{get(handles.gpibAddrPopup,'Value')});

    contents = cellstr(get(handles.sampleFreqPopup,'String'));
    sampleFreq = contents{get(handles.sampleFreqPopup,'Value')};
    if (strcmp(sampleFreq, '400 MHz'))
        sampleFreq = 400e6;
    elseif (strcmp(sampleFreq, '40 MHz'))
        sampleFreq = 40e6;
    elseif (strcmp(sampleFreq, '4 MHz'))
        sampleFreq = 4e6;
    elseif (strcmp(sampleFreq, '40 kHz'))
        sampleFreq = 40e3;
    else
        sampleFreq = 4e6;           % DEFAULT VALUE
    end
    
    set(handles.statusText, 'String', horzcat('Looking for LW410 at GPIB addr ', num2str(gpib_addr)));
    pause(0.3);

    handles.lw410.gpib_addr = gpib_addr;
    handles.lw410.fs = sampleFreq;
    id = handles.lw410.idn();
    set(handles.statusText, 'String', id);

guidata(hObject, handles);


% --- Executes on button press in numPeriodsEdit.
function numPeriodsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to numPeriodsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function numPeriodsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numPeriodsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in modulationTypePopup.
function modulationTypePopup_Callback(hObject, eventdata, handles)
% hObject    handle to modulationTypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns modulationTypePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from modulationTypePopup
changeHintsText(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function modulationTypePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modulationTypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in numOfCarriersPopup.
function numOfCarriersPopup_Callback(hObject, eventdata, handles)
% hObject    handle to numOfCarriersPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns numOfCarriersPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from numOfCarriersPopup


% --- Executes during object creation, after setting all properties.
function numOfCarriersPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numOfCarriersPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dataInputFileRadio.
function dataInputFileRadio_Callback(hObject, eventdata, handles)
% hObject    handle to dataInputFileRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dataInputFileRadio

%   it's triggered two times first this after selection changed, needs resolve!
% dataInputGroup_SelectionChangeFcn(hObject, eventdata, handles);
set(handles.seqInput,'Enable','off');
hintsText = {...
    'Input data source: file' ...
    'Input Data:' ...
    '    data would be read as bits from file' ...
};
set(handles.hintsText,'String',hintsText);

guidata(hObject, handles);


% --- Executes when selected object is changed in dataInputGroup.
function dataInputGroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in dataInputGroup 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if (get(handles.dataInputRandomRadio, 'Value') == 1) handles.ofdm.dataInputMethod =  'randsrc';
elseif (get(handles.dataInputUserRadio, 'Value') == 1) handles.ofdm.dataInputMethod = 'user';
elseif (get(handles.dataInputFileRadio, 'Value') == 1)
    handles.ofdm.dataInputMethod = 'file';
    [fileName,path] = uigetfile('*','Select file with samples');
    if (fileName ~= 0)
        handles.ofdm.file = horzcat(path, fileName);
        set(handles.seqFileInput, 'string', horzcat(path, fileName));
        set(handles.statusText, 'string', horzcat('Data input from file ', path, fileName));
    else
        set(handles.dataInputRandomRadio, 'Value', 1);
        handles.ofdm.dataInputMethod =  'randsrc';
        msgbox('File hasn''t been choosed, input changed to random sequence.','Error!');        
    end
end

disp(horzcat('Data input method changed to ', handles.ofdm.dataInputMethod));
guidata(hObject, handles);

% --- Executes on button press in butExportOFDMData.
function butExportOFDMData_Callback(hObject, eventdata, handles)
% hObject    handle to butExportOFDMData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hChildren = get(handles.ofdmWaveAxes,'Children');
sigReal = get(hChildren(1),'YData');

[filePath dir] = uiputfile('*.*', 'Choose output file');
if (exist('filepath'))
    f = fopen(horzcat(dir,filePath),'w');
    fwrite(f, sigReal, 'float32');
    fclose(f);
end

guidata(hObject, handles);


% --- Executes on selection change in sampleFreqPopup.
function sampleFreqPopup_Callback(hObject, eventdata, handles)
% hObject    handle to sampleFreqPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sampleFreqPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sampleFreqPopup


% --- Executes during object creation, after setting all properties.
function sampleFreqPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sampleFreqPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ifftLen_Callback(hObject, eventdata, handles)
% hObject    handle to ifftLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ifftLen as text
%        str2double(get(hObject,'String')) returns contents of ifftLen as a double
changeHintsText(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function ifftLen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ifftLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function prefixLengthInput_Callback(hObject, eventdata, handles)
% hObject    handle to prefixLengthInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prefixLengthInput as text
%        str2double(get(hObject,'String')) returns contents of prefixLengthInput as a double
ifftLen = get(handles.ifftLen,'String');
cpLen = get(handles.prefixLengthInput,'String');
try
    cpLen = eval(cpLen);
catch ME
    msgbox('Cyclic prefix is not number! Changed to 8','Error!');
    set(handles.prefixLengthInput,'String','8');
end
if (cpLen > ifftLen)
    msgbox(['Cyclic prefix must be smaller than ifft length! Changed to ' ifftLen],'Error!');
    set(handles.prefixLengthInput,'String',ifftLen);
end


% --- Executes during object creation, after setting all properties.
function prefixLengthInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prefixLengthInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in intervalRadio.
function intervalRadio_Callback(hObject, eventdata, handles)
% hObject    handle to intervalRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of intervalRadio
set(handles.prefixLengthInput,'Enable','on');
set(handles.guardLengthInput,'Enable','off');
guidata(hObject, handles);

% --- Executes on button press in intervalGuardIntervalRadio.
function intervalGuardIntervalRadio_Callback(hObject, eventdata, handles)
% hObject    handle to intervalGuardIntervalRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of intervalGuardIntervalRadio



function oversamplingInput_Callback(hObject, eventdata, handles)
% hObject    handle to oversamplingInputText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of oversamplingInputText as text
%        str2double(get(hObject,'String')) returns contents of oversamplingInputText as a double


% --- Executes during object creation, after setting all properties.
function oversamplingInputText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to oversamplingInputText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function dataInputGroup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataInputGroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function numOfSymbolsInput_Callback(hObject, eventdata, handles)
% hObject    handle to numOfSymbolsInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numOfSymbolsInput as text
%        str2double(get(hObject,'String')) returns contents of numOfSymbolsInput as a double
changeHintsText(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function numOfSymbolsInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numOfSymbolsInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function carrierFreqInput_Callback(hObject, eventdata, handles)
% hObject    handle to carrierFreqInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of carrierFreqInput as text
%        str2double(get(hObject,'String')) returns contents of carrierFreqInput as a double


% --- Executes during object creation, after setting all properties.
function carrierFreqInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to carrierFreqInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to oversamplingInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of oversamplingInput as text
%        str2double(get(hObject,'String')) returns contents of oversamplingInput as a double


% --- Executes during object creation, after setting all properties.
function oversamplingInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to oversamplingInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function seqInput_Callback(hObject, eventdata, handles)
% hObject    handle to seqInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seqInput as text
%        str2double(get(hObject,'String')) returns contents of seqInput as a double
enablePackets = get(handles.enablePackets, 'Value');
seqInput = get(handles.seqInput,'String');
if (enablePackets)
    try
        seqInput = eval(seqInput);
    catch ME
        msgbox('Input sequence is not string, changed to ''Hello world!''','Error!');
        set(handles.seqInput,'String','''Hello world!''');
    end
else
    try
        seqInput = eval(seqInput);
    catch ME
        msgbox('Input sequence is not array, changed to default []','Error!');
        set(handles.seqInput,'String','[]');
    end
end

% --- Executes during object creation, after setting all properties.
function seqInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seqInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in showOfdmWaveAxesBut.
function showOfdmWaveAxesBut_Callback(hObject, eventdata, handles)
% hObject    handle to showOfdmWaveAxesBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
scrsz = get(0,'ScreenSize');
h = figure('Tag','time_zoom','Name','OFDM signal','Position',[scrsz(3)/2-500 scrsz(4)/2-250 1000 500]);
newaxes = copyobj(handles.ofdmWaveAxes, h);
set(newaxes,'Units','normalized','ActivePositionProperty','position','Position',[0.1 0.1 0.85 0.85]);

% --- Executes on button press in showOfdmWaveAxes4But.
function showOfdmWaveAxes4But_Callback(hObject, eventdata, handles)
% hObject    handle to showOfdmWaveAxes4But (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
scrsz = get(0,'ScreenSize');
h = figure('Tag','time_zoom','Name','OFDM phase spectrum','Position',[scrsz(3)/2-500 scrsz(4)/2-250 1000 500]);
newaxes = copyobj(handles.ofdmWaveAxes4, h);
set(newaxes,'Units','normalized','ActivePositionProperty','position','Position',[0.1 0.1 0.85 0.85]);

function awgnLevel_Callback(hObject, eventdata, handles)
% hObject    handle to awgnLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of awgnLevel as text
%        str2double(get(hObject,'String')) returns contents of awgnLevel as a double

% --- Executes during object creation, after setting all properties.
function awgnLevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to awgnLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editSignalAmplitude_Callback(hObject, eventdata, handles)
% hObject    handle to editSignalAmplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSignalAmplitude as text
%        str2double(get(hObject,'String')) returns contents of editSignalAmplitude as a double


% --- Executes during object creation, after setting all properties.
function editSignalAmplitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSignalAmplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editTxCycleCount_Callback(hObject, eventdata, handles)
% hObject    handle to editTxCycleCount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTxCycleCount as text
%        str2double(get(hObject,'String')) returns contents of editTxCycleCount as a double


% --- Executes during object creation, after setting all properties.
function editTxCycleCount_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTxCycleCount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editUsrpInterp_Callback(hObject, eventdata, handles)
% hObject    handle to editUsrpInterp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editUsrpInterp as text
%        str2double(get(hObject,'String')) returns contents of editUsrpInterp as a double


% --- Executes during object creation, after setting all properties.
function editUsrpInterp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editUsrpInterp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function headDataIn_Callback(hObject, eventdata, handles)
% hObject    handle to headDataIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of headDataIn as text
%        str2double(get(hObject,'String')) returns contents of headDataIn as a double


% --- Executes during object creation, after setting all properties.
function headDataIn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to headDataIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pilotTonesPositions_Callback(hObject, eventdata, handles)
% hObject    handle to pilotTonesPositions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pilotTonesPositions as text
%        str2double(get(hObject,'String')) returns contents of pilotTonesPositions as a double
carriers = eval(get(handles.occupiedCarriers, 'String'));
pilotTonesPositions = get(handles.pilotTonesPositions, 'String');
try
    pilotTonesPositions = eval(pilotTonesPositions);
catch ME
    msgbox('Subcarriers mapping doesn''t contain array!','Error!');
end
if (sum(ismember(pilotTonesPositions, carriers)))
    msgbox('There is pilot tone at subcarrier position!','Error!');
end


% --- Executes during object creation, after setting all properties.
function pilotTonesPositions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pilotTonesPositions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pilotSymbols_Callback(hObject, eventdata, handles)
% hObject    handle to pilotSymbols (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pilotSymbols as text
%        str2double(get(hObject,'String')) returns contents of pilotSymbols as a double
pilotTonesAmplitudes = get(handles.pilotSymbols, 'String');
try
    pilotTonesAmplitudes = eval(pilotTonesAmplitudes);
catch ME
    msgbox('Pilot tones amplitudes is not array!','Error!');
end

% --- Executes during object creation, after setting all properties.
function pilotSymbols_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pilotSymbols (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveCurrentSettings.
function saveCurrentSettings_Callback(hObject, eventdata, handles)
% hObject    handle to saveCurrentSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
%   ToDo:
%       - generating file names by date
%

modTypeList = get(handles.modulationTypePopup,'String');
payloadModType = modTypeList{get(handles.modulationTypePopup, 'Value')}

contents = cellstr(get(handles.gpibAddrPopup,'String'));
gpib_addr = str2num(contents{get(handles.gpibAddrPopup,'Value')});
gpib_addr = cellstr(num2str(gpib_addr));

contents = cellstr(get(handles.sampleFreqPopup,'String'));
sampleFreq = contents{get(handles.sampleFreqPopup,'Value')};

paramsToSave = [ ...
    cellstr(num2str(get(handles.dataInputRandomRadio,'Value')))   ...
    cellstr(num2str(get(handles.dataInputUserRadio, 'Value')))    ...
    cellstr(num2str(get(handles.dataInputFileRadio, 'Value')))    ...
    cellstr(get(handles.seqInput, 'String')) ...
    cellstr(get(handles.seqFileInput, 'String')) ...
    cellstr(payloadModType) ...
    cellstr(get(handles.ifftLen, 'String')) ...
    cellstr(get(handles.editSignalAmplitude, 'String')) ...    
    cellstr(get(handles.numOfSymbolsInput, 'String')) ...
    cellstr(num2str(get(handles.intervalRadio,'Value')))  ...
    cellstr(num2str(get(handles.intervalNotusedRadio,'Value')))   ...
    cellstr(get(handles.prefixLengthInput, 'String')) ...
    cellstr(get(handles.guardLengthInput, 'String')) ...
    cellstr(get(handles.occupiedCarriers, 'String')) ...
    cellstr(get(handles.pilotTonesPositions, 'String')) ...
    cellstr(get(handles.pilotSymbols, 'String')) ...
    cellstr(num2str(get(handles.enablePackets,'Value')))   ...
    cellstr(get(handles.packetLen,'String'))   ...
    cellstr(get(handles.sync1,'String'))   ...
    cellstr(get(handles.sync2,'String'))   ...
    cellstr(num2str(get(handles.noChannelRadio,'Value')))   ...
    cellstr(num2str(get(handles.awgnRadio,'Value')))   ...
    cellstr(get(handles.awgnLevel,'String'))   ...
    cellstr(num2str(get(handles.noDevice,'Value')))   ...
    cellstr(num2str(get(handles.usrpDevice,'Value')))   ...
    cellstr(num2str(get(handles.lw410Device,'Value')))   ...
    cellstr(get(handles.noDeviceSampleRate,'String'))   ...
    cellstr(get(handles.USRPSampleRate,'String'))   ...
    cellstr(get(handles.oversamplingInput,'String'))   ...
    cellstr(get(handles.editTxCycleCount,'String'))   ...
    cellstr(get(handles.carrierFreqInput,'String'))   ...
    cellstr(get(handles.txGain,'String'))   ...
    cellstr(get(handles.USRPIPaddr,'String'))   ...
    cellstr(get(handles.lw410SampleRate,'String'))   ...
    cellstr(get(handles.lw410carrier,'String'))   ...
    cellstr(gpib_addr)   ...
    cellstr(sampleFreq)  ...
    cellstr(num2str(get(handles.cyclicGeneration,'Value')))   ...    
];

[fileName,path] = uiputfile('*','Select file to save settings');
f = fopen(horzcat(path,fileName),'w');
fprintf(f,repmat('%s\n',1,length(paramsToSave)),paramsToSave{1:end});
fclose(f);

% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
paramsToRead = [];

[fileName,path] = uigetfile('*','Select file with head data');

f = fopen(horzcat(path,fileName),'r');
for k = 1:37
    paramsToRead = [paramsToRead; cellstr(fgets(f))];
end
fclose(f);

set(handles.dataInputRandomRadio,'Value', str2num(paramsToRead{1}))
set(handles.dataInputUserRadio, 'Value', str2num(paramsToRead{2}))
set(handles.dataInputFileRadio, 'Value', str2num(paramsToRead{3}))
set(handles.seqInput, 'String', paramsToRead{4})
set(handles.seqFileInput, 'String', paramsToRead{5})
%     str2num(payloadModType) 
set(handles.ifftLen, 'String', paramsToRead{7})
set(handles.editSignalAmplitude, 'String', paramsToRead{8})    
set(handles.numOfSymbolsInput, 'String', paramsToRead{9})
set(handles.intervalRadio,'Value', str2num(paramsToRead{10}))  
set(handles.intervalNotusedRadio,'Value', str2num(paramsToRead{11})) 
set(handles.prefixLengthInput, 'String', paramsToRead{12})
set(handles.guardLengthInput, 'String', paramsToRead{13})
set(handles.occupiedCarriers, 'String', paramsToRead{14})
set(handles.pilotTonesPositions, 'String', paramsToRead{15})
set(handles.pilotSymbols, 'String', paramsToRead{16})
set(handles.enablePackets,'Value', str2num(paramsToRead{17}))
set(handles.packetLen,'String', paramsToRead{18}) 
set(handles.sync1,'String', paramsToRead{19})
set(handles.sync2,'String', paramsToRead{20})
set(handles.noChannelRadio,'Value', str2num(paramsToRead{21}))
set(handles.awgnRadio,'Value', str2num(paramsToRead{22}))
set(handles.awgnLevel,'String', paramsToRead{23})
set(handles.noDevice,'Value', str2num(paramsToRead{24}))
set(handles.usrpDevice,'Value', str2num(paramsToRead{25}))
set(handles.lw410Device,'Value', str2num(paramsToRead{26}))
set(handles.noDeviceSampleRate,'String', paramsToRead{27})
set(handles.USRPSampleRate,'String', paramsToRead{28})
set(handles.oversamplingInput,'String', paramsToRead{29})
set(handles.editTxCycleCount,'String', paramsToRead{30})
set(handles.carrierFreqInput,'String', paramsToRead{31}) 
set(handles.txGain,'String', paramsToRead{32})
set(handles.USRPIPaddr,'String', paramsToRead{33}) 
set(handles.lw410SampleRate,'String', paramsToRead{34}) 
set(handles.lw410SampleRate,'String', paramsToRead{35}) 
%     gpib_addr  
%     sampleFreq 
% set(handles.cyclicGeneration,'Value', str2num(paramsToRead{38}))

% change deviceOutputImg
connectedDeviceGroup_SelectionChangedFcn(hObject, eventdata, handles)

guidata(hObject, handles);


function subcarriersPos_Callback(hObject, eventdata, handles)
% hObject    handle to subcarriersPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subcarriersPos as text
%        str2double(get(hObject,'String')) returns contents of subcarriersPos as a double


% --- Executes during object creation, after setting all properties.
function subcarriersPos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subcarriersPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function subcarriersAmpl_Callback(hObject, eventdata, handles)
% hObject    handle to subcarriersAmpl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subcarriersAmpl as text
%        str2double(get(hObject,'String')) returns contents of subcarriersAmpl as a double


% --- Executes during object creation, after setting all properties.
function subcarriersAmpl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subcarriersAmpl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to editUsrpInterp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editUsrpInterp as text
%        str2double(get(hObject,'String')) returns contents of editUsrpInterp as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editUsrpInterp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to editTxCycleCount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTxCycleCount as text
%        str2double(get(hObject,'String')) returns contents of editTxCycleCount as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTxCycleCount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to carrierFreqInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of carrierFreqInput as text
%        str2double(get(hObject,'String')) returns contents of carrierFreqInput as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to carrierFreqInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to oversamplingInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of oversamplingInput as text
%        str2double(get(hObject,'String')) returns contents of oversamplingInput as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to oversamplingInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function doplerShiftIn_Callback(hObject, eventdata, handles)
% hObject    handle to doplerShiftIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of doplerShiftIn as text
%        str2double(get(hObject,'String')) returns contents of doplerShiftIn as a double


% --- Executes during object creation, after setting all properties.
function doplerShiftIn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to doplerShiftIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in usrpContinuous.
function usrpContinuous_Callback(hObject, eventdata, handles)
% hObject    handle to usrpContinuous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of usrpContinuous


% --- Executes on key press with focus on dataInputFileRadio and none of its controls.
function dataInputFileRadio_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to dataInputFileRadio (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function dataInputGroup_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to dataInputGroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes when selected object is changed in headData.
function headData_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in headData 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (get(handles.noHeadDataRadio, 'Value') == 1)
    handles.ofdm.headData = 0;
elseif (get(handles.headDataRadio, 'Value') == 1)
    handles.ofdm.headData = 1;
    [fileName,path] = uigetfile('*','Select file with head data');
    handles.ofdm.headFile = horzcat(path, fileName);
    set(handles.statusText, 'String', horzcat('Header data read from file ', handles.ofdm.headFile));
end
guidata(hObject, handles);


% --- Executes on button press in writeHeadAndDataBut.
function writeHeadAndDataBut_Callback(hObject, eventdata, handles)
% hObject    handle to writeHeadAndDataBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('WriteBinaryFile start.');
WriteBinaryFile(handles.modM,64,20e3);
set(handles.statusText, 'String', 'Head and data write as binary files OFDMHead.bin and OFDMData.bin');
disp('WriteBinary file end.');



function guardLengthInput_Callback(hObject, eventdata, handles)
% hObject    handle to guardLengthInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of guardLengthInput as text
%        str2double(get(hObject,'String')) returns contents of guardLengthInput as a double


% --- Executes during object creation, after setting all properties.
function guardLengthInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to guardLengthInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in signalOutputBut.
function signalOutputBut_Callback(hObject, eventdata, handles)
% hObject    handle to signalOutputBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data.activePanel = 'signalOutputPanel';   %set semaphor
set(handles.deviceSettingsPanel, 'visible','off');
set(handles.ofdmWavePanel,'visible','off');
set(handles.signalOutputPanel,'visible','on');

% --- Executes on button press in showOfdmWaveAxes3But.
function showOfdmWaveAxes3But_Callback(hObject, eventdata, handles)
% hObject    handle to showOfdmWaveAxes3But (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
scrsz = get(0,'ScreenSize');
h = figure('Tag','time_zoom','Name','OFDM packets spectrum','Position',[scrsz(3)/2-500 scrsz(4)/2-250 1000 500]);
newaxes = copyobj(handles.ofdmWaveAxes3, h);
set(newaxes,'Units','normalized','ActivePositionProperty','position','Position',[0.1 0.1 0.85 0.85]);


% --- Executes on button press in dataInputUserRadio.
function dataInputUserRadio_Callback(hObject, eventdata, handles)
% hObject    handle to dataInputUserRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dataInputUserRadio
set(handles.seqInput,'Enable','on');

occupiedCarriers = eval(get(handles.occupiedCarriers, 'string'));
nCarriers = length(occupiedCarriers);
nSymbols = str2num(get(handles.numOfSymbolsInput, 'string'));

enablePackets = get(handles.enablePackets, 'Value');
if (enablePackets)
    packetLen = str2double(get(handles.packetLen,'string'));
    hintsText = {...
        'Input data source: user-defined, packet transfer' ...
        'Expected Input Data:' ...
        '    string closed by quotes, ie. ''Hello world'' ' ...
    };
else
    modTypeList = get(handles.modulationTypePopup,'String');
    payloadModType = modTypeList{get(handles.modulationTypePopup, 'Value')}
    if (strcmp(payloadModType, 'QPSK'))
        M = 4
    elseif (strcmp(payloadModType, 'BPSK'))
        M = 2
    elseif (strcmp(payloadModType, 'PSK-8'))
        M = 8
    end

    hintsText = {...
        'Input data source: user-defined, ' ...
        ['Data modulation type: ' payloadModType] ...
        ['Number of Carriers: ' num2str(nCarriers)]...
        ['Expected Input Data Length (bits 0,1): ' num2str(nCarriers) '*log2(' num2str(M) ')*' num2str(nSymbols) ' = ' num2str(nCarriers*log2(M)*nSymbols)] ...
        '' ...
        'Input data are generated automatically as random sequence:' ...
        '    randi([0 1],1,dataLen)' ...
        '' ...
        'Available variables:' ...
        '    dataLen = count of bits needed to generate' ...
    };
end
set(handles.hintsText,'String',hintsText);

guidata(hObject, handles);


function occupiedCarriers_Callback(hObject, eventdata, handles)
% hObject    handle to occupiedCarriers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of occupiedCarriers as text
%        str2double(get(hObject,'String')) returns contents of occupiedCarriers as a double
ifftLen = str2num(get(handles.ifftLen, 'String'));
carriers = get(handles.occupiedCarriers, 'String');
pilotTonesPositions = eval(get(handles.pilotTonesPositions, 'String'));
try
    carriers = eval(carriers);
catch ME
    msgbox('Subcarriers mapping doesn''t contain array!','Error!');
end
if (sum(ismember(pilotTonesPositions, carriers)))
    msgbox('There is pilot tone at subcarrier position!','Error!');
end
if (~isempty(carriers))
    if ((max(carriers) > ifftLen/2-1) || (min(carriers) < -ifftLen/2))
        msgbox('Allocated carriers values should be inside <-ifftLen/2, ifftLen/2-1) interval.','Error!');
    end
end

changeHintsText(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function occupiedCarriers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to occupiedCarriers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sync1_Callback(hObject, eventdata, handles)
% hObject    handle to sync1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sync1 as text
%        str2double(get(hObject,'String')) returns contents of sync1 as a double


% --- Executes during object creation, after setting all properties.
function sync1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sync1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sync2_Callback(hObject, eventdata, handles)
% hObject    handle to sync2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sync2 as text
%        str2double(get(hObject,'String')) returns contents of sync2 as a double


% --- Executes during object creation, after setting all properties.
function sync2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sync2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in enablePackets.
function enablePackets_Callback(hObject, eventdata, handles)
% hObject    handle to enablePackets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of enablePackets
if (get(hObject,'Value') == 1)
    set(handles.symbolLengthLabel,'String', 'Number of packets');
    set(handles.dataInputUserRadio,'String', 'user string');
    set(handles.seqInput,'String', '''''');
    
    set(handles.packetHeadModType,'Enable', 'on');
    set(handles.packetLen,'Enable', 'on');
    set(handles.sync1,'Enable', 'on');
    set(handles.sync2,'Enable', 'on');
else
    set(handles.symbolLengthLabel,'String', 'Number of symbols');
    set(handles.dataInputUserRadio,'String', 'user data');
    set(handles.seqInput,'String', '[]');

    set(handles.packetHeadModType,'Enable', 'off');
    set(handles.packetLen,'Enable', 'off');
    set(handles.sync1,'Enable', 'off');
    set(handles.sync2,'Enable', 'off');
end
changeHintsText(hObject, eventdata, handles);

function packetLen_Callback(hObject, eventdata, handles)
% hObject    handle to packetLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of packetLen as text
%        str2double(get(hObject,'String')) returns contents of packetLen as a double
changeHintsText(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function packetLen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to packetLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function USRPSampleRate_Callback(hObject, eventdata, handles)
% hObject    handle to USRPSampleRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of USRPSampleRate as text
%        str2double(get(hObject,'String')) returns contents of USRPSampleRate as a double


% --- Executes during object creation, after setting all properties.
function USRPSampleRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to USRPSampleRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function noDeviceSampleRate_Callback(hObject, eventdata, handles)
% hObject    handle to noDeviceSampleRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noDeviceSampleRate as text
%        str2double(get(hObject,'String')) returns contents of noDeviceSampleRate as a double


% --- Executes during object creation, after setting all properties.
function noDeviceSampleRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noDeviceSampleRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in packetHeadModType.
function packetHeadModType_Callback(hObject, eventdata, handles)
% hObject    handle to packetHeadModType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns packetHeadModType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from packetHeadModType


% --- Executes during object creation, after setting all properties.
function packetHeadModType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to packetHeadModType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txGain_Callback(hObject, eventdata, handles)
% hObject    handle to txGain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txGain as text
%        str2double(get(hObject,'String')) returns contents of txGain as a double


% --- Executes during object creation, after setting all properties.
function txGain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txGain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in intervalNotusedRadio.
function intervalNotusedRadio_Callback(hObject, eventdata, handles)
% hObject    handle to intervalNotusedRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of intervalNotusedRadio
set(handles.prefixLengthInput,'Enable','off');
set(handles.guardLengthInput,'Enable','off');
guidata(hObject, handles);

function USRPIPaddr_Callback(hObject, eventdata, handles)
% hObject    handle to USRPIPaddr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of USRPIPaddr as text
%        str2double(get(hObject,'String')) returns contents of USRPIPaddr as a double


% --- Executes during object creation, after setting all properties.
function USRPIPaddr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to USRPIPaddr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function seqFileInput_Callback(hObject, eventdata, handles)
% hObject    handle to seqFileInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seqFileInput as text
%        str2double(get(hObject,'String')) returns contents of seqFileInput as a double


% --- Executes during object creation, after setting all properties.
function seqFileInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seqFileInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit48_Callback(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit48 as text
%        str2double(get(hObject,'String')) returns contents of edit48 as a double


% --- Executes during object creation, after setting all properties.
function edit48_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fs = str2double(get(handles.noDeviceSampleRate, 'string'))*1000;
fc = 7000000;
if (get(handles.usrpDevice, 'Value') == 1)
    oversampling = str2double(get(handles.oversamplingInput, 'string'));
    fs = str2double(get(handles.USRPSampleRate, 'string'))*1000/oversampling;
    fc = str2double(get(handles.carrierFreqInput, 'string'))*1000;
    deviceType = 1;
end;
if (get(handles.lw410Device, 'Value') == 1)
    fs = str2double(get(handles.lw410SampleRate, 'string'))*1000;
    fc = str2double(get(handles.lw410carrier, 'string'))*1000;
    deviceType = 2;
end;

sync1 = eval(get(handles.sync1, 'string'));
sync2 = eval(get(handles.sync2, 'string'));
ifftLen = str2num(get(handles.ifftLen, 'string'));
cpLen = str2double(get(handles.prefixLengthInput, 'string'));
packetLen = str2num(get(handles.packetLen, 'string'));
pilotSymbols = eval(get(handles.pilotSymbols, 'string'));
pilotCarriers = eval(get(handles.pilotTonesPositions, 'string'));
occupiedCarriers = eval(get(handles.occupiedCarriers, 'string'));

packetHeadModTypeList = get(handles.packetHeadModType,'String');
modTypeList = get(handles.modulationTypePopup,'String');
payloadModType = modTypeList{get(handles.modulationTypePopup, 'Value')};
packetHeadModType = packetHeadModTypeList(get(handles.packetHeadModType, 'Value'),:);

txTemplate = regexp( fileread('templates/OFDM_Receiver.grc'),'\n','split');

txTemplate{567} = ['<value>' num2str(fs) '</value>'];
txTemplate{108} = ['<value>' num2str(fc) '</value>'];
txTemplate{216} = ['<value>' num2str(cpLen) '</value>'];

tsync1 = '<value>(';
for k = 1:length(sync1)
    if (k ~= 1)
         tsync1 = [tsync1 ',' num2str(sync1(k))]; 
    else
         tsync1 = [tsync1 num2str(sync1(k))]; 
    end
end
tsync1 = [tsync1 ')</value>'];
txTemplate{594} = tsync1;

tsync2 = '<value>(';
for k = 1:length(sync2)
    if (k ~= 1)
         tsync2 = [tsync2 ',' num2str(sync2(k))]; 
    else
         tsync2 = [tsync2 num2str(sync2(k))]; 
    end
end
tsync2 = [tsync2 ')</value>'];
txTemplate{621} = tsync2;

txTemplate{243} = ['<value>' num2str(ifftLen) '</value>'];
txTemplate{405} = ['<value>' num2str(packetLen) '</value>'];

tpilotCarriers = '<value>((';
for k = 1:length(pilotCarriers)
    if (k ~= 1)
         tpilotCarriers = [tpilotCarriers ',' num2str(pilotCarriers(k))]; 
    else
         tpilotCarriers = [tpilotCarriers num2str(pilotCarriers(k))]; 
    end
end
tpilotCarriers = [tpilotCarriers '),)</value>'];
txTemplate{513} = tpilotCarriers;

tpilotSym = '<value>((';
for k = 1:length(pilotSymbols)
    if (k ~= 1)
         tpilotSym = [tpilotSym ',' num2str(pilotSymbols(k))]; 
    else
         tpilotSym = [tpilotSym num2str(pilotSymbols(k))]; 
    end
end
tpilotSym = [tpilotSym '),)</value>'];
txTemplate{540} = tpilotSym;

toccupiedCarriers = '<value>([';
for k = 1:length(occupiedCarriers)
    if (k ~= 1)
         toccupiedCarriers = [toccupiedCarriers ',' num2str(occupiedCarriers(k))];
    else
         toccupiedCarriers = [toccupiedCarriers num2str(occupiedCarriers(k))];
    end
end
toccupiedCarriers = [toccupiedCarriers '],)</value>'];
txTemplate{378} = toccupiedCarriers;

txTemplate{813} = ['<value>"' packetHeadModType '"</value>'];
txTemplate{841} = ['<value>"' payloadModType '"</value>'];

[fileName,path] = uiputfile('*','Select file to save debug TCVR');
if (fileName ~= 0)
    f = fopen(horzcat(path,fileName),'w');
    fprintf(f,'%s\n',txTemplate{:});
    fclose(f);
end


% --- Executes on key press with focus on butSetOutput and none of its controls.
function butSetOutput_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to butSetOutput (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in debugCompareSignalsBut.
function debugCompareSignalsBut_Callback(hObject, eventdata, handles)
% hObject    handle to debugCompareSignalsBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUIvsGnuradio;



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over dataInputFileRadio.
function dataInputFileRadio_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to dataInputFileRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fs = str2double(get(handles.noDeviceSampleRate, 'string'))*1000;
fc = 7000000;
txGain = 0;
if (get(handles.usrpDevice, 'Value') == 1)
    oversampling = str2double(get(handles.oversamplingInput, 'string'));
    fs = str2double(get(handles.USRPSampleRate, 'string'))*1000/oversampling;
    fc = str2double(get(handles.carrierFreqInput, 'string'))*1000;
    txGain = str2double(get(handles.txGain, 'string'));
    deviceType = 1;
end;
if (get(handles.lw410Device, 'Value') == 1)
    fs = str2double(get(handles.lw410SampleRate, 'string'))*1000;
    fc = str2double(get(handles.lw410carrier, 'string'))*1000;
    deviceType = 2;
end;

sync1 = eval(get(handles.sync1, 'string'));
sync2 = eval(get(handles.sync2, 'string'));
ifftLen = str2num(get(handles.ifftLen, 'string'));
cpLen = str2double(get(handles.prefixLengthInput, 'string'));
packetLen = str2num(get(handles.packetLen, 'string'));
pilotSymbols = eval(get(handles.pilotSymbols, 'string'));
pilotCarriers = eval(get(handles.pilotTonesPositions, 'string'));
occupiedCarriers = eval(get(handles.occupiedCarriers, 'string'));

packetHeadModTypeList = get(handles.packetHeadModType,'String');
modTypeList = get(handles.modulationTypePopup,'String');
payloadModType = modTypeList{get(handles.modulationTypePopup, 'Value')};
packetHeadModType = packetHeadModTypeList(get(handles.packetHeadModType, 'Value'),:);

txTemplate = regexp( fileread('templates/OFDM_Transmitter.grc'),'\n','split');

txTemplate{567} = ['<value>' num2str(fs) '</value>'];
txTemplate{108} = ['<value>' num2str(fc) '</value>'];
txTemplate{189} = ['<value>' num2str(txGain) '</value>'];
txTemplate{216} = ['<value>' num2str(cpLen) '</value>'];

tsync1 = '<value>(';
for k = 1:length(sync1)
    if (k ~= 1)
         tsync1 = [tsync1 ',' num2str(sync1(k))]; 
    else
         tsync1 = [tsync1 num2str(sync1(k))]; 
    end
end
tsync1 = [tsync1 ')</value>'];
txTemplate{594} = tsync1;

tsync2 = '<value>(';
for k = 1:length(sync2)
    if (k ~= 1)
         tsync2 = [tsync2 ',' num2str(sync2(k))]; 
    else
         tsync2 = [tsync2 num2str(sync2(k))];
    end
end
tsync2 = [tsync2 ')</value>'];
txTemplate{621} = tsync2;

txTemplate{243} = ['<value>' num2str(ifftLen) '</value>'];
txTemplate{405} = ['<value>' num2str(packetLen) '</value>'];

tpilotCarriers = '<value>((';
for k = 1:length(pilotCarriers)
    if (k ~= 1)
         tpilotCarriers = [tpilotCarriers ',' num2str(pilotCarriers(k))]; 
    else
         tpilotCarriers = [tpilotCarriers num2str(pilotCarriers(k))]; 
    end
end
tpilotCarriers = [tpilotCarriers '),)</value>'];
txTemplate{513} = tpilotCarriers;

tpilotSym = '<value>((';
for k = 1:length(pilotSymbols)
    if (k ~= 1)
         tpilotSym = [tpilotSym ',' num2str(pilotSymbols(k))]; 
    else
         tpilotSym = [tpilotSym num2str(pilotSymbols(k))];
    end
end
tpilotSym = [tpilotSym '),)</value>'];
txTemplate{540} = tpilotSym;

toccupiedCarriers = '<value>([';
for k = 1:length(occupiedCarriers)
    if (k ~= 1)
         toccupiedCarriers = [toccupiedCarriers ',' num2str(occupiedCarriers(k))];
    else
         toccupiedCarriers = [toccupiedCarriers num2str(occupiedCarriers(k))];
    end
end
toccupiedCarriers = [toccupiedCarriers '],)</value>'];
txTemplate{378} = toccupiedCarriers;

txTemplate{880} = ['<value>"' packetHeadModType '"</value>'];
txTemplate{908} = ['<value>"' payloadModType '"</value>'];

[fileName,path] = uiputfile('*','Select file to save debug TCVR');
f = fopen(horzcat(path,fileName),'w');
fprintf(f,'%s\n',txTemplate{:});
fclose(f);


% --- Executes on button press in renewSyncBut.
function renewSyncBut_Callback(hObject, eventdata, handles)
% hObject    handle to renewSyncBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ifftLen = str2double(get(handles.ifftLen, 'string'));
occupiedCarriers = eval(get(handles.occupiedCarriers, 'string'));
pilotCarriers = eval(get(handles.pilotTonesPositions, 'string'));
generateSync(ifftLen,occupiedCarriers,pilotCarriers);

sync1 = '';
f = fopen('sync1.txt','r');
sync1 = fscanf(f,'%s');
fclose(f);

sync2 = '';
f = fopen('sync2.txt','r');
sync2 = fscanf(f,'%s');
fclose(f);

sync1(1) = '[';
sync1(end) = ']';
sync2(1) = '[';
sync2(end) = ']';

set(handles.sync1,'String',sync1);
set(handles.sync2,'String',sync2);

guidata(hObject, handles);


% --- Executes on button press in genLoopbackBut.
function genLoopbackBut_Callback(hObject, eventdata, handles)
% hObject    handle to genLoopbackBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oversampling = str2double(get(handles.noDeviceOversamplingInput, 'string'));
fs = str2double(get(handles.noDeviceSampleRate, 'string'))*1000/oversampling;
if (get(handles.usrpDevice, 'Value') == 1)
    oversampling = str2double(get(handles.oversamplingInput, 'string'));
    fs = str2double(get(handles.USRPSampleRate, 'string'))*1000/oversampling;
    deviceType = 1;
end;
if (get(handles.lw410Device, 'Value') == 1)
    fs = str2double(get(handles.lw410SampleRate, 'string'))*1000;
    deviceType = 2;
end;

sync1 = eval(get(handles.sync1, 'string'));
sync2 = eval(get(handles.sync2, 'string'));
ifftLen = str2num(get(handles.ifftLen, 'string'));
cpLen = str2double(get(handles.prefixLengthInput, 'string'));
packetLen = str2num(get(handles.packetLen, 'string'));
pilotSymbols = eval(get(handles.pilotSymbols, 'string'));
pilotCarriers = eval(get(handles.pilotTonesPositions, 'string'));
occupiedCarriers = eval(get(handles.occupiedCarriers, 'string'));

packetHeadModTypeList = get(handles.packetHeadModType,'String');
modTypeList = get(handles.modulationTypePopup,'String');
payloadModType = modTypeList{get(handles.modulationTypePopup, 'Value')};
packetHeadModType = packetHeadModTypeList(get(handles.packetHeadModType, 'Value'),:);

txTemplate = regexp( fileread('templates/OFDM_Test_Loopback.grc'),'\n','split');

txTemplate{567} = ['<value>' num2str(fs) ' </value>'];          %sample rate
txTemplate{216} = ['<value>' num2str(cpLen) ' </value>'];       %sample rate

tsync1 = '<value>(';
for k = 1:length(sync1)
    if (k ~= 1)
         tsync1 = [tsync1 ',' num2str(sync1(k))]; 
    else
         tsync1 = [tsync1 num2str(sync1(k))]; 
    end
end
tsync1 = [tsync1 ')</value>'];
txTemplate{594} = tsync1;                                       %sync1

tsync2 = '<value>(';
for k = 1:length(sync2)
    if (k ~= 1)
         tsync2 = [tsync2 ',' num2str(sync2(k))]; 
    else
         tsync2 = [tsync2 num2str(sync2(k))]; 
    end
end
tsync2 = [tsync2 ')</value>'];                                  
txTemplate{621} = tsync2;                                       %sync2

txTemplate{243} = ['<value>' num2str(ifftLen) '</value>'];
txTemplate{405} = ['<value>' num2str(packetLen) '</value>'];

tpilotSym = '<value>((';
for k = 1:length(pilotCarriers)
    if (k ~= 1)
         tpilotSym = [tpilotSym ',' num2str(pilotCarriers(k))]; 
    else
         tpilotSym = [tpilotSym num2str(pilotCarriers(k))]; 
    end
end
tpilotSym = [tpilotSym '),)</value>'];
txTemplate{513} = tpilotSym;

tpilotCarriers = '<value>((';
for k = 1:length(pilotSymbols)
    if (k ~= 1)
         tpilotCarriers = [tpilotCarriers ',' num2str(pilotSymbols(k))]; 
    else
         tpilotCarriers = [tpilotCarriers num2str(pilotSymbols(k))]; 
    end
end
tpilotCarriers = [tpilotCarriers '),)</value>'];
txTemplate{540} = tpilotCarriers;

toccupiedCarriers = '<value>([';
for k = 1:length(occupiedCarriers)
    if (k ~= 1)
         toccupiedCarriers = [toccupiedCarriers ',' num2str(occupiedCarriers(k))];
    else
         toccupiedCarriers = [toccupiedCarriers num2str(occupiedCarriers(k))];
    end
end
toccupiedCarriers = [toccupiedCarriers '],)</value>'];
txTemplate{378} = toccupiedCarriers;

txTemplate{982} = ['<value>"' packetHeadModType '"</value>'];
txTemplate{1069} = ['<value>"' packetHeadModType '"</value>'];
txTemplate{1010} = ['<value>"' payloadModType '"</value>'];
txTemplate{1097} = ['<value>"' payloadModType '"</value>'];

[fileName,path] = uiputfile('*','Select file to save debug TCVR');
f = fopen(horzcat(path,fileName),'w');
fprintf(f,'%s\n',txTemplate{:});
fclose(f);


% --- Executes on button press in translateBut.
function translateBut_Callback(hObject, eventdata, handles)
% hObject    handle to translateBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.langENG)
    handles.langENG = 0
    set(handles.dataInputGroup ,'Title','Vstupni data');
    set(handles.basicSettings ,'Title','Zakladni nastaveni');
    set(handles.prefixRadioGroup ,'Title','Prefix/Ochranny interval');
    set(handles.subcarriersMapping ,'Title','Mapovani subnosnych a pilotnich tonu');
    set(handles.packetPanel ,'Title','Nastaveni paketu');
    set(handles.chanelParamsGroup ,'Title','Nastaveni kanalu');
    set(handles.hintsPanel ,'Title','Napoveda');
    set(handles.deviceSettingsPanel ,'Title','Nastaveni vystupniho zarizeni');
    set(handles.connectedDeviceGroup ,'Title','pripojene zarizeni');
    set(handles.noDeviceSettings ,'Title','Nastaveni bez zarizeni');
    set(handles.usrpSettings ,'Title','USRP parametry');
    set(handles.lw410Params ,'Title','LW410 nastaveni');
    set(handles.signalOutputPanel ,'Title','Vystupni signal');
    set(handles.dataInputRandomRadio ,'String','nahodni sekvence');
    set(handles.dataInputUserRadio ,'String','uziv. retezec');
    set(handles.dataInputFileRadio ,'String','ze souboru');
    set(handles.modulationTypeText ,'String','Typ modulace dat');
    set(handles.numIFFTText ,'String','IFFT delka');
    set(handles.textSignalAmpl ,'String','Amplituda signalu');
    set(handles.symbolLengthLabel ,'String','Pocet symbolu');
    set(handles.intervalRadio ,'String','Cyklicky prefix');
    set(handles.prefixLengthText ,'String','Delka prefixu');
    set(handles.guardRadio ,'String','Ochranny interval');
    set(handles.guardIntervalText ,'String','Ochranny int. [vlevo vpravo]');
    set(handles.intervalNotusedRadio ,'String','Nepouzit');
    set(handles.subcarrierMappingText ,'String','Mapovani subnosnych');
    set(handles.pilotPositionsText ,'String','Pozice pilotnich tonu');
    set(handles.pilotAmplitudesText ,'String','Amplituda pilotnich tonu');
    set(handles.packetHeadModText ,'String','Modulace hlavicek paketu');
    set(handles.packetLenText ,'String','Velikost paketu (byte)');
    set(handles.enablePackets ,'String','pouzit pakety');
    set(handles.sync1Text ,'String','1. synchronizacni slovo');
    set(handles.sync2Text ,'String','2. synchronizacni slovo');
    set(handles.renewSyncBut ,'String','generovat');
    set(handles.noChannelRadio ,'String','idealni prenos');
    set(handles.cyclicGeneration ,'String','generovat cyklicky');
    set(handles.noDevice ,'String','bez zarizeni');
    set(handles.sampleRateText ,'String','Vzorkovaci kmitocet (kHz)');
    set(handles.noDeviceOversamplingText ,'String','Prevzorkovani');
    set(handles.usrpSampleRateText ,'String','Vzorkovaci kmitocet (kHz)');
    set(handles.txCycleCountText ,'String','Pocet opakovani');
    set(handles.txGainText ,'String','TX zisk');
    set(handles.lw410SampleRateText ,'String','Vzork. kmitocet (kHz)');
    set(handles.gpibAddrLabel,'String','GPIB adresa');
    set(handles.sampleFreqLabel ,'String','Vzorkovaci frekvence');
    set(handles.carrierFreqText ,'String','Kmitocet nosne (kHz)');
    set(handles.usrpOversamplingText ,'String','Prevzorkovani');    
    set(handles.lw410carrierText ,'String','Kmitocet nosne (kHz)');
else
    handles.langENG = 1;
    set(handles.dataInputGroup ,'Title','Data Input');
    set(handles.basicSettings ,'Title','Basic Settings');
    set(handles.prefixRadioGroup ,'Title','Prefix/Guard Band');
    set(handles.subcarriersMapping ,'Title','Subcarriers and pilots mapping');
    set(handles.packetPanel ,'Title','Packet Settings');
    set(handles.chanelParamsGroup ,'Title','Packet Settings');
    set(handles.hintsPanel ,'Title','Hints');
    set(handles.deviceSettingsPanel ,'Title','Device Settings');
    set(handles.connectedDeviceGroup ,'Title','Connected Device');
    set(handles.noDeviceSettings ,'Title','No device settings');
    set(handles.usrpSettings ,'Title','USRP params');
    set(handles.lw410Params ,'Title','LW410 Settings');
    set(handles.signalOutputPanel ,'Title','Signal output');
    set(handles.dataInputRandomRadio ,'String','random seq.');
    set(handles.dataInputUserRadio ,'String','user string');
    set(handles.dataInputFileRadio ,'String','from file');
    set(handles.modulationTypeText ,'String','Data Modulation Type');
    set(handles.numIFFTText ,'String','IFFT length');
    set(handles.textSignalAmpl ,'String','Signal amplitude');
    set(handles.symbolLengthLabel ,'String','Number of symbols');
    set(handles.intervalRadio ,'String','Cyclic prefix');
    set(handles.prefixLengthText ,'String','Prefix length');
    set(handles.guardRadio ,'String','Guard interval');
    set(handles.guardIntervalText ,'String','Guard int. [left right]');
    set(handles.intervalNotusedRadio ,'String','Not used');
    set(handles.subcarrierMappingText ,'String','Subcarriers mapping');
    set(handles.pilotPositionsText ,'String','Pilot tones positions');
    set(handles.pilotAmplitudesText ,'String','Pilot symbols');
    set(handles.packetHeadModText ,'String','Packet Head Modulation');
    set(handles.packetLenText ,'String','Packet length (bytes)');
    set(handles.enablePackets ,'String','enable packets');
    set(handles.sync1Text ,'String','1st synchronization word');
    set(handles.sync2Text ,'String','2nd synchronization word');
    set(handles.renewSyncBut ,'String','renew');
    set(handles.noChannelRadio ,'String','no channel');
    set(handles.cyclicGeneration ,'String','cyclic generation');
    set(handles.noDevice ,'String','no device');
    set(handles.sampleRateText ,'String','Sample rate (kHz)');
    set(handles.noDeviceOversamplingText ,'String','Oversampling');
    set(handles.usrpSampleRateText ,'String','Sample rate(kHz)');
    set(handles.usrpOversamplingText ,'String','Oversampling');    
    set(handles.txCycleCountText ,'String','Tx cycle count');
    set(handles.txGainText ,'String','TX gain');
    set(handles.lw410SampleRateText ,'String','Sample Rate (kHz)');
    set(handles.gpibAddrLabel,'String','GPIB addres');
    set(handles.sampleFreqLabel ,'String','Sample freq.');
    set(handles.carrierFreqText ,'String','Carrier freq. (kHz)');
    set(handles.lw410carrierText ,'String','Carrier freq. (kHz)');
end
guidata(hObject, handles);


% --- Executes on button press in dataInputRandomRadio.
function dataInputRandomRadio_Callback(hObject, eventdata, handles)
% hObject    handle to dataInputRandomRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dataInputRandomRadio
set(handles.seqInput,'Enable','off');

enablePackets = get(handles.enablePackets, 'Value');
packetLen = str2double(get(handles.packetLen,'string'));
modTypeList = get(handles.modulationTypePopup,'String');
payloadModType = modTypeList{get(handles.modulationTypePopup, 'Value')}
occupiedCarriers = eval(get(handles.occupiedCarriers, 'string'));
nSymbols = str2num(get(handles.numOfSymbolsInput, 'string'));
nCarriers = length(occupiedCarriers);
if (enablePackets)
    packetLen = get(handles.packetLen, 'string');
    hintsText = {...
        'Input data source: random sequence, ' ...
        ['Data modulation type: ' payloadModType] ...
        ['Number of Carriers: ' num2str(nCarriers)]...
        ['Packet Length: ' packetLen] ...
        ['Packet Count: ' num2str(nSymbols)] ...
        '' ...  
        'Input data are generated automatically as random bytes:' ...
        ['    rand(255,1,' packetLen '*' num2str(nSymbols) ')'] ...
    };
else
    if (strcmp(payloadModType, 'QPSK'))
        M = 4
    elseif (strcmp(payloadModType, 'BPSK'))
        M = 2
    elseif (strcmp(payloadModType, 'PSK-8'))
        M = 8
    end
    hintsText = {...
        'Input data source: random sequence, ' ...
        ['Data modulation type: ' payloadModType] ...
        ['Number of Carriers: ' num2str(nCarriers)]...
        ['Generated Input Data Length (bits 0,1): ' num2str(nCarriers) '*log2(' num2str(M) ')*' num2str(nSymbols) ' = ' num2str(nCarriers*log2(M)*nSymbols)] ...
        '' ...
        'Input data are generated automatically as random sequence:' ...
        '    randi([0 1],1,dataLen)' ...
    };
end
set(handles.hintsText,'String',hintsText);

guidata(hObject, handles);


% --- Executes on button press in awgnRadio.
function awgnRadio_Callback(hObject, eventdata, handles)
% hObject    handle to awgnRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of awgnRadio
set(handles.awgnLevel,'Enable','on');
guidata(hObject, handles);


% --- Executes on button press in noChannelRadio.
function noChannelRadio_Callback(hObject, eventdata, handles)
% hObject    handle to noChannelRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noChannelRadio
set(handles.awgnLevel,'Enable','off');
guidata(hObject, handles);


% --- Executes on button press in guardRadio.
function guardRadio_Callback(hObject, eventdata, handles)
% hObject    handle to guardRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of guardRadio
set(handles.prefixLengthInput,'Enable','off');
set(handles.guardLengthInput,'Enable','on');
guidata(hObject, handles);

%
%   Change hints text to have it actual, this is called inside functions
%   which change informations.
%
function changeHintsText(hObject, eventdata, handles)
if (get(handles.dataInputRandomRadio, 'Value') == 1)
    dataInputRandomRadio_Callback(hObject, eventdata, handles)
elseif (get(handles.dataInputUserRadio, 'Value') == 1)
    dataInputUserRadio_Callback(hObject, eventdata, handles)
elseif (get(handles.dataInputFileRadio, 'Value') == 1)
    dataInputFileRadio_Callback(hObject, eventdata, handles)
end
guidata(hObject, handles);    

% --- Executes on button press in showFramesConstellationBut.
function showFramesConstellationBut_Callback(hObject, eventdata, handles)
% hObject    handle to showFramesConstellationBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hChildren = get(handles.ofdmWaveAxes3,'Children');
sigReal = get(hChildren(2),'YData');
sigImag = get(hChildren(1),'YData');

handles.payloadMod.constellation;
hold on;
plot(sigReal, sigImag);
hold off;



function noDeviceOversamplingInput_Callback(hObject, eventdata, handles)
% hObject    handle to noDeviceOversamplingInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noDeviceOversamplingInput as text
%        str2double(get(hObject,'String')) returns contents of noDeviceOversamplingInput as a double


% --- Executes during object creation, after setting all properties.
function noDeviceOversamplingInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noDeviceOversamplingInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lw410SampleRate_Callback(hObject, eventdata, handles)
% hObject    handle to lw410SampleRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lw410SampleRate as text
%        str2double(get(hObject,'String')) returns contents of lw410SampleRate as a double


% --- Executes during object creation, after setting all properties.
function lw410SampleRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lw410SampleRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lw410carrier_Callback(hObject, eventdata, handles)
% hObject    handle to lw410carrier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lw410carrier as text
%        str2double(get(hObject,'String')) returns contents of lw410carrier as a double


% --- Executes during object creation, after setting all properties.
function lw410carrier_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lw410carrier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in connectedDeviceGroup.
function connectedDeviceGroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in connectedDeviceGroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.outputDeviceImg);
if (get(handles.usrpDevice, 'Value') == 1)
    set(handles.outputDeviceImg,'Visible','on');
    image(handles.usrpDeviceImg);
    set(handles.outputDeviceImg,'Visible','off');
elseif (get(handles.lw410Device, 'Value') == 1)
    set(handles.outputDeviceImg,'Visible','on');
    image(handles.lw410Img);
    set(handles.outputDeviceImg,'Visible','off');
elseif (get(handles.noDevice, 'Value') == 1)
    set(handles.outputDeviceImg,'Visible','off');
    set(get(handles.outputDeviceImg,'children'),'visible','off')
end

% --- Executes on button press in lw410Device.
function usrpDevice_Callback(hObject, eventdata, handles)
% hObject    handle to lw410Device (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lw410Device

% --- Executes on button press in lw410Device.
function lw410Device_Callback(hObject, eventdata, handles)
% hObject    handle to lw410Device (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lw410Device
