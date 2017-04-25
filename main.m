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

% Last Modified by GUIDE v2.5 26-Apr-2017 00:38:13

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
set(handles.ofdmWavePanel,'visible','on')
set(handles.deviceSettingsPanel, 'visible','off');

handles.lw410 = LW410Interface();
handles.ofdm.dataInputMethod = 'randsrc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%26.4.2017
handles.ofdm.headData = 0;

% Update handles structure
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

guidata(hObject, handles);  %how to obtain data


% --- Executes on button press in deviceSetingsBut.
function deviceSetingsBut_Callback(hObject, eventdata, handles)
% hObject    handle to deviceSetingsBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data.activePanel = 'deviceSettingsPanel';   %set semaphor
set(handles.deviceSettingsPanel, 'visible','on');
set(handles.ofdmWavePanel,'visible','off')

guidata(hObject, handles);  %how to obtain data


% --- Executes on button press in generateOFDMWaveBut.
function generateOFDMWaveBut_Callback(hObject, eventdata, handles)
% hObject    handle to generateOFDMWaveBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USRP Transmission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % 
% % for counter = 1:20
% %   data = randi([0 1], 30, 1);
% %   modSignal = step(hMod, data);
% % end
% % 

handles.usrp.interp = str2double(get(handles.editUsrpInterp, 'string'))
handles.lw410.fs = 100e6/handles.usrp.interp; %zatial napevno!
fs = handles.lw410.fs;      %NEZABUDNUT NASTAVOVAT!

%--------------------------------------------------------------------------
% READING DATA FOR GENERATING WAVES FROM GUI by LuboJ.
%--------------------------------------------------------------------------
signalAmplitude = str2double(get(handles.editSignalAmplitude, 'string'))
'Signal amplitude set to: '
signalAmplitude

numCarriers = str2double(get(handles.numOfCarriersInput, 'string'))
numSymbols = str2double(get(handles.numOfSymbolsInput, 'string'))

cpLen = str2double(get(handles.intervalLengthInput, 'string'))
oversampleFactor = str2double(get(handles.oversamplingInput, 'string'));
fc = str2double(get(handles.carrierFreqInput, 'string'))*1e3
                                                             
%
%   Typ modulacie
%

modTypeList = get(handles.modulationTypePopup,'String');
modType = modTypeList{get(handles.modulationTypePopup, 'Value')}
if (strcmp(modType, 'QPSK'))
    M = 4
elseif (strcmp(modType, 'pi/4 DQPSK'))
    M = 4
elseif (strcmp(modType, '16QAM'))
    M = 16
else
    M = 0
end
numBitsInGroup = log10(M)/log10(2)
numDataPoints = numCarriers * numSymbols * numBitsInGroup;

%
%   Zdroj vstupnych dat
%
if (strcmp(handles.ofdm.dataInputMethod,'file'))
    data_source = dlmread(handles.ofdm.file);
elseif (strcmp(handles.ofdm.dataInputMethod,'randsrc'))
    data_source = randsrc(1, numDataPoints, [0 1]);
elseif (strcmp(handles.ofdm.dataInputMethod,'user'))
    data_source = eval(get(handles.seqInput, 'string'));
    data_source_aux = [];
    for i = 1:length(data_source)
        if (data_source(i) == 0) data_source_aux = [data_source_aux 0 0];
        elseif (data_source(i) == 1) data_source_aux = [data_source_aux 0 1];
        elseif (data_source(i) == 2) data_source_aux = [data_source_aux 1 0];
        elseif (data_source(i) == 3) data_source_aux = [data_source_aux 1 1];
        end
    end
    data_source = data_source_aux;
end
%handles.data.seqData = data_source;

%
%   Nastavenie pilotnych tonov.
%
pilotTonesPositions = eval(get(handles.pilotTonesPositions, 'string'));
pilotTonesAmplitudes = eval(get(handles.pilotTonesAmplitudes, 'string'));
if (length(pilotTonesPositions) > 0)
    usePilots = 1;
else
    usePilots = 0;
end

%
%   Head data, reading as binary from file.
%
if (handles.ofdm.headData == 1)
    f = fopen(handles.ofdm.headFile,'r')
    headData = abs(fread(f,'bit1'));
    fclose(f);
    numHeadData = length(headData)
    % headData
end

%
%   Vytvorenie modulatoru
%       handles.hModulator
%
clear modulated_data;
if (strcmp(modType, 'QPSK'))
    handles.hModulator =  comm.QPSKModulator('BitInput',true,'SymbolMapping','Binary');
elseif (strcmp(modType, 'pi/4 DQPSK'))
    handles.hModulator =  comm.DQPSKModulator(pi/4, 'BitInput', true);
elseif (strcmp(modType, '16QAM'))
    handles.hModulator = comm.RectangularQAMModulator(16,'BitInput',true);
end

%  Vypocet pociatku cyklickeho prefixu
cpStart = (numCarriers-cpLen)*oversampleFactor

%   Nastavenie parametrov podla GUI.

insertGuardInterval = get(handles.intervalGuardIntervalRadio, 'Value');

awgnLevel = str2double(get(handles.awgnLevel, 'string'));
if (get(handles.noChannelRadio, 'Value') == 1)
    useChannelModel = 0;
elseif (get(handles.awgnRadio, 'Value') == 1)
    useChannelModel = 1;
elseif (get(handles.rayleighRadio, 'Value') == 1)
    maxDoplerShift = str2double(get(handles.doplerShiftIn, 'string'));
    useChannelModel = 2;
    channelModel = rayleighchan(1/fs, maxDoplerShift);
elseif (get(handles.ricianRadio, 'Value') == 1)
    useChannelModel = 3;    
end

% USRP inicializacia

% % % % % % Tx = comm.SDRuTransmitter(...
% % % % % %   'Platform','N200/N210/USRP2',...
% % % % % %   'IPAddress', '192.168.10.2', ...
% % % % % %   'CenterFrequency', fc, ...
% % % % % %   'InterpolationFactor', handles.usrp.interp);
% % % % % % hMod = comm.DPSKModulator('BitInput',true);

handles.cyclicGeneration

nfor = str2double(get(handles.editTxCycleCount, 'string'))      %defaultne 100krat prebehne tx
benchmark = zeros(1,nfor);         %meranie rychlosti
cyclicGeneration = 1;
while (cyclicGeneration == 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Jadro, OFDM modulator
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    for cycleCnt = 0:1:nfor
        disp(horzcat('generate samples, step ',num2str(cycleCnt)));
        tic;

        %%%%%%
        % OFDM modulator
        %   data_source <--- VSTUPNE DATA
        %
        %   hodnota vstupnych vzoriek pre USRP je napatie
        %
        if (handles.ofdm.headData == 1)
            data_source = randi([0 1], numCarriers, numSymbols);    
            data_source = vertcat(repmat(headData,1,numSymbols), data_source);
            data_source = reshape(data_source, (numCarriers+numHeadData)*numSymbols, 1);                
            modulated_data = signalAmplitude * step(handles.hModulator, data_source);
            % have to divide M, because binary data get compressed by
            % modulation
            data_matrix = reshape(modulated_data, (numCarriers+numHeadData)/log2(M), numSymbols);        
        else
            data_source = randi([0 1], length(data_source), 1);    
            modulated_data = signalAmplitude * step(handles.hModulator, data_source);
            data_matrix = reshape(modulated_data, numCarriers, numSymbols);        
        end
        
        %
        % Pilot tones           ZLE, SKONTROLOVAT CI SU VYPOZICIOVANE!
        %
        if (usePilots)
            for k = 1:length(pilotTonesPositions)
                    data_matrix = vertcat( ...
                        data_matrix(1:pilotTonesPositions(k)-1+k-1,:), ...
                        pilotTonesAmplitudes(k)*ones(1,numSymbols), ...
                        data_matrix(pilotTonesPositions(k)+k-1:end,:) ...
                    );
            end
        end
        
        %
        % Upsampling, in the end data are MULTIPLIED TO CONSERVE ENERGY
        %
        if (oversampleFactor ~= 1)
            ifft_data = ifft( ...
                    [data_matrix(1:end/2,:); ...
                    zeros((oversampleFactor-1)*numCarriers, numSymbols); ...
                    data_matrix(end/2+1:end,:)], ...
                    oversampleFactor*numCarriers) * oversampleFactor;
        else
            ifft_data = ifft(data_matrix,numCarriers);
        end
                
        %
        % Guard interval
        %
        if ( insertGuardInterval == 1)
            ifft_data = vertcat( ...
                zeros(cpLen/2*oversampleFactor, numSymbols), ...
                ifft_data, ...
                zeros(cpLen/2*oversampleFactor, numSymbols));
        elseif (get(handles.intervalCyclicPrefixRadio, 'Value') == 1)
            ifft_data = vertcat(ifft_data(cpStart+1:end,:), ifft_data);
        else
            %%%%NOTHING%%%%%
        end
    
        %Serialization
        ofdm_signal = reshape(ifft_data, numel(ifft_data), 1);

        %
        %   USRP Tx
        %
        if (useChannelModel == 1)
            ofdm_signal_noisy = awgn(ofdm_signal, awgnLevel, 20*log10(sum(abs(ofdm_signal))));
        elseif (useChannelModel == 2)
            ofdm_signal_noisy = filter(channelModel, ofdm_signal);
        elseif (useChannelModel == 3)
    %         ofdm_signal_noisy = filter(channelModel, ofdm_signal);
        else                                  %useChannel == 0, default, no noise
                                              %no change of signal, ideal case
              ofdm_signal_noisy = ofdm_signal;
        end
    % % % % % %     step(Tx, ofdm_signal_noisy);
        benchmark(cycleCnt+1) = toc;       %meranie casu vypoctu
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Interaktivna cast
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pause(0.2);
    cyclicGeneration = get(handles.cyclicGeneration,'value');

    signalAmplitude = str2double(get(handles.editSignalAmplitude, 'string'));
    awgnLevel = str2double(get(handles.awgnLevel, 'string'));
    if (get(handles.noChannelRadio, 'Value') == 1)
        useChannelModel = 0;
    elseif (get(handles.awgnRadio, 'Value') == 1)
        useChannelModel = 1;
    elseif (get(handles.rayleighRadio, 'Value') == 1)
        maxDoplerShift = str2double(get(handles.doplerShiftIn, 'string'));
        useChannelModel = 2;
        channelModel = rayleighchan(1/fs, maxDoplerShift);
    elseif (get(handles.ricianRadio, 'Value') == 1)
        useChannelModel = 3;    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Signal in time plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    axes(handles.ofdmWaveAxes);
    plot(zeros(1,numCarriers));
    title('OFDM signal');
    xlabel('sample[-]');
    ylabel('A[V]');
    plot(linspace(0,handles.lw410.fs,length(ofdm_signal)), real(ofdm_signal));
    hold on; plot(linspace(0,handles.lw410.fs,length(ofdm_signal)), imag(ofdm_signal),'r'); hold off;
    
    plot(linspace(0,handles.lw410.fs,length(ofdm_signal)), real(ofdm_signal));
    title('OFDM signal in time domain');
    hold on; plot(linspace(0,handles.lw410.fs,length(ofdm_signal)), imag(ofdm_signal),'r'); hold off;
    grid on;
    xlabel('t[s]'); ylabel('A[V]');

    disp('Energia OFDM signalu (dBm)');
    10*log10(sum(abs(ofdm_signal).^2))

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Signal after channel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sigSpec = 20*log10(abs(fft(ofdm_signal, 2048))./2048);
    minX = -handles.lw410.fs/2;               %VYPOCET FREKVENCNEJ OSI PRE ZOBRAZENIE!
                                                %SKONTROLOVAT!
    maxX = -minX;
    xAxis = linspace(minX, maxX, 2048);

    axes(handles.ofdmWaveAxes4);
    sigSpecNoisy = 20*log10(abs(fft(ofdm_signal_noisy, 2048))./2048);
    plot(xAxis, fftshift(sigSpecNoisy));
    ylim([-100 -40]);
    xlim([-fs/2 fs/2]);
    hold on;
    plot(xAxis,fftshift(sigSpec), '--r')
    hold off;
    if (useChannelModel == 1)
        title('Signal after AWGN channel'); xlabel('f[kHz]'); ylabel('A[dBm]');
    elseif (useChannelModel == 2)
        title('Signal after rayleigh channel'); xlabel('f[Hz]'); ylabel('A[dBm]');
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
                'USRP generating samples. no channel. ', ...
                'Signal energy ', ...
                num2str(ofdm_signal_energy),'dBm.', ...
                'Total channel energy ', ...
                num2str(ofdm_signal_energy_noisy),'dBm.' ...
            ));
    elseif (useChannelModel == 1)   % AWGN kanal parametre vypis
        set(handles.statusText, 'String', ...
            horzcat(    ...
                'USRP generating samples. Using AWGN channel. ', ...
                'Signal energy ', ...
                num2str(ofdm_signal_energy),'dBm.', ...
                'Total channel energy ', ...
                num2str(ofdm_signal_energy_noisy),'dBm.' ...
            ));
    elseif (useChannelModel == 2)   % rayleigh kanal parametre vypis
        set(handles.statusText, 'String', ...
            horzcat(    ...
                'USRP generating samples. Using rayleigh channel. ', ...
                'Signal energy ', ...
                num2str(ofdm_signal_energy),'dBm.', ...
                'Total channel energy ', ...
                num2str(ofdm_signal_energy_noisy),'dBm.' ...
            ));
    else
    end
end     %koniec modulovania

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       KONIEC GENEROVANIA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dlmwrite('benchmark.mat',benchmark,';');    %zapis do suboru
handles.benchmark = benchmark;
% plot(dlmread('benchmark.mat',';'))        %vykreslenie

% 
% % PLOT - OFDM Phase spectrum
% phaseSpec = 180/pi*angle(fft(real(ofdm_signal), 2048));
% length(phaseSpec)
% axes(handles.ofdmWaveAxes4);
% plot(xAxis, phaseSpec);
% title('OFDM phase spectrum'); xlabel('f[kHz]'); ylabel('A[deg]');

axes(handles.benchmarkAxes);
plot(benchmark(5:end));
%title('Generating cycle run time (s)');
xlabel('cycle [-]'); ylabel('t[s]');
grid on;

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
    [fileName,path] = uigetfile('','Select file with samples');
    handles.ofdm.file = horzcat(path, fileName);
    set(handles.seqInput, 'string', horzcat(path, fileName));
end

disp(horzcat('Data input method changed to ', handles.ofdm.dataInputMethod));
guidata(hObject, handles);

% --- Executes on button press in butExportOFDMData.
function butExportOFDMData_Callback(hObject, eventdata, handles)
% hObject    handle to butExportOFDMData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filePath dir] = uiputfile('mat', 'Choose output file')
dlmwrite(horzcat(dir,filePath), handles.data.seqData, ',');

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



function numOfCarriersInput_Callback(hObject, eventdata, handles)
% hObject    handle to numOfCarriersInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numOfCarriersInput as text
%        str2double(get(hObject,'String')) returns contents of numOfCarriersInput as a double


% --- Executes during object creation, after setting all properties.
function numOfCarriersInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numOfCarriersInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function intervalLengthInput_Callback(hObject, eventdata, handles)
% hObject    handle to intervalLengthInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of intervalLengthInput as text
%        str2double(get(hObject,'String')) returns contents of intervalLengthInput as a double


% --- Executes during object creation, after setting all properties.
function intervalLengthInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to intervalLengthInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in intervalCyclicPrefixRadio.
function intervalCyclicPrefixRadio_Callback(hObject, eventdata, handles)
% hObject    handle to intervalCyclicPrefixRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of intervalCyclicPrefixRadio


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


% --- Executes on button press in showConstellationBut.
function showConstellationBut_Callback(hObject, eventdata, handles)
% hObject    handle to showConstellationBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% PLOT - Constellation diagram      FOR NOW TURN OFF!
modTypeList = get(handles.modulationTypePopup,'String');
modType = modTypeList{get(handles.modulationTypePopup, 'Value')}
if (strcmp(modType, 'QPSK'))
    handles.hModulator.constellation;
elseif (strcmp(modType, 'pi/4 DQPSK'))
    % NO CONSTELLATION
elseif (strcmp(modType, '16QAM'))
    handles.hModulator.constellation;
end



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


% --- Executes on button press in benchmarkButton.
function benchmarkButton_Callback(hObject, eventdata, handles)
% hObject    handle to benchmarkButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
scrsz = get(0,'ScreenSize');
h = figure('Tag','Benchmark','Name','OFDM signal','Position',[scrsz(3)/2-500 scrsz(4)/2-250 1000 500]);
plot(handles.benchmark(2:end));



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



function pilotTonesAmplitudes_Callback(hObject, eventdata, handles)
% hObject    handle to pilotTonesAmplitudes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pilotTonesAmplitudes as text
%        str2double(get(hObject,'String')) returns contents of pilotTonesAmplitudes as a double


% --- Executes during object creation, after setting all properties.
function pilotTonesAmplitudes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pilotTonesAmplitudes (see GCBO)
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
dlmwrite('OFDMGeneratorSettings.txt', handles.cyclicGeneration);

% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



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


% --- Executes on button press in benchmarkButton.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to benchmarkButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



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
    [fileName,path] = uigetfile('','Select file with head data');
    handles.ofdm.headFile = horzcat(path, fileName);
end
guidata(hObject, handles);
