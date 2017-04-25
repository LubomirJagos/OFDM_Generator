function varargout = LW410_Basic_generation(varargin)
% LW410_BASIC_GENERATION MATLAB code for LW410_Basic_generation.fig
%      LW410_BASIC_GENERATION, by itself, creates a new LW410_BASIC_GENERATION or raises the existing
%      singleton*.
%
%      H = LW410_BASIC_GENERATION returns the handle to a new LW410_BASIC_GENERATION or the handle to
%      the existing singleton*.
%
%      LW410_BASIC_GENERATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LW410_BASIC_GENERATION.M with the given input arguments.
%
%      LW410_BASIC_GENERATION('Property','Value',...) creates a new LW410_BASIC_GENERATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LW410_Basic_generation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LW410_Basic_generation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LW410_Basic_generation

% Last Modified by GUIDE v2.5 14-Dec-2016 11:01:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LW410_Basic_generation_OpeningFcn, ...
                   'gui_OutputFcn',  @LW410_Basic_generation_OutputFcn, ...
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

% --- Executes just before LW410_Basic_generation is made visible.
function LW410_Basic_generation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LW410_Basic_generation (see VARARGIN)

% Choose default command line output for LW410_Basic_generation
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


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LW410_Basic_generation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LW410_Basic_generation_OutputFcn(hObject, eventdata, handles) 
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

%--------------------------------------------------------------------------
% READING DATA FOR GENERATING WAVES FROM GUI by LuboJ.
%--------------------------------------------------------------------------
numCarriers = str2double(get(handles.numOfCarriersInput, 'string'))
numSymbols = str2double(get(handles.numOfSymbolsInput, 'string'))
no_of_ifft_points = numCarriers
cp_len = str2double(get(handles.intervalLengthInput, 'string'))
oversampling = str2double(get(handles.oversamplingInput, 'string'))
fc = str2double(get(handles.carrierFreqInput, 'string'))*1e3
                                                             
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
no_of_data_points = numCarriers * numSymbols * numBitsInGroup;

if (strcmp(handles.ofdm.dataInputMethod,'file'))
    data_source = dlmread(handles.ofdm.file);
elseif (strcmp(handles.ofdm.dataInputMethod,'randsrc'))
    data_source = randsrc(1, no_of_data_points, [0 1]);
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
handles.data.seqData = data_source;

%
%   2.  Perform modulation
%
clear modulated_data;
if (strcmp(modType, 'QPSK'))
    handles.hModulator =  comm.QPSKModulator('BitInput',true,'SymbolMapping','Binary');
elseif (strcmp(modType, 'pi/4 DQPSK'))
    handles.hModulator =  comm.DQPSKModulator(pi/4, 'BitInput', true);
elseif (strcmp(modType, '16QAM'))
    handles.hModulator = comm.RectangularQAMModulator(16,'BitInput',true);
end
modulated_data = step(handles.hModulator, data_source.');

%   Write modulated data into file for further comparison
dlmwrite('sendData.dat',modulated_data,',');

%   3.  Do IFFT on each block
%   Make the serial stream a matrix where each column represents a pre-OFDM
%   block (w/o cyclic prefixing)
data_matrix = reshape(modulated_data, numCarriers, numSymbols);

% Compute from where to take cyclic prefix
cp_start = numCarriers-cp_len; 
cp_end = numCarriers;

%   3. Columnwise IFFT and copy cyclic prefix.
%ifft_data = ifft(data_matrix,no_of_ifft_points*oversampling);
ifft_data = ifft(data_matrix,no_of_ifft_points);

if (get(handles.intervalGuardIntervalRadio, 'Value') == 1)
    ifft_data = vertcat(zeros(cp_len, numSymbols), ifft_data);
elseif (get(handles.intervalCyclicPrefixRadio, 'Value') == 1)
    ifft_data = vertcat(ifft_data(cp_start+1:end,:), ifft_data);
else
    % DO NOTHING!
end

%   4. Convert to serial stream for transmission
[rows_ifft_data cols_ifft_data]=size(ifft_data);
len_ofdm_data = rows_ifft_data*cols_ifft_data;

%   Actual OFDM signal to be transmitted
ofdm_signal = reshape(ifft_data, 1, len_ofdm_data);

% % % % % % % % % % % % % % % % % % figure;
% % % % % % % % % % % % % % % % % % plot(real(ofdm_signal), 'b-x'); hold on;
% % % % % % % % % % % % % % % % % % plot(imag(ofdm_signal), 'r-x'); hold off;


% IQ modulation
ofdm_signal = resample(ofdm_signal,oversampling,1);   %resample to have more samples per symbol
tc = [0:handles.lw410.sampletime:(length(ofdm_signal)-1)*handles.lw410.sampletime];

ofdm_signal = ofdm_signal .* (cos(2*pi*fc*tc) + j*sin(2*pi*fc*tc));
%i = 1.0i; ofdm_signal = ofdm_signal .* exp(i*2*pi*fc*tc);              %also IQ modulation

% This should insert square wave into signal begin to know at receiver where
% to start demodulate.
% Output level divided by 2 because for some reason levels from generator
% are higher, when measure INCREASE OUTPUT BANDWIDTH FILTER (CHAN 1 button)

dlmwrite('sendOFDMSignal.dat', ofdm_signal,',');
%ofdm_signal = ofdm_signal/max(abs(ofdm_signal))/2;          %NORMALIZATION
%ofdm_signal = real(ofdm_signal);
handles.data.seq = 2*real(ofdm_signal);     %multiply by 2 because there is just one part of compelx signal
dlmwrite('sendOFDMSeq.dat', handles.data.seq,',');

% PLOT - Real part of OFDM signal
axes(handles.ofdmWaveAxes);
    % For plotting there is multiply by 2, because it's only real part
    % without imaginary part.
    plot(4*real(ofdm_signal)); title('OFDM signal'); %LuboJ. NEED CORRECT check!
    xlabel('sample[-]'); ylabel('A[V]');

% PLOT - OFDM signal spectrum
axes(handles.ofdmWaveAxes3);
xAxis = linspace(0, 1/handles.lw410.sampletime/1e3, 2048);
plot(xAxis, 20*log10(abs(fft(real(ofdm_signal), 2048))));
title('OFDM magnitude spectrum'); xlabel('f[kHz]'); ylabel('A[-]');

% PLOT - OFDM Phase spectrum
phaseSpec = 180/pi*angle(fft(real(ofdm_signal), 2048));
length(phaseSpec)
axes(handles.ofdmWaveAxes4);
plot(xAxis, phaseSpec);
title('OFDM phase spectrum'); xlabel('f[kHz]'); ylabel('A[deg]');

set(handles.statusText, 'String', 'Sequence set to OFDM.');
guidata(hObject, handles);


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

% --- Executes on button press in showOfdmWaveAxes3But.
function showOfdmWaveAxes3But_Callback(hObject, eventdata, handles)
% hObject    handle to showOfdmWaveAxes3But (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
scrsz = get(0,'ScreenSize');
h = figure('Tag','time_zoom','Name','OFDM magnitude spectrum','Position',[scrsz(3)/2-500 scrsz(4)/2-250 1000 500]);
newaxes = copyobj(handles.ofdmWaveAxes3, h);
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
