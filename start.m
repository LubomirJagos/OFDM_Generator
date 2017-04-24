function varargout = start(varargin)
% START M-file for start.fig
%      START, by itself, creates a new START or raises the existing
%      singleton*.
%
%      H = START returns the handle to a new START or the handle to
%      the existing singleton*.
%
%      START('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in START.M with the given input arguments.
%
%      START('Property','Value',...) creates a new START or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before start_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to start_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help start

% Last Modified by GUIDE v2.5 03-Nov-2016 19:30:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @start_OpeningFcn, ...
                   'gui_OutputFcn',  @start_OutputFcn, ...
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


% --- Executes just before start is made visible.
function start_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to start (see VARARGIN)

% Choose default command line output for start
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
text2 = [{'Level [V]','device not found','Error','Couldn''t connect to the device.','Number of samples is outside the interval.','Work directory not wound.','Select folder containing program LWConv.exe.','Sequence cannot be used.','Wrong set of values.'};
{'Úroveò [V]','zaøízení nenalezeno','Chyba','Nepodaøilo se spojit se zaøízením.','Poèet vzorkù je mimo povolený rozsah.','Nenalezena pracovní složka.','Vyberte složku obsahující program LWConv.exe.','Sekvence nelze použít.','Neplatné zadání.'}];
set(gcf,'UserData',text2);
scrsz = get(0,'ScreenSize');
pos = get(gcf,'Position');
set(gcf,'Position',[(scrsz(3)-pos(3))/2 (scrsz(4)-pos(4))/2 pos(3:4)]);

% UIWAIT makes start wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = start_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in ppm_settings_lang.
function ppm_settings_lang_Callback(hObject, eventdata, handles)
% hObject    handle to ppm_settings_lang (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ppm_settings_lang contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppm_settings_lang
lang = get(hObject,'Value');
text = [{'Line codes','Sequence','Source','random sequence of length','optional assignment','Apply','Repair','Settings','Language :',{'english';'czech'},'GPIB address :','Connection test :','Line code','valid sequence','Bit period [ms] :','Frequency [kHz] :','Samples/bit :','Level L [V] :','Level H [V] :','Line code :','not avaliable','Time domain','Frequency domain','Filter','Filter selection','none','Output','Info','Line code [S/s] :','After filtration [S/s] :','Output [S/s] :','Samples [S] :','Apply'};
{'Linkové kódy','Sekvence','Zdroj','náhodná sekvence délky','vlastní zadání','Potvrdit','Opravit','Nastavení','Jazyk :',{'anglicky';'èesky'},'Adresa GPIB :','Test spojení :','Linkový kód','platná sekvence','Bitová perioda [ms] :','Kmitoèet [kHz] :','Vzorkù/bit :','Úroveò L [V] :','Úroveò H [V] :','Linkový kód :','není dostupné','Èasová oblast','Kmitoètová oblast','Filtr','Výbìr filtru','žádný','Výstup','Informace','Linkový kód [S/s] :','Po filtraci [S/s] :','Výstup [S/s] :','Vzorkù [S] :','Nastavit'}];
set(gcf,'Name',text{lang,1});
set(handles.pnl_seq,'Title',text{lang,2});
set(handles.btn_seq_source,'Title',text{lang,3});
set(handles.seq_source_random,'String',text{lang,4});
set(handles.seq_source_own,'String',text{lang,5});
set(handles.seq_source_valid,'String',text{lang,14});
set(handles.btn_seq_accept,'String',text{lang,6});
set(handles.btn_seq_repair,'String',text{lang,7});
% 
set(handles.pnl_settings,'Title',text{lang,8});
set(handles.txt_settings_lang,'String',text{lang,9});
set(handles.ppm_settings_lang,'String',text{lang,10});
set(handles.txt_settings_gpibaddr,'String',text{lang,11});
set(handles.txt_settings_conntest,'String',text{lang,12});
% 
set(handles.pnl_link,'Title',text{lang,13});
set(handles.pnl_link_settings,'Title',text{lang,8});
set(handles.txt_link_settings_bitper,'String',text{lang,15});
set(handles.txt_link_settings_bitfreq,'String',text{lang,16});
set(handles.txt_link_settings_spb,'String',text{lang,17});
set(handles.txt_link_settings_levelL,'String',text{lang,18});
set(handles.txt_link_settings_levelH,'String',text{lang,19});
set(handles.txt_link_settings_code,'String',text{lang,20});
set(handles.pnl_link_time,'Title',text{lang,22});
%
set(handles.axes_link_na,'String',text{lang,21});
set(handles.axes_filter_time_na,'String',text{lang,21});
set(handles.axes_filter_freq_na,'String',text{lang,21});
set(handles.axes_output_time_na,'String',text{lang,21});
set(handles.axes_output_freq_na,'String',text{lang,21});
set(handles.pnl_filter_settings_na,'String',text{lang,21});
%
set(handles.pnl_filter_time,'Title',text{lang,22});
set(handles.pnl_filter_freq,'Title',text{lang,23});
set(handles.pnl_output_time,'Title',text{lang,22});
set(handles.pnl_output_freq,'Title',text{lang,23});
set(handles.pnl_filter,'Title',text{lang,24});
set(handles.pnl_filter_select,'Title',text{lang,25});
set(handles.filter_none,'String',text{lang,26});
set(handles.pnl_filter_settings,'Title',text{lang,8});
%
set(handles.pnl_output,'Title',text{lang,27});
set(handles.pnl_output_info,'Title',text{lang,28});
set(handles.txt_output_info_lc,'String',text{lang,29});
set(handles.txt_output_info_af,'String',text{lang,30});
set(handles.txt_output_info_out,'String',text{lang,31});
set(handles.txt_output_info_smp,'String',text{lang,32});
set(handles.btn_output_enable,'String',text{lang,33});

text2 = get(gcf,'UserData');
temp = setdiff(1:size(text,1),lang);

%comment by LuboJ.
%if sum(arrayfun(@(t) strcmp(get(handles.edt_settings_test,'String'),text2{t,2}), temp))>0
%    set(handles.edt_settings_test,'String',text2{lang,2});
%end

start('display_link',gcbo,[],guidata(gcbo));

% --- Executes during object creation, after setting all properties.
function ppm_settings_lang_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppm_settings_lang (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edt_seq_disp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_seq_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_seq_accept.
function btn_seq_accept_Callback(hObject, eventdata, handles)
% hObject    handle to btn_seq_accept (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
max_seqlen = 300;
% display_link = true;
if get(handles.seq_source_random,'Value') == 1
    % ***** Random seqence *********************************************
    seqlen = str2double(get(handles.edt_seq_seqlen,'String'));
    if seqlen > max_seqlen
        seqlen = max_seqlen;
        set(handles.edt_seq_seqlen,'String',num2str(max_seqlen));
    end    
    if (seqlen < 1) || isnan(seqlen)
        seqlen = 0;
    end
    
    seq = rand([1 seqlen])<0.5;
    set(handles.edt_seq_disp,'String',strrep(int2str(seq),' ',''));
elseif get(handles.seq_source_own,'Value') == 1
    % ***** User defined sequence **************************************
    temp = get(handles.edt_seq_disp,'String');
    temp = temp(regexp(temp,'[01]'));
    if numel(temp) > max_seqlen
        temp = temp(1:max_seqlen);
    end
    set(handles.edt_seq_disp,'String',temp);
end
start('display_link',gcbo,[],guidata(gcbo))

function edt_seq_seqlen_Callback(hObject, eventdata, handles)
% hObject    handle to edt_seq_seqlen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_seq_seqlen as text
%        str2double(get(hObject,'String')) returns contents of edt_seq_seqlen as a double


% --- Executes during object creation, after setting all properties.
function edt_seq_seqlen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_seq_seqlen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ppm_settings_address.
function ppm_settings_address_Callback(hObject, eventdata, handles)
% hObject    handle to ppm_settings_address (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ppm_settings_address contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppm_settings_address


% --- Executes during object creation, after setting all properties.
function ppm_settings_address_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppm_settings_address (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
if license('test','instr_control_toolbox')
    set(hObject,'Enable','on');
end

% --- Executes on selection change in ppm_link_code.
function ppm_link_code_Callback(hObject, eventdata, handles)
% hObject    handle to ppm_link_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ppm_link_code contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppm_link_code
start('display_link',gcbo,[],guidata(gcbo))

% --- Executes during object creation, after setting all properties.
function ppm_link_code_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppm_link_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_link_bitper_Callback(hObject, eventdata, handles)
% hObject    handle to edt_link_bitper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_link_bitper as text
%        str2double(get(hObject,'String')) returns contents of edt_link_bitper as a double
temp = str2double(get(hObject,'String'));
set(handles.edt_link_freq,'String',num2str(1/temp));
start('display_link',gcbo,[],guidata(gcbo));

% --- Executes during object creation, after setting all properties.
function edt_link_bitper_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_link_bitper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_link_level_L_Callback(hObject, eventdata, handles)
% hObject    handle to edt_link_level_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_link_level_L as text
%        str2double(get(hObject,'String')) returns contents of edt_link_level_L as a double
start('display_link',gcbo,[],guidata(gcbo))

% --- Executes during object creation, after setting all properties.
function edt_link_level_L_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_link_level_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_settings_test.
function btn_settings_test_Callback(hObject, eventdata, handles)
% hObject    handle to btn_settings_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gpib_addr = get(handles.ppm_settings_address,'Value');
g = visa('AGILENT', ['GPIB0::', num2str(gpib_addr), '::INSTR']);
%g = gpib('AGILENT', 0, 1);

% Set the property values.
set(g, 'BoardIndex', 0);
set(g, 'ByteOrder', 'littleEndian');
set(g, 'BytesAvailableFcn', '');
set(g, 'BytesAvailableFcnCount', 48);
set(g, 'BytesAvailableFcnMode', 'eosCharCode');
set(g, 'EOIMode', 'on');
set(g, 'EOSCharCode', 'LF');
set(g, 'EOSMode', 'none');
set(g, 'ErrorFcn', '');
set(g, 'InputBufferSize', 1.2e6);
set(g, 'Name', ['VISA-GPIB0-', num2str(gpib_addr)]);
set(g, 'OutputBufferSize', 1e6);
set(g, 'OutputEmptyFcn', '');
set(g, 'PrimaryAddress', gpib_addr);
set(g, 'RecordDetail', 'compact');
set(g, 'RecordMode', 'overwrite');
set(g, 'RecordName', 'record.txt');
set(g, 'SecondaryAddress', 0);
set(g, 'Tag', '');
set(g, 'Timeout', 50);
set(g, 'TimerFcn', '');
set(g, 'TimerPeriod', 2);
set(g, 'UserData', []);
try
    fopen(g);
    % test IDN, set OUT ON
    fprintf(g, '%s\n', '*IDN?');
    id = fscanf(g, '%s');
    fclose(g);
    set(handles.edt_settings_test,'String',id,'FontWeight','bold','ForegroundColor',[0 0 0]);
    pause(0.8);
    set(handles.edt_settings_test,'FontWeight','normal','ForegroundColor',[0.38 0.38 0.38]);
catch
    text2 = get(gcf,'UserData');
    lang = get(handles.ppm_settings_lang,'Value');
    set(handles.edt_settings_test,'String',text2{lang,2},'FontWeight','bold','ForegroundColor',[1 0 0]);
    pause(0.8);
    set(handles.edt_settings_test,'FontWeight','normal','ForegroundColor',[0.38 0.38 0.38]);
end

function edt_settings_test_Callback(hObject, eventdata, handles)
% hObject    handle to edt_settings_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_settings_test as text
%        str2double(get(hObject,'String')) returns contents of edt_settings_test as a double


% --- Executes during object creation, after setting all properties.
function edt_settings_test_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_settings_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_output_enable.
function btn_output_enable_Callback(hObject, eventdata, handles)
% hObject    handle to btn_output_enable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
text2 = get(gcf,'UserData');
lang = get(handles.ppm_settings_lang,'Value');
temp = get(handles.axes_output_time,'UserData');
sampletime = temp(1,2);
data = temp(2,:);
if (size(data,2) >= 64) && (size(data,2) <= 1e6)
    try
        gpib_addr = get(handles.ppm_settings_address,'Value');

        % save asc file
        path = get(hObject,'UserData');
        if ~exist([path '\LWConv.exe'],'file')
            path = uigetdir('c:\Users\Lubomir Jagos\Documents\VUT_Brno_skola\Diplomka 2\LW410',text2{lang,7});
            set(hObject,'UserData',path);
        end
        if exist([path '\LWConv.exe'],'file')
            fid = fopen([path '\temp.asc'],'w');
            fprintf(fid, '%f\n', data);
            fclose(fid);

            % konverze *.asc -> *.dif a vytvoreni souboru *.dif
            command = [path, '\LWConv.exe ', path, '\temp.asc ', path, '\temp.dif 0 USERSEQ USERSEQ ', num2str(sampletime)];
            dos(command);

            % otevreni souboru *.dif
            fid = fopen([path, '\temp.dif'],'r');
            z = fread(fid, [1 inf], 'uchar');   
            fclose(fid);

            T = strcat('GPIB0::', num2str(gpib_addr), '::INSTR');
            % Create the instrument object.
            g = visa('AGILENT', T);
            %g = gpib('AGILENT', 0, 1);

            % Set the property values.
            set(g, 'BoardIndex', 0);
            set(g, 'ByteOrder', 'littleEndian');
            set(g, 'BytesAvailableFcn', '');
            set(g, 'BytesAvailableFcnCount', 48);
            set(g, 'BytesAvailableFcnMode', 'eosCharCode');
            set(g, 'EOIMode', 'on');
            set(g, 'EOSCharCode', 'LF');
            set(g, 'EOSMode', 'none');
            set(g, 'ErrorFcn', '');
            set(g, 'InputBufferSize', 1.2e6);
            set(g, 'Name', ['VISA-GPIB0-', num2str(gpib_addr)]);
            set(g, 'OutputBufferSize', 1e6);
            set(g, 'OutputEmptyFcn', '');
            set(g, 'PrimaryAddress', gpib_addr);
            set(g, 'RecordDetail', 'compact');
            set(g, 'RecordMode', 'overwrite');
            set(g, 'RecordName', 'record.txt');
            set(g, 'SecondaryAddress', 0);
            set(g, 'Tag', '');
            set(g, 'Timeout', 50);
            set(g, 'TimerFcn', '');
            set(g, 'TimerPeriod', 2);
            set(g, 'UserData', []);
            fopen(g);
            % test IDN, set OUT ON
            fprintf(g, '%s\n', '*IDN?');
            fscanf(g, '%s');
            fprintf(g, '%s\n', 'OUTPut1 on');
            % write data
            data = (z);
            fwrite(g, ['WAVE:DATA ', data], 'char')
            fclose(g);
        else
            msgbox(text2{lang,6},text2{lang,3},'error');
        end
    catch
        msgbox(text2{lang,4},text2{lang,3},'error');
    end
else
    msgbox(text2{lang,5},text2{lang,3},'error');    
end
% --- Executes when selected object is changed in btn_seq_source.
function btn_seq_source_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in btn_seq_source 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(hObject,'Tag')
    case 'seq_source_random'
        set(handles.edt_seq_disp,'Enable','off','String','');
        set(handles.edt_seq_seqlen,'Enable','on','String','20');
    case 'seq_source_own'
        set(handles.edt_seq_seqlen,'Enable','off','String','');
        set(handles.edt_seq_disp,'Enable','on','String','');
    case 'seq_source_valid'
        set(handles.edt_seq_seqlen,'Enable','off','String','');
        set(handles.edt_seq_disp,'Enable','off','String','1011110000100001110000110100001010');
end

function display_link(hObject, eventdata, handles)
text2 = get(gcf,'UserData');
lang = get(handles.ppm_settings_lang,'Value');
seq = str2num(regexprep(get(handles.edt_seq_disp,'String'), '([01])', '$1 '));
seqlen = numel(seq);
T = str2double(get(handles.edt_link_bitper,'String'));
L = str2double(get(handles.edt_link_level_L,'String'));
if L < -5
    L = -5;
    set(handles.edt_link_level_L,'String','-5');
end
H = str2double(get(handles.edt_link_level_H,'String'));
if H > 5
    H = 5;
    set(handles.edt_link_level_H,'String','5');
end
link_id = get(handles.ppm_link_code,'Value');
if (seqlen > 0) && (T > 0) && (H > L) && isnumeric(L) && isnumeric(H) && isnumeric(T)
    T = str2double(get(handles.edt_link_bitper,'String'))/1000;
    if T > 0.1
        T = 0.1;
        set(handles.edt_link_bitper,'String',num2str(T*1000));
        set(handles.edt_link_freq,'String',num2str(1/T/1000));
    end
    t = 0:T:seqlen*T;
    tx = t(1:2:end);
    txx = 0:T/2:seqlen*T;
    axes(handles.axes_link);
    all_ok = true;
    repair = true;
    O = abs(H-L)/2+L;
    A = abs(H-L)/2;
    spbval = str2double(get(handles.ppm_spb,'String'));
    spbindex = get(handles.ppm_spb,'Value');
    mx = spbval(spbindex); % samples per period, minimum = 2
    link_names = get(handles.ppm_link_code,'String');
    switch link_names{link_id}
        case '2B1Q'
            set(handles.edt_link_level_L,'Enable','on');
            set(handles.edt_link_level_H,'Enable','on');
            if ~mod(seqlen,2)
                temp = reshape(seq,2,seqlen/2);
                C = (3-2*temp(2,:)).*(2*temp(1,:)-1);
                C = C*A/3+O;
                stairs(tx,[C C(end)]);
                C = reshape(repmat(C,2*mx,1),1,2*mx*numel(C));
            else
                repair = false;
                all_ok = false;
            end
        case 'AMI'
            set(handles.edt_link_level_L,'Enable','on');
            set(handles.edt_link_level_H,'Enable','on');
            if ~mod(sum(seq),2)
                C = cumprod(1-2*seq).*seq;
                C = -C*A+O;
                stairs(t,[C C(end)]);
                C = reshape(repmat(C,mx,1),1,mx*numel(C));
            else
                all_ok = false;
            end
        case 'Bi-Phi-L'
            set(handles.edt_link_level_L,'Enable','on');
            set(handles.edt_link_level_H,'Enable','on');
            temp = reshape(repmat(seq,2,1),1,2*numel(seq));
            C = repmat([-1 1],1,seqlen)+repmat([2 -2],1,seqlen).*temp;
            C = C*A+O;
            stairs(txx,[C C(end)]);
            C = reshape(repmat(C,mx/2,1),1,mx/2*numel(C));
            
        case 'Bi-Phi-M'
            set(handles.edt_link_level_L,'Enable','on');
            set(handles.edt_link_level_H,'Enable','on');
            if ~mod(sum(~seq),2)
                C = zeros(1,2*length(seq));
                C(1) = 0;
                C(2) = seq(1);
                for k = 2:length(seq)
                    if seq(k) == 1
                        C(2*k-1:2*k) = [1-C(2*k-2) C(2*k-2)];
                    else
                        C(2*k-1:2*k) = [1-C(2*k-2) 1-C(2*k-2)];
                    end
                end
                C = (2*C-1)*A+O;
                stairs(txx,[C C(end)]);
                C = reshape(repmat(C,mx/2,1),1,mx/2*numel(C));          
            else
                all_ok = false;
            end
            
        case 'Bi-Phi-S'
            set(handles.edt_link_level_L,'Enable','on');
            set(handles.edt_link_level_H,'Enable','on');
            if ~mod(sum(seq),2)            
                C = zeros(1,2*length(seq));
                C(1) = 0;
                C(2) = 1-seq(1);
                for k = 2:length(seq)
                    if seq(k) == 0
                        C(2*k-1:2*k) = [1-C(2*k-2) C(2*k-2)];
                    else
                        C(2*k-1:2*k) = [1-C(2*k-2) 1-C(2*k-2)];
                    end
                end
                C = (2*C-1)*A+O;
                stairs(txx,[C C(end)]);
                C = reshape(repmat(C,mx/2,1),1,mx/2*numel(C));          
            else
                all_ok = false;
            end
            
        case 'Bipolar RZ'
            set(handles.edt_link_level_L,'Enable','on');
            set(handles.edt_link_level_H,'Enable','on');
            C = zeros(1,2*seqlen);
            C(1:2:2*seqlen) = 2*seq-1;
            C = C*A+O;
            stairs(txx,[C C(end)]);
            C = reshape(repmat(C,mx/2,1),1,mx/2*numel(C));
            
        case 'Delay modulation'
            set(handles.edt_link_level_L,'Enable','on');
            set(handles.edt_link_level_H,'Enable','on');
            C = zeros(1,6*length(seq));
            C(1) = seq(1);
            temp = [seq seq seq];
            for k = 2:length(temp)
                if temp(k) == 1
                    C(2*k-1:2*k) = [C(2*k-2) 1-C(2*k-2)];
                else
                    if temp(k-1) == 0
                        C(2*k-1:2*k) = [1-C(2*k-2) 1-C(2*k-2)];
                    else
                        C(2*k-1:2*k) = [C(2*k-2) C(2*k-2)];                    
                    end
                end
            end
            if all(C(2*seqlen+1:2*seqlen+2) == C(4*seqlen+1:4*seqlen+2))
                C = C(2*seqlen+1:4*seqlen);
                C = (2*C-1)*A+O;
                stairs(txx,[C C(end)]);
                C = reshape(repmat(C,mx/2,1),1,mx/2*numel(C));          
            else
                all_ok = false;
                repair = false;
            end
            
        case 'Dicode NRZ'
            set(handles.edt_link_level_L,'Enable','on');
            set(handles.edt_link_level_H,'Enable','on');
            C = (1-2*seq).*abs(diff([seq(end) seq]));
            C = C*A+O;
            stairs(t,[C C(end)]);
            C = reshape(repmat(C,mx,1),1,mx*numel(C));
            
        case 'Dicode RZ'
            set(handles.edt_link_level_L,'Enable','on');
            set(handles.edt_link_level_H,'Enable','on');
            C = zeros(1,2*seqlen);
            C(1:2:2*seqlen) = (1-2*seq).*abs(diff([seq(end) seq]));
            C = C*A+O;
            stairs(txx,[C C(end)]);
            C = reshape(repmat(C,mx/2,1),1,mx/2*numel(C));            
            
        case 'HDB3'
            set(handles.edt_link_level_L,'Enable','on');
            set(handles.edt_link_level_H,'Enable','on');
            if ~mod(sum(seq),2)
                C = cumprod(1-2*seq).*seq;
                C = [C, C, C];
                substs = [];
                k = 1;
                while(k<=length(C)-3)
                    if sum(abs(C(k:k+3))) == 0
                        substs = [substs k];
                        k = k + 4;
                    else
                        k = k + 1;
                    end
                end
                if ~mod(sum(substs>seqlen & substs<=2*seqlen),2)
                    for k = 2:length(substs)
                        if mod(sum(abs(C(substs(k-1)+4:substs(k)-1))),2)
                            % odd ones
                            C(substs(k)+3) = C(substs(k)-1);
                        else
                            % even ones
                            C(substs(k)) = -C(substs(k)-1);
                            C(substs(k)+3) = -C(substs(k)-1);
                            C(substs(k)+4:end) = -C(substs(k)+4:end);
                        end
                    end
                    C = C(seqlen+1:2*seqlen);
                    C = C*A+O;
                    stairs(t,[C C(end)]);
                    C = reshape(repmat(C,mx,1),1,mx*numel(C));
                else
                    all_ok = false;
                    repair = false;
                end
            else
                all_ok = false;
                repair = false;
            end   
            
        case 'NRZ-L'
            set(handles.edt_link_level_L,'Enable','on');
            set(handles.edt_link_level_H,'Enable','on');
            C = 2*seq-1;
            C = C*A+O;
            stairs(t,[C C(end)]);
            C = reshape(repmat(C,mx,1),1,mx*numel(C));
            
        case 'NRZ-M'
            set(handles.edt_link_level_L,'Enable','on');
            set(handles.edt_link_level_H,'Enable','on');
            if ~mod(sum(seq),2)
                C = 2*mod(cumsum(seq),2)-1;
                C = C*A+O;
                stairs(t,[C C(end)]);
                C = reshape(repmat(C,mx,1),1,mx*numel(C));
            else
                all_ok = false;
            end
            
        case 'NRZ-S'
            set(handles.edt_link_level_L,'Enable','on');
            set(handles.edt_link_level_H,'Enable','on');
            if ~mod(sum(~seq),2)
                C = 2*mod(cumsum(abs(seq-1)),2)-1;
                C = C*A+O;
                stairs(t,[C C(end)]);
                C = reshape(repmat(C,mx,1),1,mx*numel(C));
            else
                all_ok = false;
            end
            
        case 'RZ-AMI'
            set(handles.edt_link_level_L,'Enable','on');
            set(handles.edt_link_level_H,'Enable','on');
            if ~mod(sum(seq),2)
                C = zeros(1,2*seqlen);
                C(1:2:2*seqlen) = cumprod(1-2*seq).*seq;
                C = -C*A+O;
                stairs(txx,[C C(end)]);
                C = reshape(repmat(C,mx/2,1),1,mx/2*numel(C));
            else
                all_ok = false;
            end

        case 'RZ-HDB3'
            set(handles.edt_link_level_L,'Enable','on');
            set(handles.edt_link_level_H,'Enable','on');
            if ~mod(sum(seq),2)
                C = cumprod(1-2*seq).*seq;
                C = [C, C, C];
                substs = [];
                k = 1;
                while(k<=length(C)-3)
                    if sum(abs(C(k:k+3))) == 0
                        substs = [substs k];
                        k = k + 4;
                    else
                        k = k + 1;
                    end
                end
                if ~mod(sum(substs>seqlen & substs<=2*seqlen),2)
                    for k = 2:length(substs)
                        if mod(sum(abs(C(substs(k-1)+4:substs(k)-1))),2)
                            % odd ones
                            C(substs(k)+3) = C(substs(k)-1);
                        else
                            % even ones
                            C(substs(k)) = -C(substs(k)-1);
                            C(substs(k)+3) = -C(substs(k)-1);
                            C(substs(k)+4:end) = -C(substs(k)+4:end);
                        end
                    end
                    C = C(seqlen+1:2*seqlen);
                    C = C*A+O;
                    C(2,:) = zeros(1,size(C,2));
                    C = reshape(C,1,numel(C));
                    stairs(txx,[C C(end)]);
                    C = reshape(repmat(C,mx/2,1),1,mx/2*numel(C));
                else
                    all_ok = false;
                    repair = false;
                end
            else
                all_ok = false;
                repair = false;
            end   
            
        case 'Unipolar RZ'
            set(handles.edt_link_level_L,'Enable','off');
            set(handles.edt_link_level_H,'Enable','on');
            C = zeros(1,2*seqlen);
            C(1:2:2*seqlen) = seq;
            C = C*H;
            L = 0;
            O = H/2;
            stairs(txx,[C C(end)]);
            C = reshape(repmat(C,mx/2,1),1,mx/2*numel(C)); 
            
        otherwise
            all_ok = false;
            repair = false;
    end
    if all_ok
        xlim([0 t(end)]);
        if strcmp(get(handles.edt_link_level_L,'Enable'),'on') || strcmp(get(handles.edt_link_level_H,'Enable'),'on')
            if max(abs(C)) == 0
                ylim([-1 1]);
            else
                ylim([L H]*1.2+[-0.5 0.5]);
            end
        end
        seqn = hggroup;
        for k = 1:seqlen
            text(t(k)+T/2,O,num2str(seq(k)),'HorizontalAlignment','center','VerticalAlignment','middle','Parent',seqn,'FontSize',4+round(200/(seqlen+1)),'Color',[0.85 0.85 0.85],'FontWeight','bold');
        end
        set(handles.axes_link,'Children',flipud(get(handles.axes_link,'Children')));
        set(handles.axes_link,'TickLength',[0.003 0.1],'XTick',linspace(0,t(end),21))
        xlabel('t [s]');
        text2 = get(gcf,'UserData');
        lang = get(handles.ppm_settings_lang,'Value');
        ylabel(text2{lang,1});
        set(handles.btn_seq_repair,'Enable','off');
        set(handles.axes_link_na,'Visible','off');
        
        t = linspace(0,seqlen*T-T/mx,numel(C));
        set(handles.edt_output_info_fvz_link,'String',strrep(num2str(1/t(2),'%0.2e'), 'e+0', 'e+'));
        set(handles.axes_link,'UserData',[t; C]);
        set(handles.filter_none,'Enable','on');
        if license('test','communication_toolbox') &&  license('test','signal_toolbox')
            set(handles.filter_rc,'Enable','on');
            set(handles.filter_srrc,'Enable','on');
        end
        if license('test','instr_control_toolbox')
            set(handles.btn_output_enable,'Enable','on');
            set(handles.btn_settings_test,'Enable','on');
        end
        temp = get(handles.filter_none,'Value')*handles.filter_none+get(handles.filter_rc,'Value')*handles.filter_rc+get(handles.filter_srrc,'Value')*handles.filter_srrc;
        pnl_filter_select_SelectionChangeFcn(temp, eventdata, handles);
    else
        if strcmp(get(handles.edt_seq_seqlen,'Enable'),'on') && repair
            start('btn_seq_repair_Callback',handles.btn_seq_repair,[],guidata(gcbo));   % repair random sequence
        else
            % axes
            set(handles.axes_link,'UserData',[]);
            delete([get(handles.axes_link,'Children'); get(handles.axes_filter_time,'Children'); get(handles.axes_filter_freq,'Children'); get(handles.axes_output_time,'Children'); get(handles.axes_output_freq,'Children')]);
            set([handles.axes_link handles.axes_filter_time handles.axes_filter_freq handles.axes_output_time handles.btn_output_time_zoom handles.axes_output_freq handles.btn_output_freq_zoom],'Visible','off');
            set([handles.axes_link_na handles.axes_filter_time_na handles.axes_filter_freq_na handles.axes_output_time_na handles.axes_output_freq_na],'Visible','on');

            % filter
            set(get(handles.pnl_filter_select,'Children'),'Enable','off');
            set(get(handles.pnl_filter_settings,'Children'),'Enable','off');
            
            % output
            set(handles.btn_output_enable,'Enable','off');

            if repair
                set(handles.btn_seq_repair,'Enable','on');                                  % let user to repair sequence manually
            end
            msgbox(text2{lang,8},text2{lang,3},'error');
        end
    end
else
    % axes
    set(handles.axes_link,'UserData',[]);
    delete([get(handles.axes_link,'Children'); get(handles.axes_filter_time,'Children'); get(handles.axes_filter_freq,'Children'); get(handles.axes_output_time,'Children'); get(handles.axes_output_freq,'Children')]);
    set([handles.axes_link handles.axes_filter_time handles.axes_filter_freq handles.axes_output_time handles.btn_output_time_zoom handles.axes_output_freq handles.btn_output_freq_zoom],'Visible','off');
    set([handles.axes_link_na handles.axes_filter_time_na handles.axes_filter_freq_na handles.axes_output_time_na handles.axes_output_freq_na],'Visible','on');
    % filter
    set(get(handles.pnl_filter_select,'Children'),'Enable','off');
    set(get(handles.pnl_filter_settings,'Children'),'Enable','off');
    % output
    set(handles.btn_output_enable,'Enable','off');
    msgbox(text2{lang,9},text2{lang,3},'error');
end



function edt_link_freq_Callback(hObject, eventdata, handles)
% hObject    handle to edt_link_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_link_freq as text
%        str2double(get(hObject,'String')) returns contents of edt_link_freq as a double
temp = str2double(get(hObject,'String'));
set(handles.edt_link_bitper,'String',num2str(1/temp));
start('display_link',gcbo,[],guidata(gcbo));

% --- Executes during object creation, after setting all properties.
function edt_link_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_link_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_link_level_H_Callback(hObject, eventdata, handles)
% hObject    handle to edt_link_level_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_link_level_H as text
%        str2double(get(hObject,'String')) returns contents of edt_link_level_H as a double
start('display_link',gcbo,[],guidata(gcbo));

% --- Executes during object creation, after setting all properties.
function edt_link_level_H_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_link_level_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_seq_repair.
function btn_seq_repair_Callback(hObject, eventdata, handles)
% hObject    handle to btn_seq_repair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
max_seqlen = 300;
temp = get(handles.edt_seq_disp,'String');
temp = temp(regexp(temp,'[01]'));
if numel(temp) > max_seqlen
    temp = temp(1:max_seqlen);
end
seq = str2num(regexprep(temp, '([01])', '$1 '));
% link_names = get(handles.ppm_link_code,'String');
% if link_names{link_id} == 'AMI' || link_names{link_id} == 'RZ-AMI' || link_names{link_id} == 'NRZ-M' || link_names{link_id} == 'NRZ-S' || link_names{link_id} == 'NRZ-S'
    seq(end) = ~seq(end);
% end
set(handles.edt_seq_disp,'String',strrep(int2str(seq),' ',''));
start('display_link',gcbo,[],guidata(gcbo));
set(hObject,'Enable','off');


% --- Executes when selected object is changed in pnl_filter_select.
function pnl_filter_select_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in pnl_filter_select 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
set(handles.axes_output_time_na,'Visible','off');
set(handles.axes_output_freq_na,'Visible','off');
seq = get(handles.axes_link,'UserData');
seq_t = seq(1,:);
seq = seq(2,:);
text2 = get(gcf,'UserData');
lang = get(handles.ppm_settings_lang,'Value');
switch get(hObject,'Tag')
    case 'filter_none'
        delete(get(handles.axes_filter_time,'Children'));
        delete(get(handles.axes_filter_freq,'Children'));
        set(handles.axes_filter_time,'Visible','off');
        set(handles.axes_filter_freq,'Visible','off');
        set(handles.axes_filter_time_na,'Visible','on');        
        set(handles.axes_filter_freq_na,'Visible','on');
        set(handles.pnl_filter_settings_na,'Visible','on');
        set([handles.ppm_filter_settings_fsfd handles.ppm_filter_settings_rolloff handles.ppm_filter_settings_delay handles.text_fsfd handles.text_rolloff handles.text_delay],'Visible','off');
        fsfdval = 1;
        fsfdindex = 1;
        hFilter = [];
        fvz = 1/seq_t(2);
        fvz_av = [4e4, 4e5, 4e6, 4e7, 4e8];
        fvz_av = fvz_av(fvz<=fvz_av);
        filtered_seq_t = (0:length(seq)*fvz_av(1)/fvz-1)*1/fvz_av(1);
        %commented by LuboJ.
        %filtered_seq = seq(arrayfun(@(x)sum(round(x*1e9)>=round(seq_t*1e9)), filtered_seq_t));
        set(handles.edt_output_info_fvz_filter,'String',get(handles.edt_output_info_fvz_link,'String'));
    case 'filter_rc'
        set(handles.pnl_filter_settings_na,'Visible','off');
        set(handles.axes_filter_time_na,'Visible','off');        
        set(handles.axes_filter_freq_na,'Visible','off');        
        set([handles.ppm_filter_settings_fsfd handles.ppm_filter_settings_rolloff handles.ppm_filter_settings_delay handles.text_fsfd handles.text_rolloff handles.text_delay],'Visible','on');
        set([handles.ppm_filter_settings_fsfd handles.ppm_filter_settings_rolloff handles.ppm_filter_settings_delay handles.text_fsfd handles.text_rolloff handles.text_delay],'Enable','on');
        set(handles.text_fsfd,'Visible','on');
        set(handles.text_rolloff,'Visible','on');
        set(handles.text_delay,'Visible','on');
        axes(handles.axes_filter_time);
        fsfdval = str2double(get(handles.ppm_filter_settings_fsfd,'String'));
        rolloffval = str2double(get(handles.ppm_filter_settings_rolloff,'String'));
        delayval = str2double(get(handles.ppm_filter_settings_delay,'String'));
        fsfdindex = get(handles.ppm_filter_settings_fsfd,'Value');
        rolloffindex = get(handles.ppm_filter_settings_rolloff,'Value');
        delayindex = get(handles.ppm_filter_settings_delay,'Value');
        hFilter = rcosine(1,fsfdval(fsfdindex),'normal',rolloffval(rolloffindex),delayval(delayindex));
        dt = seq_t(2)/fsfdval(fsfdindex);
        hFilter_t = ((1:length(hFilter))-ceil(length(hFilter)/2))*dt;
        stem(hFilter_t,hFilter,'MarkerSize',3);
        xlim([hFilter_t(1) hFilter_t(end)]);
        ylim([-0.5 1.5]);
        xlabel('t [s]');
        ylabel(text2{lang,1});
        
        axes(handles.axes_filter_freq);
        fFilter = fftshift(abs(fft(hFilter)))/fsfdval(fsfdindex);
        f = linspace(-1/dt/2,1/dt/2,length(fFilter));
        plot(f,fFilter);
        xlim([f(1) f(end)]);
        ylim([0 1.4]);        
        xlabel('f [Hz]');
        ylabel('A [-]');
        filtered_seq = rcosflt(seq,1,fsfdval(fsfdindex),'normal',rolloffval(rolloffindex),delayval(delayindex))';
        temp = floor(length(hFilter)/2);
        repetices = ceil(temp/length(seq))*2+1;
        filtered_seqx = zeros(1,length(seq)*fsfdval(fsfdindex)*repetices+2*temp);
        for r = 1:repetices    % repetices
            filtered_seqx(length(seq)*fsfdval(fsfdindex)*(r-1)+1:length(seq)*fsfdval(fsfdindex)*(r-1)+length(filtered_seq)) = filtered_seqx(length(seq)*fsfdval(fsfdindex)*(r-1)+1:length(seq)*fsfdval(fsfdindex)*(r-1)+length(filtered_seq)) + filtered_seq;
        end
        index1 = floor(repetices/2)*length(seq)*fsfdval(fsfdindex)+temp+1;
        index2 = (floor(repetices/2)+1)*length(seq)*fsfdval(fsfdindex)+temp;
        fvz = 1/seq_t(2)*fsfdval(fsfdindex);
        set(handles.edt_output_info_fvz_filter,'String',strrep(num2str(fvz,'%0.2e'), 'e+0', 'e+'));
        fvz_av = [4e4, 4e5, 4e6, 4e7, 4e8];
        fvz_av = fvz_av(fvz<=fvz_av);
        [n,d] = rat(fvz_av(1)/fvz);
        filtered_seq = resample(filtered_seqx,n,d);
        filtered_seq = filtered_seq(round((index1-1)*n/d)+1:round(index2*n/d));

    case 'filter_srrc'
        set(handles.pnl_filter_settings_na,'Visible','off');
        set(handles.axes_filter_time_na,'Visible','off');        
        set(handles.axes_filter_freq_na,'Visible','off');        
        set([handles.ppm_filter_settings_fsfd handles.ppm_filter_settings_rolloff handles.ppm_filter_settings_delay handles.text_fsfd handles.text_rolloff handles.text_delay],'Visible','on');
        set([handles.ppm_filter_settings_fsfd handles.ppm_filter_settings_rolloff handles.ppm_filter_settings_delay handles.text_fsfd handles.text_rolloff handles.text_delay],'Enable','on');
        set(handles.text_fsfd,'Visible','on');
        set(handles.text_rolloff,'Visible','on');
        set(handles.text_delay,'Visible','on');
        axes(handles.axes_filter_time);
        fsfdval = str2double(get(handles.ppm_filter_settings_fsfd,'String'));
        rolloffval = str2double(get(handles.ppm_filter_settings_rolloff,'String'));
        delayval = str2double(get(handles.ppm_filter_settings_delay,'String'));
        fsfdindex = get(handles.ppm_filter_settings_fsfd,'Value');
        rolloffindex = get(handles.ppm_filter_settings_rolloff,'Value');
        delayindex = get(handles.ppm_filter_settings_delay,'Value');
        hFilter = rcosine(1,fsfdval(fsfdindex),'sqrt',rolloffval(rolloffindex),delayval(delayindex));
        dt = seq_t(2)/fsfdval(fsfdindex);
        hFilter_t = ((1:length(hFilter))-ceil(length(hFilter)/2))*dt;
        stem(hFilter_t,hFilter,'MarkerSize',3);
        xlim([hFilter_t(1) hFilter_t(end)]);
        ylim([-0.25 1]);
        xlabel('t [s]');
        ylabel(text2{lang,1});
        
        axes(handles.axes_filter_freq);
        fFilter = fftshift(abs(fft(hFilter)))/fsfdval(fsfdindex);
        f = linspace(-1/dt/2,1/dt/2,length(fFilter));
        plot(f,fFilter);
        xlim([f(1) f(end)]);
        ylim([0 1.4]);
        xlabel('f [Hz]');
        ylabel('A [-]');        
        filtered_seq = rcosflt(seq,1,fsfdval(fsfdindex),'sqrt',rolloffval(rolloffindex),delayval(delayindex))';
        temp = floor(length(hFilter)/2);
        repetices = ceil(temp/length(seq))*2+1;
        filtered_seqx = zeros(1,length(seq)*fsfdval(fsfdindex)*repetices+2*temp);
        for r = 1:repetices    % repetices
            filtered_seqx(length(seq)*fsfdval(fsfdindex)*(r-1)+1:length(seq)*fsfdval(fsfdindex)*(r-1)+length(filtered_seq)) = filtered_seqx(length(seq)*fsfdval(fsfdindex)*(r-1)+1:length(seq)*fsfdval(fsfdindex)*(r-1)+length(filtered_seq)) + filtered_seq;
        end
        index1 = floor(repetices/2)*length(seq)*fsfdval(fsfdindex)+temp+1;
        index2 = (floor(repetices/2)+1)*length(seq)*fsfdval(fsfdindex)+temp;
        fvz = 1/seq_t(2)*fsfdval(fsfdindex);
        set(handles.edt_output_info_fvz_filter,'String',strrep(num2str(fvz,'%0.2e'), 'e+0', 'e+'));
        fvz_av = [4e4, 4e5, 4e6, 4e7, 4e8];
        fvz_av = fvz_av(fvz<=fvz_av);
        [n,d] = rat(fvz_av(1)/fvz);
        filtered_seq = resample(filtered_seqx,n,d);
        filtered_seq = filtered_seq(round((index1-1)*n/d)+1:round(index2*n/d));
end
set(handles.btn_output_time_zoom,'Visible','on');
axes(handles.axes_output_time);
filtered_seq_t = (0:length(filtered_seq)-1)*1/fvz_av(1);
stem(filtered_seq_t,filtered_seq,'Marker','.','MarkerSize',1,'MarkerFaceColor','b','MarkerEdgeColor','b');
temp = get(gca,'Children');
set(temp(1),'Color',[0.9 0.9 0.9]);
uistack(temp(1),'down');
xlim([0 filtered_seq_t(end)]);
ylim([min([filtered_seq 0]) max([filtered_seq 0])]+[-0.5 0.5]);
xlabel('t [s]');
ylabel(text2{lang,1});
set(handles.axes_output_time,'UserData',[filtered_seq_t;filtered_seq]);
set(handles.edt_output_info_fvz_output,'String',strrep(num2str(1/filtered_seq_t(2),'%0.2e'), 'e+0', 'e+'));
set(handles.edt_output_info_smp,'String',strrep(num2str(length(filtered_seq),'%0.2e'), 'e+0', 'e+'));

set(handles.btn_output_freq_zoom,'Visible','on');
axes(handles.axes_output_freq);
Fs = fvz_av(1);
nfft = 2^(nextpow2(length(filtered_seq))); 
fftx = fft(filtered_seq,nfft); 
NumUniquePts = ceil((nfft+1)/2); 
fftx = fftx(1:NumUniquePts); 
seq_psd = abs(fftx)/length(filtered_seq); 
seq_psd = seq_psd.^2; 
% Since we dropped half the FFT, we multiply mx by 2 to keep the same energy.
% The DC component and Nyquist component, if it exists, are unique and should not be multiplied by 2.
if rem(nfft, 2) % odd nfft excludes Nyquist point
  seq_psd(2:end) = seq_psd(2:end)*2;
else
  seq_psd(2:end -1) = seq_psd(2:end -1)*2;
end
f = (0:NumUniquePts-1)*Fs/nfft; 
plot(f(f<=(1/seq_t(2))/2),seq_psd(f<=(1/seq_t(2))/2));
xlabel('f [Hz]');
ylabel('PSD [W/Hz]');
xlim([0 (1/seq_t(2))/2]);
ylim([0 max(seq_psd)*1.5]);

function apply_filter(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function pnl_filter_select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pnl_filter_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function btn_settings_test_CreateFcn(hObject, eventdata, handles)
% hObject    handle to btn_settings_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if license('test','instr_control_toolbox')
    set(hObject,'Enable','on');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function ppm_filter_settings_rolloff_Callback(hObject, eventdata, handles)
% hObject    handle to ppm_filter_settings_rolloff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ppm_filter_settings_rolloff as text
%        str2double(get(hObject,'String')) returns contents of ppm_filter_settings_rolloff as a double
if get(handles.filter_rc,'Value') == 1
    pnl_filter_select_SelectionChangeFcn(handles.filter_rc, eventdata, handles);
else
    pnl_filter_select_SelectionChangeFcn(handles.filter_srrc, eventdata, handles);
end


% --- Executes during object creation, after setting all properties.
function ppm_filter_settings_rolloff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppm_filter_settings_rolloff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ppm_filter_settings_delay.
function ppm_filter_settings_delay_Callback(hObject, eventdata, handles)
% hObject    handle to ppm_filter_settings_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ppm_filter_settings_delay contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppm_filter_settings_delay
if get(handles.filter_rc,'Value') == 1
    pnl_filter_select_SelectionChangeFcn(handles.filter_rc, eventdata, handles);
else
    pnl_filter_select_SelectionChangeFcn(handles.filter_srrc, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function ppm_filter_settings_delay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppm_filter_settings_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ppm_filter_settings_fsfd.
function ppm_filter_settings_fsfd_Callback(hObject, eventdata, handles)
% hObject    handle to ppm_filter_settings_fsfd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ppm_filter_settings_fsfd contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppm_filter_settings_fsfd
if get(handles.filter_rc,'Value') == 1
    pnl_filter_select_SelectionChangeFcn(handles.filter_rc, eventdata, handles);
else
    pnl_filter_select_SelectionChangeFcn(handles.filter_srrc, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function ppm_filter_settings_fsfd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppm_filter_settings_fsfd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ppm_spb.
function ppm_spb_Callback(hObject, eventdata, handles)
% hObject    handle to ppm_spb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ppm_spb contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ppm_spb
start('display_link',gcbo,[],guidata(gcbo));

% --- Executes during object creation, after setting all properties.
function ppm_spb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppm_spb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_output_info_fvz_link_Callback(hObject, eventdata, handles)
% hObject    handle to edt_output_info_fvz_link (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_output_info_fvz_link as text
%        str2double(get(hObject,'String')) returns contents of edt_output_info_fvz_link as a double


% --- Executes during object creation, after setting all properties.
function edt_output_info_fvz_link_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_output_info_fvz_link (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_output_info_fvz_filter_Callback(hObject, eventdata, handles)
% hObject    handle to edt_output_info_fvz_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_output_info_fvz_filter as text
%        str2double(get(hObject,'String')) returns contents of edt_output_info_fvz_filter as a double


% --- Executes during object creation, after setting all properties.
function edt_output_info_fvz_filter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_output_info_fvz_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_output_info_fvz_output_Callback(hObject, eventdata, handles)
% hObject    handle to edt_output_info_fvz_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_output_info_fvz_output as text
%        str2double(get(hObject,'String')) returns contents of edt_output_info_fvz_output as a double


% --- Executes during object creation, after setting all properties.
function edt_output_info_fvz_output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_output_info_fvz_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_output_info_smp_Callback(hObject, eventdata, handles)
% hObject    handle to edt_output_info_smp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_output_info_smp as text
%        str2double(get(hObject,'String')) returns contents of edt_output_info_smp as a double


% --- Executes during object creation, after setting all properties.
function edt_output_info_smp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_output_info_smp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_output_time_zoom.
function btn_output_time_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to btn_output_time_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
scrsz = get(0,'ScreenSize');
h = figure('Tag','time_zoom','Name',get(handles.pnl_output_time,'Title'),'Position',[scrsz(3)/2-500 scrsz(4)/2-250 1000 500]);
newaxes = copyobj(handles.axes_output_time,h);
set(newaxes,'Units','normalized','ActivePositionProperty','position','Position',[0.1 0.1 0.85 0.85]);

% --- Executes on button press in btn_output_freq_zoom.
function btn_output_freq_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to btn_output_freq_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
scrsz = get(0,'ScreenSize');
h = figure('Tag','freq_zoom','Name',get(handles.pnl_output_freq,'Title'),'Position',[scrsz(3)/2-500 scrsz(4)/2-250 1000 500]);
newaxes = copyobj(handles.axes_output_freq,h);
set(newaxes,'Units','normalized','ActivePositionProperty','position','Position',[0.1 0.1 0.85 0.85]);
