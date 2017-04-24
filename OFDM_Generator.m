clear all;

gpib_addr = 1;

% Open serial communication to device.
%g = visa('AGILENT', ['GPIB0::', num2str(gpib_addr), '::INSTR']);
g = visa('AGILENT', ['GPIB0::1::INSTR']);
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

%communication with device is treated like writing into file
try
    disp('Getting device ID.');
    fopen(g);
    % test IDN, set OUT ON
    fprintf(g, '%s\n', '*IDN?');
    id = fscanf(g, '%s');
    fclose(g);
    disp(['Device ID: ' id]);
catch
    text2 = get(gcf,'UserData');
    disp(['Got exception: ' text2]);
end

%data generation
%f = 10e3;
%sampletime = 2.5e-5;

    %precision of f depends how many times is oversampled 20 looks OK.
   f = 7.42e3;
   sampletime = 1/(20 * f);

path = 'c:\LeCroy_LW410';
command = [path, '\LWConv.exe ', path, '\temp.asc ', path, '\temp.dif 0 USERSEQ USERSEQ ', num2str(sampletime)];

if ~exist([path '\LWConv.exe'],'file')
    path = uigetdir('c:\Users\Lubomir Jagos\Documents\VUT_Brno_skola\Diplomka 2\LW410',text2{lang,7});
end
if exist([path '\LWConv.exe'],'file')
        
    t = sampletime:sampletime:(800*sampletime);

    data = sin(2*pi*f*t);
    %data = sawtooth(2*pi*10*t);    %10kHz
    %data = [linspace(0,2,300) linspace(2,0.3,200) linspace(1,0.3,300)];
    
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
end  