clear all;

f = 1e3;
maxSamples = 260.000e3;

clear obj.g
obj.g = visa('AGILENT', ['GPIB0::1::INSTR'])
set(obj.g, 'BoardIndex', 0);
set(obj.g, 'ByteOrder', 'littleEndian');
set(obj.g, 'BytesAvailableFcn', '');
set(obj.g, 'BytesAvailableFcnCount', 48);
set(obj.g, 'BytesAvailableFcnMode', 'eosCharCode');
set(obj.g, 'EOIMode', 'on');
set(obj.g, 'EOSCharCode', 'LF');
set(obj.g, 'EOSMode', 'none');
set(obj.g, 'ErrorFcn', '');
set(obj.g, 'InputBufferSize', 1.2e6);
set(obj.g, 'Name', ['VISA-GPIB0-', '1']);
set(obj.g, 'OutputBufferSize', 10e6);
set(obj.g, 'OutputEmptyFcn', '');
set(obj.g, 'PrimaryAddress', 1);
set(obj.g, 'RecordDetail', 'compact');
set(obj.g, 'RecordMode', 'overwrite');
set(obj.g, 'RecordName', 'record.txt');
set(obj.g, 'SecondaryAddress', 0);
set(obj.g, 'Tag', '');
set(obj.g, 'Timeout', 50);
set(obj.g, 'TimerFcn', '');
set(obj.g, 'TimerPeriod', 2);
set(obj.g, 'UserData', []);

obj.g

try
    disp('LW410Interface.waveData(): Writing user defined sequence into LW410.');
    
lw410 = LW410Interface();

%sampletime = 2.5e-9;
sampletime = lw410.getSampletime();

T = 1/f
t = [0:sampletime:maxSamples*sampletime];
seq = 4.2*sin(2*pi*f.*t);
seq(end) = 0;               %end hard set to zero

    plot(seq);
    data = lw410.convertData(seq);
    
    outputEnabled = 1;
    if outputEnabled == 1
        outputEnabled = 'on';
    else
        outputEnabled = 'false';
    end
    disp('wave_data ----------------------');
    disp(obj.g);
    fopen(obj.g);

    fprintf(obj.g, '%s\n', horzcat('OUTPut1 ', outputEnabled));
    fwrite(obj.g, ['WAVE:DATA ', data], 'char')
    fclose(obj.g);
catch exception
    disp(horzcat('LW410Interface.wave_data(): ', getReport(exception)));
    delete(obj.g);      %delete visa obj from memory in case there is error (too many samples, ...)
end
