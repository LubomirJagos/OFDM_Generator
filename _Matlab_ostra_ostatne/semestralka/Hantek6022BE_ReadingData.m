function [ch1Data, ch2Data, triggerPointIndex] = Hantek6022BE_ReadingData( ...
    deviceIndex, ...
    channel,    ...
    nLength, ...
    voltCh1Div, ...
    voltCh2Div, ...
    sampleFreq, ...
    nTrigSweep, ...
    nTrigSrc, ...
    nTrigLevel, ...
    nSlope, ...
    nHTrigPos, ...
    nDisLen, ...
    nTrigPoint, ...                 %pointer to store index of triggered point
    nInsertMode)

    disp('Matlab Hatek 6022BE scope control')

%     deviceIndex = 0

%     channel = 0
        % 0 = Ch1
        % 1 = Ch2

%     nLength = 100
        % samples to read

%     voltCh1Div = 4
%     voltCh2Div = 4
        % 0 = 20mV
        % 1 = 50mV
        % 2 = 100mV
        % 3 = 200mV
        % 4 = 500mV
        % 5 = 1V
        % 6 = 2V
        % 7 = 5V
    if (voltCh1Div == 0) voltUnit1 = 0.00078125;
    elseif (voltCh1Div == 1) voltUnit1 = 0.0015625;
    elseif (voltCh1Div == 2) voltUnit1 = 0.003125;
    elseif (voltCh1Div == 3) voltUnit1 = 0.00625;
    elseif (voltCh1Div == 4) voltUnit1 = 0.01625;
    elseif (voltCh1Div == 5) voltUnit1 = 0.01625; %0.0325;  %newValue correct??
    elseif (voltCh1Div == 6) voltUnit1 = 0.065;
    elseif (voltCh1Div == 7) voltUnit1 = 0.16;
    end
    if (voltCh2Div == 0) voltUnit2 = 0.00078125;
    elseif (voltCh2Div == 1) voltUnit2 = 0.0015625;
    elseif (voltCh2Div == 2) voltUnit2 = 0.003125;
    elseif (voltCh2Div == 3) voltUnit2 = 0.00625;
    elseif (voltCh2Div == 4) voltUnit2 = 0.01625;
    elseif (voltCh2Div == 5) voltUnit2 = 0.01625; %0.0325;  %newValue correct??
    elseif (voltCh2Div == 6) voltUnit2 = 0.065;
    elseif (voltCh2Div == 7) voltUnit2 = 0.16;
    end
    
%     sampleFreq = 13
        % 0 - 10 = 48MHz
        % 11 = 16MHz
        % 12 = 8MHz
        % 13 = 4MHz
        % 14 - 24 = 1MHz
        % 25 = 500kHz
        % 26 = 200kHz
        % 27 = 100kHz

    % generating time axis for plot
    sampleFreqValue = Hantek6022BE_GetSampleFreq(sampleFreq)        
    sampleTimeValue = 1/sampleFreqValue
    t = [0:sampleTimeValue:(nLength-1)*sampleTimeValue];

    calLevel = libpointer('int16Ptr',zeros(1, nLength));

    pCh1Data = libpointer('int16Ptr', zeros(1, nLength));
    pCh2Data = libpointer('int16Ptr', zeros(1, nLength));

    %   PARAMS PASSED AS ARGUMENTS!!!
    %nTrigSweep = 0
    %nTrigSrc = 0
    %nTrigLevel = 0 
    %nSlope = 0
        % 0 = raise slope
        % 1 = fall slope
    %nHTrigPos = 0

    %nDisLen = 0            %data to display

    nTrigPoint = libpointer('uint32Ptr', 0);
    %nInsertMode = 0
        % 0 = step
        % 1 = line
        % 2 = sin(x)/x


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [w n] = loadlibrary('c:\_hantek\HTMarch.dll','c:\_hantek\HTMarch.h','alias','hantek')
    libfunctions hantek

    r = calllib('hantek', 'dsoOpenDevice', deviceIndex)
    r = calllib('hantek', 'dsoGetCalLevel', deviceIndex, calLevel, nLength)

    %stem(calLevel.Value)

    r = calllib('hantek', 'dsoSetVoltDIV', deviceIndex, 0, int32(voltCh1Div))
    r = calllib('hantek', 'dsoSetVoltDIV', deviceIndex, 1, int32(voltCh2Div))
    r = calllib('hantek', 'dsoSetTimeDIV', deviceIndex, int32(sampleFreq))

        nTrigSweep
        nTrigSrc
        nTrigLevel
        nSlope
        
    r = calllib('hantek', 'dsoReadHardData',    ...
        deviceIndex,    ...
        pCh1Data,       ...
        pCh2Data,       ...
        nLength,        ...
        calLevel,       ...
        int32(voltCh1Div),     ...
        int32(voltCh2Div),     ...
        int16(nTrigSweep),     ...
        int16(nTrigSrc),       ...
        int16(nTrigLevel),     ...
        int16(nSlope),         ...
        int32(sampleFreq),     ...
        nHTrigPos,      ...
        nDisLen,        ...
        nTrigPoint,     ...
        nInsertMode)

    %r = calllib('hantek', 'dsoReadHardData', deviceIndex, pCh1Data, pCh2Data, nLength, voltCh1Div, voltCh2Div, nTrigSweep, nTrigSrc, nTrigLevel, nSlope, sampleFreq, nHTrigPos, nDisLen, nTrigPoint, nInsertMode)

    % y1 = double(pCh1Data.Value) .* voltUnit1;       %conversion needed, otherwise it's rounding value => NOT GOOD!
    % y2 = double(pCh2Data.Value) .* voltUnit2;
%     y1 = pCh1Data.Value;                          %original sample values
%     y2 = pCh2Data.Value;

%     subplot(211); stem(t, y1); grid on; hold on;
%     subplot(212); stem(t, y2, 'redx'); grid on; hold off;

    ch1Data = pCh1Data.Value;
    ch2Data = pCh2Data.Value;
%     ch1Data = double(pCh1Data.Value) .* voltUnit1;       %conversion needed, otherwise it's rounding value => NOT GOOD!
%     ch2Data = double(pCh2Data.Value) .* voltUnit2;

    %Hantek store index of point which triggered sensing into output param
    triggerPointIndex = nTrigPoint.Value;
end













