%inerface for LW410
%   by LuboJ.
%
%   Experiment to define interface as handle class, but it's useless since,
%   expect there would be more LW410 connected to the same PC. In that case
%   instances share same properties so it's NOT POSSIBLE!

classdef LW410InterfaceHandle < handle
    properties
        g
        gpib_addr = 1
        
        %Maximum samples to store into LW410:
        %(1MB memory / 4B (per double sample) = 262 083 samples (header also take some space))
        %experimentally tested
        
        %for now figure out:
        %   250us   (4kHz)      NO!
        %   25us    (40kHz)     YES!
        %   2.5us   (400kHz)    NOT TESTED!
        %   250ns   (4MHz)      NOT TESTED!
        %   25ns    (40MHz)     YES!
        %   2.5ns   (400MHz)   YES!
        %   25ps                NO!
        sampletime = 2.5e-6; %2.5e-9
        freq = 400e3;
        lwConvPath = 'c:\LeCroy_LW410'
    end
    
    methods
        %constructor
        function obj = LW410Interface(varargin)
            if nargin == 0
            elseif nargin == 1
                gpib_addr = varargin{1};
            elseif nargin == 2
                gpib_addr = varargin{1};
                sampletime = varargin{2};
            elseif nargin == 3
                gpib_addr = varargin{1};
                sampletime = varargin{2};
                lwConvPath = varargin{3};
            else
            end
        end
            
        function [sampletime] = get.sampletime(obj)
            sampletime = obj.sampletime;
        end
        
        function obj = set.sampletime(obj, newValue)
            if (newValue == 25e-6 ||  ...
                newValue == 25e-9 ||  ...
                newValue == 2.5e-9)
                obj.sampletime = newValue;
            else
                disp('Sampletime has to be 25us, 25ns or 2.5ns');
            end
        end
        
        %parameters
        %   gpib_addr (default 1)
        %returns
        %   void
        function obj = init(obj)
            clear g;
            disp('init(): init g');
            obj.g = visa('AGILENT', ['GPIB0::', num2str(obj.gpib_addr), '::INSTR']);

            % Set the property values.
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
            set(obj.g, 'Name', ['VISA-GPIB0-', num2str(obj.gpib_addr)]);
            set(obj.g, 'OutputBufferSize', 1100e3);
            set(obj.g, 'OutputEmptyFcn', '');
            set(obj.g, 'PrimaryAddress', obj.gpib_addr);
            set(obj.g, 'RecordDetail', 'compact');
            set(obj.g, 'RecordMode', 'overwrite');
            set(obj.g, 'RecordName', 'record.txt');
            set(obj.g, 'SecondaryAddress', 0);
            set(obj.g, 'Tag', '');
            set(obj.g, 'Timeout', 50);
            set(obj.g, 'TimerFcn', '');
            set(obj.g, 'TimerPeriod', 2);
            set(obj.g, 'UserData', []);     
        end
        
        %parameters
        %returns
        %   device ID
        function id = idn(obj)
            %this is value class, means I need to save copy to change
            %values
            obj.init();
            
            try
                disp('LW410Interface.idn(): Getting device ID.');
                %communication with device is treated like writing into file
                fopen(obj.g);
                % test IDN, set OUT ON
                fprintf(obj.g, '%s\n', '*IDN?');
                id = fscanf(obj.g, '%s');
                fclose(obj.g);
                disp(['LW410Interface.idn(): Device ID: ' id]);
            catch exception
                id = '';
                disp(horzcat('LW410Interface.idn(): ', getReport(exception)));
            end
        end
        
        %parameters
        %   wave_data
        %   outputEnabled
        %       1 == 'ON'
        %       0 == 'OFF'
        %returns
        %   void
        function [] = wave_data(obj, varargin)
            %this is value class, means I need to save copy to change
            %values
            obj.init();
            
            try
                data = convertData(varargin{1});
                outputEnabled = varargin{2};
                if outputEnabled == 1
                    outputEnabled = 'on';
                else
                    outputEnabled = 'false';
                end
                disp('LW410Interface.waveData(): Writing user defined sequence into LW410.');
                disp('----------------------- wave_data ----------------------');
                disp(obj.g);
                fopen(obj.g);
                
                fprintf(obj.g, '%s\n', horzcat('OUTPut1 ', outputEnabled));
                fwrite(obj.g, ['WAVE:DATA ', data], 'char')
                fclose(obj.g);
            catch exception
                disp(horzcat('LW410Interface.wave_data(): ', getReport(exception)));
                delete(obj.g);      %delete visa obj from memory in case there is error (too many samples, ...)
            end
        end
        
        %parameters
        %   sampletime
        %returns
        %   [seq] converted data
        function seq = convertData(obj, varargin)
            disp('LW410Interface.waveData(): Converting sequence into DIFF format.');
            seq = varargin{1};
            command = [obj.lwConvPath, '\LWConv.exe ', obj.lwConvPath, '\temp.asc ', obj.lwConvPath, '\temp.dif 0 USERSEQ USERSEQ ', num2str(obj.sampletime)];

            %if path is not corrent user can choose path
            %if choosed path is not correct original one is
            %set back
            
            tmpPath = '';
            if ~exist([obj.lwConvPath '\LWConv.exe'],'file')
                tmpPath = obj.lwConvPath;
                obj.lwConvPath = uigetdir('','Looking for LWConv.exe');
            end
            if exist([obj.lwConvPath '\LWConv.exe'],'file')
                fid = fopen([obj.lwConvPath '\temp.asc'],'w');
                fprintf(fid, '%f\n', seq);
                fclose(fid);

                % temp.asc -> temp.dif
                command = [obj.lwConvPath, '\LWConv.exe ', obj.lwConvPath, '\temp.asc ', obj.lwConvPath, '\temp.dif 0 USERSEQ USERSEQ ', num2str(obj.sampletime)]
                dos(command);

                fid = fopen([obj.lwConvPath, '\temp.dif'],'r');
                seq = fread(fid, [1 inf], 'uchar');   
                fclose(fid);            
            else
                %setting path back to default
                if (~empty(tmpPath))
                    obj.lwConvPath = tmpPath;
                end
            end
        end
    end
    
end

