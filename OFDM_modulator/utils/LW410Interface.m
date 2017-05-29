%inerface for LW410
%   by LuboJ.
%
%Class is 'value class', so its values are passed by value, no reference,
%ie. if there is change in init() I have to write obj = obj.init() inside
%methods!
%
%If we would have it like in other programming languages, we have to extend
%handle class, in that case class is passed by reference to it.
classdef LW410Interface
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
        %   250ns   (4MHz)      YES!
        %   25ns    (40MHz)     YES!
        %   2.5ns   (400MHz)    YES!
        %   25ps                NO!
        sampletime = 250e-9;
        fs = 4e6;
        inSetMethod = 0;
        
        path = 'LW410conv'
    end
            
    methods
        %constructor
        function obj = LW410Interface( varargin)
            if nargin == 0
            elseif nargin == 1
                obj.gpib_addr = varargin{1};
            elseif nargin == 2
                obj.gpib_addr = varargin{1};
                obj.sampletime = varargin{2};
                obj.fs = 1/varargin{2};
            elseif nargin == 3
                obj.gpib_addr = varargin{1};
                obj.sampletime = varargin{2};
                obj.fs = 1/sampletime;
                obj.path = varargin{3};
            else
            end
        end
                    
        function obj = set.gpib_addr(obj, varargin)
            obj.gpib_addr = varargin{1};
        end
        
        % if condition to avoid infinite recursion
        function obj = set.sampletime(obj, newValue)
            obj.sampletime = newValue;
            if (obj.inSetMethod == 0)
                obj.inSetMethod = 1;
                obj.fs = 1/newValue;
                obj.inSetMethod = 0;
            end
        end
        
        % if condition to avoid infinite recursion
        function obj = set.fs(obj, newValue)
            obj.fs = newValue;
            if (obj.inSetMethod == 0)
                obj.inSetMethod = 1;
                obj.sampletime = 1/newValue;
                obj.inSetMethod = 0;
            end
        end
        
        function obj = set.g(obj, newValue)
            obj.g = newValue;
        end
        function g = get.g(obj)
            g = obj.g;
        end
       
        %parameters
        %   gpib_addr (default 0)
        %returns
        %   void
        function obj = init(obj)
            clear obj.g
            disp('init(): init obj.g');
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
        function [id] = idn(obj)
            %this is value class, means I need to save copy to change
            %values
            obj = obj.init()
            
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
            obj = obj.init()
            
            try
                data = obj.convertData(varargin{1});
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
                fwrite(obj.g,'WAVE:MARK:TYPE EDGE', 'char')
                fwrite(obj.g,'WAVE:MARK:EDGE:DEF', 'char')
                fwrite(obj.g,'WAVE:MARK:EDGE:TIME 0', 'char')
                %fwrite(obj.g,'WAVE:MARK:CLOCK:FREQ 1000', 'char')
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
        function [seq] = convertData(obj, varargin)
            disp('LW410Interface.waveData(): Converting sequence into DIFF format.');
            seq = varargin{1};
            command = [obj.path, 'LW410conv\LWConv.exe ', obj.path, 'LW410conv\temp.asc ', obj.path, '\temp.dif 0 USERSEQ USERSEQ ', num2str(obj.sampletime)];

            %if obj.path is not corrent user can choose path
            %if choosed path is not correct original one is
            %set back
            
            tmpPath = '';
            if ~exist([obj.path '\LWConv.exe'],'file')
                tmpPath = obj.path;
                obj.path = uigetdir('','Looking for LWConv.exe');
            end
            if exist([obj.path '\LWConv.exe'],'file')
                fid = fopen([obj.path '\temp.asc'],'w');
                fprintf(fid, '%f\n', seq);
                fclose(fid);

                % temp.asc -> temp.dif
                command = [obj.path, '\LWConv.exe ', obj.path, '\temp.asc ', obj.path, '\temp.dif 0 USERSEQ USERSEQ ', num2str(obj.sampletime)]
                dos(command);

                fid = fopen([obj.path, '\temp.dif'],'r');
                seq = fread(fid, [1 inf], 'uchar');   
                fclose(fid);            
            else
                %setting path back to default
                if (~empty(tmpPath))
                    obj.path = tmpPath;
                end
            end
        end
    end
    
end

