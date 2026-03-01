classdef ProgressBar < handle
% class for command-line progress-bar notification.
    properties
        strPercentageLength;
        strDotsMaximum;
    end
    methods
        %--- constructor
        function this = ProgressBar()
             %% Initialization
            % Vizualization parameters
            this.strPercentageLength = 10;   %   Length of percentage string (must be >5)
            this.strDotsMaximum      = 10;   %   The total number of dots in a progress bar
        end
        %--- print method
        function run(this, msg)
            % This function creates a text progress bar. It should be called with a 
            % STRING argument to initialize and terminate. Otherwise the number corresponding 
            % to progress in % should be supplied.
            % INPUTS:   C   Either: Text string to initialize or terminate 
            %                       Percentage number to show progress 
            % OUTPUTS:  N/A
            % Example:  Please refer to demo_textprogressbar.m
            % Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
            % Version: 1.0
            % Changes tracker:  29.06.2010  - First version
            % Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/
            %% Main 
            persistent strCR;                %   Carriage return pesistent variable
            if isempty(strCR) && ~ischar(msg)
                % Progress bar must be initialized with a string
                error('The text progress must be initialized with a string!');
            elseif isempty(strCR) && ischar(msg)
                % Progress bar - initialization
                fprintf('\n%s', msg);
                strCR = -1;
            elseif ~isempty(strCR) && ischar(msg)
                % Progress bar  - termination
                strCR = [];  
                fprintf([msg '\n']);
            elseif isnumeric(msg)
                % Progress bar - normal progress
                msg = floor(msg);
                percentageOut = [num2str(msg) '%%'];
                percentageOut = [percentageOut repmat(' ',1,this.strPercentageLength-length(percentageOut)-1)];
                nDots = floor(msg/100*this.strDotsMaximum);
                dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,this.strDotsMaximum-nDots) ']'];
                strOut = [percentageOut dotOut];
                
                % Print it on the screen
                if strCR == -1
                    % Don't do carriage return during first run
                    fprintf(strOut);
                else
                    % Do it during all the other runs
                    fprintf([strCR strOut]);
                end
                
                % Update carriage return
                strCR = repmat('\b',1,length(strOut)-1);
                
            else
                % Any other unexpected input
                error('Unsupported argument type');
            end
        end
    end
end