%SIMULATOR a simple program with gui for experiments with Turing patterns.
% ... write here a short manual, especially for the nogui mode.

function  simulator(flag)
 
    if nargin == 0
        flag = '';
        % if no flag given, use GUI mode. Start it here
    end
    
    switch lower(flag)
        case {'-gui',''}
            fprintf('\tGUI Mode\n')
            % Start here GUI program.
        case '-nogui'
            fprintf('\n\tsimulator in NOGUI Mode\n')
            fprintf('----------------------------------------------\n');
            while true
                
                % Start here the nogui script.
                %  A user dialog would be nice. Maybe get arguments 
                % for nogui mode from varargin cell-array.
    
                % Ask the user for a dialog, implement a user dialog on
                % console.
                fprintf('Choose your model, or setup manually.\n');
                fprintf('Avaliable predefined models are:\n')
                fprintf('''Schnakenberg'' and ''Brusselator''.\n')
                fprintf('For ''Schnakenberg'' type s ')
                fprintf('and for ''Brusselator'' type b.\n')
                fprintf('To setup manually, type m.\n');
                fprintf('To leave the program type e or exit.\n');
                answer = input('>>>',"s");
               

                

                switch lower(answer)
                    % If the user decides to run the program w/o giving data, 
                    % run this:
                    
                    case {'schnakenberg' 's'}
                        model = setModel('Schnakenberg'); 
                        fprintf('You running Schnakenberg model ')
                        fprintf('with predefined setting.\n')
                        fprintf('To set up your own model use GUI or \n');
                        fprintf('choose option m in nogui mode.\n')
                        
                        % start to prepare model: Set FEM and create geometry.
                        fem = Bilinear2D;
                        geometry = RectangleR(0,100,0,50,model.h);            
                   
                        x = geometry.x; 
                        y = geometry.y;
                        % define the initial pulse.
                        p = zeros(geometry.nPoints,1); 
                        p(( (x-50).^2+(y-25).^2) <= 4) = 1;

                        % Run the model
                        % Initial condition u at t=0  and v at t=0. 
                        u0 = model.u0*ones(geometry.nPoints,1) + p;
                        v0 = model.v0*ones(geometry.nPoints,1);
                        
                        % call implizit Euler 
                        FEM.solveImplicit(fem,geometry,model,...
                            u0,v0,model.tend,model.dt); 
       
                        fprintf('Simulation finished.\n');

                    case {'brusselator' 'b'}
                        model = setModel('Brusselator'); 
                        fprintf('You running Brussselator model ')
                        fprintf('with predefined setting.\n')
                        fprintf('To set up your own model use GUI or \n');
                        fprintf('choose option m in nogui mode.\n')

                        % start to prepare model: Set FEM and create geometry.
                        fem = Bilinear2D;
                        geometry = RectangleR(0,100,0,50,model.h);            
                   
                        x = geometry.x; 
                        y = geometry.y;
                        % define the initial pulse.
                        p = zeros(geometry.nPoints,1); 
                        p(( (x-50).^2+(y-25).^2) <= 4) = -3;
                        % Run the model
                        % Initial condition u at t=0  and v at t=0. 
                        u0 = model.u0*ones(geometry.nPoints,1) + p;
                        v0 = model.v0*ones(geometry.nPoints,1);
                        
                        % call implizit Euler 
                        FEM.solveImplicit(fem,geometry,model,...
                            u0,v0,model.tend,model.dt);
                        fprintf('Simulation finished.\n');

                    case {'m' 'manually'}
                        %start dialog to assemble the arguments.
                        % Check the arguments the user gives here...
                    case {'e','exit'}
                        return
                    otherwise
                        fprintf(['\tUnknown model or option: "' answer,...
                            '". \n\tAllowed arguments are ''s'', ''g''',...
                            ' ''m'' ''e''.\n\n']);
                        % 2 seconds to read this...
                        pause(2) 
                end

                
            end
        case {'-h','-help'}
            fprintf('Print here help for this program.\n');
            % If you write a help in the header of this file, you can call
            % it here.
            help simulator
        otherwise
            fprintf(['\tUnknown argument flag: "' flag,...
                '" Allowed flags are -gui and -nogui.\n']);
            
    end     
end