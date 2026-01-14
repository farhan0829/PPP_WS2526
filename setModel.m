classdef setModel < handle
    %settings Class that stores the setting for the Brussleator 
    % and Schnakenberg model
    %  @Farhan: add detailed explanation here.    %   

    properties(Access = public)
        Au (1,1) double = 10
        Av (1,1) double = 0
        Buu (1,1) double = 1 %  The linear u-v coupling.
        Buv (1,1) double = 0 %  Can be seen as a matrix M
        Bvu (1,1) double = 0 %  with M = [Buu Buv
        Bvv (1,1) double = 0 %            Bvu Bvv]
        Cu (1,1) double = 0    
        Cv (1,1) double = 0
        Du (1,1) double {mustBePositive} = 1
        Dv (1,1) double {mustBePositive} = 1
        
        u0 (:,1) double
        v0 (:,1) double
        
        dt (1,1) double {mustBePositive} = 0.1  % time step size
        h (1,1) double {mustBePositive} = 1   % spatial step size
        tend (1,1) double {mustBePositive} = 300;   % time to run

        stepsPerFrame (1,1) double {mustBeInteger} = 10;
    end

    methods(Access = public)
        function obj = setModel(arg)  
            %setCoefficients Constructor method for the setCoefficints class.
            % Set-up the default values for coefficients of Schnakenberg 
            % and Brusselator model.
            % If you use it with parameters, 'schnakenberg' and
            %  'brusselator' are allowed. 
            % Calling w/o arguments, all coefficients are zero.
            switch nargin
                case 0
                    % default constructor
                case 1
                    switch lower(arg)
                        case 'schnakenberg'
                            obj = obj.schnakenberg();
                        case 'brusselator'
                            obj = obj.brusselator();
                        otherwise
                            MException('setCoefficients:unknownArgument',...
                                ['The argument ',char(arg),' is unknown.',...
                                ' Valid are ''schnakenberg'' or ''brusselator''.']).throw;
            
                    end
                otherwise
                    % too many argumnts: MATLAB will trow an 
                    % exception for us.
            end                    
        end
    end

    methods(Access = public,...
            Static = true)
        function obj = brusselator()
            %BRUSSELATOR sets the values for default Brusselator model.
            %   @Farhan: add detailed explanation here
            obj = setModel();
            obj.Au = 2;         obj.Av = 0; % fix
            obj.Buu = -4;       obj.Buv = 0; % Buu = -(Bvv+1)
            obj.Bvu = 3;        obj.Bvv = 0;
            obj.Cu = 1;         obj.Cv = -1;% fix
            obj.Du = 1;         obj.Dv = 8; % Du = 1 fix

            obj.u0 = obj.Au;
            obj.v0 = obj.Bvu/obj.Au;

            obj.h = 1;
            obj.dt = 1;
            obj.tend = 2000;

            obj.stepsPerFrame = 5; 
        end

        function obj = schnakenberg()
            %SCHNAKENBERG sets the values for default Schnakenberg model.
            %   @Farhan: add detailed explanation here
            obj = setModel();
            obj.Au = 0.01;      obj.Av = 2;    
            obj.Buu = -1;       obj.Buv = 0;  
            obj.Bvu = 0;        obj.Bvv = 0;
            obj.Cu = 1;         obj.Cv = -1;   % fix
            obj.Du = 1;         obj.Dv = 100;  % Du = 1 fix

            obj.u0 = obj.Au+obj.Av;
            obj.v0 = obj.Av/(obj.Au+obj.Av)^2;

            obj.h = 1;
            obj.dt = 0.1;
            obj.tend = 800;
            obj.stepsPerFrame = 2; 
        end

        function obj = custom(Du,Dv,Cu,Cv,Bu,Bv,Eu,Ev,Au,Av)
            %CUSTOM sets the values of the coefficients in the generalized
            %equn.
            %   @Farhan: add detailed explanation here. Important: Describe
            %   here the order of coefficients. Will be shown when calling
            % help setCoefficients.custom
            obj.Au = Au;
            obj.Av = Av;
            obj.Bu = Bu;
            obj.Bv = Bv;
            obj.Cu = Cu;
            obj.Cv = Cv;
            obj.Du = Du;
            obj.Dv = Dv;
            obj.Eu = Eu;
            obj.Ev = Ev;
        end
    end
end