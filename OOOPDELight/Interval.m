classdef Interval  < grid1D
     %Interval 1D FEM "Mesh"
     % Interval class implements an Interval I = [a,b] in means of p-e-t
     % structure.
    
    
    
    methods
        function obj = Interval(varargin)
            %%
            %
            % * Interval IN double[,double] OUT:self
            %
            % The method discretizes an 1D domain (interval).
            %
            % The first argument is a vector containing the boundaries of 
            % the interval, or a vector containing fixed points.
            % The second argument is the optional maximal meshwidth (hmax).
            % Default value is 0.1.
            %
            % If the distance between neighboured points is lesser than the
            % maximal meshwith, then the mesh will be refined until h <
            % hmax.
            %
            % Call
            %
            %       grid1D.interval([0,1])
            %            
            %       grid1D.interval([0,1],0.1)
            %
            if nargin ~= 0
                obj.interval(varargin{:});
            end
        end

    end
    
end

