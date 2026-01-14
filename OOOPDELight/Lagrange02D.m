%% Lagrange02D
%
% Class that implements Lagrange-0 elements. The only method is mass.

%% Inheritance
% Lagrange02D < finiteElements2D
classdef Lagrange02D < finiteElements2D
    %class stub that implements L0 elements     
    
    properties(Constant)       
        idx = [];
    end
     
    
    methods(Static,Access = public)
        %% Static methods with Access = public
        %
        function M = mass(grid)
            %%
            % * mass IN:grid2D OUT:double
            %
            % Computes the mass matrix for Lagrange-0 elements. 
            %
            % Call:
            %
            %       M = Lagrange02D.mass(grid2D)
            %
            j =  Lagrange02D.makeJ(grid);            
            M = 0.5*sparse(1:grid.nElements,...
                1:grid.nElements,j);  
        end
    end
    methods(Static,Access = public) 
        % !!!! Dummies !!!!
        function makeIndex()
        end 
    end
    
    methods(Static,Access = protected) 
        % !!!! Dummies !!!!
          
        
        function fluxThroughEdges()
        end
        function localErrorL2() 
        end
        function fluxJumps()
        end
    end
end

