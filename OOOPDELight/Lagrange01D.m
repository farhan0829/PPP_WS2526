%% lagrange01D
%
% Class that implements Lagrange-0 elements. The only method is mass.

%% Inheritance
% Lagrange01D < finiteElements1D
classdef Lagrange01D < finiteElements1D
    %class stub that implements L0 elements     
    
    properties(Constant)        
        idx = [];
    end
     
    
    methods(Static,Access = public)
        %% Static methods with Access = public
        %
        function M = mass(grid)
            %%
            % * mass IN:grid1D OUT:double
            %
            % Computes the mass matrix for Lagrange-0 elements. 
            %
            % Call:
            %
            %       M =Llagrange01D.mass(grid1D)
            %
            j =  Lagrange01D.makeJ(grid);            
            M = sparse(1:grid.nElements,...
                1:grid.nElements,j);  
        end  
        function makeIndex()
        end 
    end
    
    methods(Static,Access = protected)
        % !!!! Dummies !!!! Abstract methods, need an "implementation"
         % !!!! Dummies !!!!
        
        function fluxThroughEdges()
        end
        function localErrorL2() 
        end
        function fluxJumps()
        end
    end
end