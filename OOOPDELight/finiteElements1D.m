%% finiteElements1D
% Abstract class to build 1D Finite Elements
% 
% Provides element matrix prototypes

%% Inheritance
% finiteElements1D < finiteElements

classdef (Abstract) finiteElements1D < finiteElements    
    
    % Changes 06/29/2015
    % Bugfix: Q,G now correct initialized
    % H and R now stored in economy style as 2xn and 2x1 matrices.
    properties(Abstract,Constant)
        idx
    end
    methods(Static, Access = public)
        %% Static methods with Access = public
        %
        function [Q,G,H,R] = assemb(gridObj)
            % 
            % * assemb IN:grid1D OUT:double,double,double:double
            % 
            % Method that assembles the boundary matrices and vectors.
            %
            % Call:
            % 
            %       [Q,G,H,R] = finiteElements1D.assemb(grid1D)
            %
            
            % (c) 2013 Uwe Pruefert
            b = gridObj.b;
            rE =  gridObj.e(1,2);
            n = gridObj.nPoints;
            
            Q = sparse(n,n);
            H =  sparse(2,n);  
            R = sparse(2,1);
            G = sparse(n,1);
            
            % left BD
            m = b(2,1);
            lengthq = b(3,1);
            lengthg = b(4,1);
            
            if m == 0 % case Robin
                Q(1,1) = eval(char(b(5:5+lengthq-1,1)'));
                G(1) = eval(char(b(5+lengthq:5+lengthq+lengthg-1,1)'));
            else % case Dirichlet BC
                lengthh = b(5,1);
                lengthr = b(6,1);                 
                H(1,1) = eval(char(b(9:9+lengthh-1,1)'));
                R(1) = eval(char(b(9+lengthh:9+lengthh+lengthr-1,1)'));
            end
             
            % right BD
            
            m = b(2,2);
            lengthq = b(3,2);
            lengthg = b(4,2);           
            
            
            if m == 0 % case Robin
                Q(rE,rE) = eval(char(b(5:5+lengthq-1,2)));
                G(rE) = eval(char(b(5+lengthq:5+lengthq+lengthg-1,2)));               
            else % case   Dirichlet BC
                lengthh = b(5,2);
                lengthr = b(6,2); 
                H(2,rE) = eval(char(b(9:9+lengthh-1,2)'));                  
                R(2) = eval(char(b(9+lengthh:9+lengthh+lengthr-1,2)'));
            end              
        end  
        
        function [H] = periodic(gridObj)
            %%
            % * periodic IN:grid1D OUT:double
            %
            % Implements periodic boudary conditions in 1D.
            %
            % Call:
            %
            %       P = finiteElements1D.periodic(grid1D)
            %
            % Example:
            %
            % Add matrix multiplied by multipier >>1 
            % to linear system.
            % Let g and fem grid1D and lagrange11D objects.
            %
            %       [K,M,F] = fem.assema(g,1,1,f(g.p));
            %       P = fem.periodic(g);
            %       y = (K+M+1e3*P)\F
            %
            H = sparse(gridObj.nPoints,gridObj.nPoints);             
            
            H(1,gridObj.t(1,1)) = -1;
            H(1,gridObj.t(end-1,end)) = 1; 
            
            h = gridObj.h;
            
            H(gridObj.t(end-1,end),gridObj.t(1,1)) = -1/h(1);
            H(gridObj.t(end-1,end),gridObj.t(2,1)) = 1/h(1);
            
            H(gridObj.t(end-1,end),gridObj.t(1,end)) = 1/h(end);
            H(gridObj.t(end-1,end),gridObj.t(2,end)) = -1/h(end);    
        end
    end  
    
    methods(Access = public)
        %% Methods with Access = public
        %
        function [K,M,F] = createMatrixEntries(obj,gridObj,cf,af,ff) 
            %%
            % * createMatrixEntries IN:self,grid1D,double|char|function_handle,
            % double|char|function_handle,double|char|function_handle
            %
            % Method that computes the entries uses by assema method.
            % Arguments are the grid, and the coefficients c, a, f. 
            % Implements abtract method from superclass.
            sizeVector = sqrt(size(obj.S,1));
            sizeMatrix = size(obj.S,1);
            nt = gridObj.nElements;
            [cval,aval,fval] = obj.aCoefficients(gridObj,cf,af,ff);
            J = obj.makeJ(gridObj);  
            
            K1 = getConstantPartOfStiffnessMatrix(obj,gridObj);
           
            K = reshape(K1*sparse(1:gridObj.nElements,...
                1:gridObj.nElements,cval(1,:))...
                     ,1,gridObj.nElements*sizeMatrix);
        
            M = reshape(obj.M*(aval.*J),1,nt*sizeMatrix);            
            F = reshape(obj.F*(J.*fval),1,nt*sizeVector);            
        end
        
        function val = createConvectionEntries(obj,gridObj,b) 
            %%
            % * createConvectionEntries IN:self,grid1D,double|char|function_handle OUT:double
            %
            % Method that computes the entries uses by convection method.
            % Arguments are the grid, and the coefficient b. 
            % Implements abtract method from superclass.
            sizeMatrix = size(obj.C,1);
            J = obj.convCoefficients(gridObj,b); 
            val = reshape(obj.C*J,1,gridObj.nElements*sizeMatrix);        
        end       
    end
    
    methods(Access=public)
        function K1 = getConstantPartOfStiffnessMatrix(obj,gridObj)
            %%            
            % * getConstantPartOfStiffnessMatrix IN:self,grid1D OUT:double
            %
            % Method that computes the part of the stiffness matrix that 
            % depends not from the data.
            % Arguments is the grid. 
            %
            J = obj.makeJ(gridObj);             
            K1 = obj.S*(1./J);
        end
    end
    
    methods(Static,Access=protected)
        %% Static methods with Access = protected
        
        
        function [J] = makeJ(gridObj)  
            %
            % * makeJ IN:grid1D OUT double
            %
            % Computes the value of the Jacobian determinant.
            % Implements a abstract method from superclass.
            if ~isa(gridObj,'grid1D')
                gridObj.wrongClass.throwAsCaller;
            end
            J = gridObj.p(gridObj.t(2,:))-gridObj.p(gridObj.t(1,:));
        end
    end 
end

