  %% finiteElements
    %
    % Abstract class to build finite elements.
    %
    % Copy behavior = handle.
    % 
    %% Public, constant Properties:
    % *    wrongNumberInputs (MException)
    % *    wrongNumberOutputs (MException)
    % *    wrongFormat (MException)
    % *    wrongClass (MException)
    % *    wrongSize (MException)
    %
    %% Public, constant abstract property:    
    % * idx (double)
         
    
    
    % Change log: 
    % March 2017: Error in localErrorL2 corrected.
    % (c) 2015 by U.Prüfert
classdef (Abstract) finiteElements < handle    
    properties(Constant)
        % It defines only universal error Objects.
        wrongNumberInputs = MException('FINITEELEMENTS:WRONGNUMBERINPUTS',...
                        'The number of input arguments is wrong.');
        wrongNumberOutputs = MException('FINITEELEMENTS:WRONGNUMBERINPUTS',...
                        'The number of output arguments is wrong.');         
        wrongFormat = MException('FINITEELEMENTS:WRONGFORMAT',...
                                    'reading geometry data failed, maybe incorrect format.');
        wrongClass = MException('FINITEELEMENTS:WRONGCLASS',...
                                    'Wrong argument class.');
        wrongSize = MException('FINITEELEMENTS:WRONGSIZE',...
                                    'Wrong sized arguments.');                        
    end  
  
    properties(Abstract = true, Constant = true)
        idx 
    end
    
    methods(Access = public)
        %% Public methods
        % 
        % The methods assema and convection are actually not abstract, but
        % behave like abstract functions.
        %
        % * errorInd        
        
        %
        
        % Assembling of all Matrices valid for all ...        
        function [varargout] = assema(obj,gridObj,cf,af,ff)
            %%
            % * assema IN:self,gridd,double|char|function_handle, 
            %             double|char|function_handle, 
            %                 double|char|function_handle
            %
            % The assemble method for finite elements.
            % It assembles the Stiffness matrix K, the Mass Matrix M and the RHS 
            % vector F in sparse format. If you call it with two output arguments, then
            % S = K+M. 
            %
            % Call:
            %
            %       [K,M,F] = finiteElements.assema(grid,c,a,f)
            %       [S,F] = finiteElements.assema(grid,c,a,f)
            %
            %  
            
            
            % (c) 2013 by Uwe Prüfert
            % specialized methods dep. on element type
            scl = superclasses(gridObj);
            if isempty(scl) 
                obj.wrongClass.throw;
            end
            for k = 1:length(scl)                
                isgrid = strcmp(scl{k},'gridd');
                if isgrid
                    break;
                else
                    isgrid = false;
                end
            end
            if ~isgrid
                obj.wrongClass.throw;
            end
            [K,M,F] = obj.createMatrixEntries(gridObj,cf,af,ff);
            [idx0,idx1,idx2] =  obj.makeIndex(gridObj.t(obj.idx,:),gridObj.nElements); 
            % do the sparse magic...  
            switch nargout
                case 2  
                    varargout{1} = sparse(idx1,idx2,K+M,gridObj.nPoints,gridObj.nPoints); 
                    varargout{2} = sparse(idx0,1,F,gridObj.nPoints,1);
                case 3  
                    varargout{1} = sparse(idx1,idx2,K,gridObj.nPoints,gridObj.nPoints);  
                    varargout{2} = sparse(idx1,idx2,M,gridObj.nPoints,gridObj.nPoints);
                    varargout{3} = sparse(idx0,1,F,gridObj.nPoints,1);
                otherwise
                    throw(obj.wrongNumberOutputs)
            end
        end
        
        function B = convection(obj,gridObj,b)
            %%
            %
            % * convection IN:self,gridd,double|char|function_handle
            %
            % Method that computes the convection matrix.
            %
            % Call:
            %
            %       B = finiteElements.convection(gridd,bvec) 

            % (C) 2013 by Uwe Prüfert
        
            
            scl = superclasses(gridObj);
            for k = 1:length(scl)
                isgrid = strcmp(scl{k},'gridd');
                if isgrid
                    break;
                end
            end
            if ~isgrid
                obj.wrongClass.throwAsCaller;
            end          
            [~,idx1,idx2] = obj.makeIndex(gridObj.t(obj.idx,:),gridObj.nElements);
            val = obj.createConvectionEntries(gridObj,b);
            B = sparse(idx1,idx2,val,gridObj.nPoints,gridObj.nPoints);    
        end
        
        
        function B = sparsityPattern(obj,gridObj,varargin)
            %%
            % * sparsityPattern IN:self,grid[,double[,double,double,double,double]] OUT:double
            % 
            % Method that computes the sparsity pattern of the system
            % matrices.
            %
            % Call:
            %
            %        B = finiteElements.sparsityPattern(gridObject)
            %        B = finiteElements.sparsityPattern(gridObject,A)
            %        B = finiteElements.sparsityPattern(gridObject,K,M,C,Q,H)
            %       
            % A should be obj.A, intesionally for defining sparsity patterns
            % for systems of PDEs 
            % K,M,C,Q,H  should be the matrices from a call of obj.assema,
            % obj.assemb and obj.convection
            % The empty call computes all relevant matrices internally. 
            if ~isa(gridObj,'gridd')
                 obj.wrongClass.throw
            end
            
            switch length(varargin)
                case 0  % no matrices given
                    [K,M,~] = obj.assema(gridObj,'1','1','1');
                    [Q,~,H,~] = obj.assemb(gridObj);
                    C = obj.convection(gridObj,ones(size(gridObj.p(:,1))));
                    [indx1,indx2,s] = find(K+M+C+(H'*H)+Q);
                    B = sparse(indx1,indx2,s~=0);                    
                case 1  % whole matrix
                    [indx1,indx2,s] = find(varargin{1});
                    B = sparse(indx1,indx2,s~=0);                   
                case 5  % all relevant matrices in order K,M,C,Q,H
                    [indx1,indx2,s] = find(varargin{1}+...
                                           varargin{2}+...
                                           varargin{3}+...
                                           varargin{4}+...
                                           varargin{5}'*varargin{5});
                    B = sparse(indx1,indx2,s~=0); 
                otherwise
                    obj.wrongNumberInputs.throwAsCaller;
            end
        end
        
        
        
        
        function M = sourceTermMatrix(obj,grid) 
            %%
            % * sourceTermMatrix IN:self,gridd OUT: double
            %
            % Computes the matrix M such that F = M*f,
            % where f ist a vector of length #Elements which contrains the
            % values of the source f in the centers of each element. This
            % is usefull if you want to calculate the source from a given
            % vector, e.g.  when solving coupled systems where the solution
            % of one equation is the source for another equation.
            %
            % Call:
            %
            %        M = finiteElements.sourceTermMatrix(gridd)
            
            % (c) 2016 Uwe Pruefert.
           
            idx1 = reshape(grid.t(1:grid.nPointsInElements,:),...
                1,grid.nElements*grid.nPointsInElements);
            idx2 = reshape(repmat(1:grid.nElements,grid.nPointsInElements,1),...
                1,grid.nElements*grid.nPointsInElements); 
        
            M = sparse(idx2,idx1,reshape(obj.F*obj.makeJ(grid),1,grid.nPointsInElements*grid.nElements),grid.nElements,grid.nPoints)';
        end
        
        function errorPerElement = errorInd(obj,gridObj,y,c,a,f,alpha,beta,m)
            %%
            % * errorInd
            % IN:self,gridd,double,double,double[,double,double,double]
            % OUT:double
            %
            % Call:
            %
            %       errorPerElement = finiteElements.errorInd(grid,u,c,a,f[,alfa,beta,m])
            %
            % Method that computes the error for every element. The measure is
            % taken from pdetool. Note that errorInd works only for scalar
            % PDEs. For Systems, call it N times with the source and solution
            % component. Note that this may not really helpfull when applying
            % on coupled systems on non linear PDEs.
        
            switch nargin
                case 6
                    alpha = 1;   beta = 1;   m = 2;
                case 8
                    m = 2; 
                case 9
                    % okay
                otherwise
                    obj.wrongNumberInputs.throwAsCaller; 
            end
            % -------------------------------------------------------------
            % This may dimension be independent but may depend from
            % element-type. 
            % Compute areas and side lengths, method from grid-class
            [sideLength,~]= gridObj.sideLengthAndArea();
            
            % L2 norm of || f - a * u || over triangles f and a are element
            % data and u node. localError is abstract! The following calls
            % may depend from element-type and hence also from mesh -
            % extended or not.
            normFminusau = obj.localErrorL2(gridObj,y,a,f);
            
            % multiply by triangle's longes side 
            normFminusau = normFminusau.*max(sideLength).^m;
            
            %    flux jumps computed by assembly of ddncu into obj.nPoints
            %    x obj.nPoints sparse matrix jmps(i,j) becomes abs(jump
            %    across edge between nodes i and j). note that sparse(...)
            %    accepts duplicates of indices and performs summation !
            %    (i.e., the flux differences )
            
            
            % fluxes through edges, static and abstract method
            fluxThroughElementEdges = obj.fluxThroughEdges(gridObj,y,c);           
           
            
            % computeFluxJumps static and abstract method
            jumps = obj.fluxJumps(gridObj,fluxThroughElementEdges,m);
            % okay, all methods are abstract and static and should work for 1D--3D
            % -------------------------------------------------------------           
          
            errorPerElement = alpha*normFminusau + beta*sqrt(0.5*jumps);
        end
    end
    
    methods(Static,Access = public)
        %% Static methods with Access = public        
        
        function L = stiffSpring(M)
            %%
            % * stiffSpring IN:double OUT:double
            %
            % Computes a guess for the stiff-spring coefficient
            % L = obj.stiffSpring(M).
            % M should be the Systems Matrix, e.g. M = K+M+C
            % (C) 2013 by Uwe Prüfert
            L = 1e3*norm(M,1);
        end
    end
     
    methods(Access = public)    
        function disp(obj)
            disp(['    ',class(obj)])
        end
    end
     
    methods(Static,Access = protected) 
        %% Static methods with Access = protected
        
        function localL2 = localErrorL2(obj,u,a,f)
            %%
            % * localErrorL2 IN:self,double,double,double OUT:double
            %
            %
            % localErrorL2 is made for local error measurment. It may be
            % overwritten in higher order FE classes. However, for linear FEs
            % it works for every dimension
            % localErrorL2 computes the 
            % $L^2$
            % -norm of f-a*u
            % Evaluates f-a*u in the center of Element and multiplies with
            % area. Result is per element. This should work for all linear
            % P1 elements in 1D--3D. May be overwritten in P2 etc. classes.
            %
            % Call:
            %
            %       error = obj.localErrorL2(u,a,f).
            %
            % Note that u, a and f  must be vectors or scalars, not function
            % handles. 
            
            
            % (c) 2015 Uwe Prüfert
            % Changelog May 4th 2105: Make it dimension independent

            % This must be (abstract) in the gridxD classes. 
            % 
            % In 1D and 3D may area be
            % misleading: It is the volume of the lement and length is the
            % length or area of the edge or boundary 
            [~,area] = obj.sideLengthAndArea;

            if isscalar(f)
                f = ones(1,obj.nElements)*f;
            elseif min(size(f))==1&&max(size(f))==obj.nPoints
                f = obj.point2Center(f);
                f = f(:)';             
            elseif min(size(f))==1&&max(size(f))==obj.nElements
                f = f(:)';                
            else
                MException(obj.wrongInputFormatID,...
                    [obj.wrongInputFormatStr,...
                    ' must be vector of lenght np or ne,  or must be a  scalar.']).throwAsCaller;
            end
            if isscalar(a)
                a =ones(1,obj.nElements)*a;
            elseif min(size(a))==1&&max(size(a))==obj.nPoints
                a = obj.point2Center(a);
                a = a(:)';             
            elseif min(size(a))==1&&max(size(a))==obj.nElements
                a = a(:)';                
            else
                MException(obj.wrongInputFormatID,...
                    [obj.wrongInputFormatStr,...
                    ' must be vector of lenght np or ne,  or must be a  scalar.']).throwAsCaller;
            end 
           
            localL2 = sqrt((f-a.*obj.point2Center(u)).^2.*area);         
        end
        
        function [cval,aval,fval] =  aCoefficients(gridObj,cc,aa,ff)
            %%
            % * aCoefficients IN:grid1D,double|char|function_handle,
            % double|char|function_handle,double|char|function_handle
            % OUT:double,double,double
            %
            % Methods that compute the numerical values of the coefficients
            % c, a, and f from various formats.
            %
            [cval,aval,fval] = gridObj.aCoefficientsMpt(cc,aa,ff);
        end
        
        function bval = convCoefficients(gridObj,b)
            %
            % * convCoefficients IN:grid1D,double|char|function_handle            
            % OUT:double,double,double
            %
            % Methods that compute the numerical values of the coefficient
            % b from various formats. 
            %
            bval = gridObj.convCoefficientsMpt(b);
        end
    end
    
    
    
    % Declarations of abstract methods   
    methods(Abstract = true, Access = public) 
        %% Abtract methods with Access = public
        %
        [K,M,F] = createMatrixEntries(obj,gridObj,cf,af,ff);
        %%
        % * createMatrixEntries IN:self,grid1D,double|char|function_handle,
        % double|char|function_handle,double|char|function_handle
        %
        % Method that computes the entries uses by assema method.
        % Arguments are the grid, and the coefficients c, a, f. 
        %
        val = createConvectionEntries(obj,gridObj,b);
        %%
        % * createConvectionEntries IN:self,grid1D,double|char|function_handle OUT:double
        %
        % Method that computes the entries uses by convection method.
        % Arguments are the grid, and the coefficient b. 
        %  
    end
    
    methods(Abstract = true, Static = true, Access = public) 
        %% Static abtract methods with Access = public:    
        %        
    
        [idx0,idx1,idx2] = makeIndex(idx,nElements);  
        %%
        % * makeIndex IN:double,double OUT:double,double,double
        %
        % Method that computes a index structure to be uses to create
        % sparse matrices within assema, etc.
        %
        [Q,G,H,R] = assemb(gridObj); 
        %% 
        % * assemb IN:grid1D OUT:double,double,double:double
        % 
        % Method that assembles the boundary matrices and vectors.
        %
        % Call:
        % 
        %       [Q,G,H,R] = finiteElements1D.assemb(grid1D)
            %
    end
    
    methods(Abstract,Static,Access = protected)
        %% Static abtract methods with Access = protected:    
        % * fluxThroughEdges IN:gridd,double,double OUT:double
        % * fluxJumps IN:gridd,double,double OUT:double       
        
        % Functions for discretization error handling: adaptive error
        % indicator etc. NOTE: Tricky! These methods are static. The argument
        % gridObbj may differ from element-type to element-type.
        fluxThrougElementEdges= fluxThroughEdges(gridObj,u,c) 
        jumps = fluxJumps(gridObj,fluxThroughElementEdges,order)
        J = makeJ(gridObj);
        %%
        % * makeJ IN:grid1D OUT double
        %
        % Computes the value of the Jacobian determinant.
        %
    end    
end

