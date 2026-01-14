%% grid1D
% Class definition for 1D "grids".   
%
% Copy style = handle

%% Inheritance
%
% grid1D < gridd < handle 
%
%%
classdef(Abstract = true) grid1D < gridd
    % grid1D class definition for 1D "grids"   
     
    % Changes.
    % Constructor must be called withpout arguments.
    % 
    % interval replaces initMesh method
    % Prefered call now:
    % >> g = grid1D
    % >> g.interval([a,b])
    % >> g.interval([a,b],hmax)
    % Copy works now correct.
    % Arguments in plot work now correct.
    % Sept 2015: Some Correction when handle arguments in interval
    %            Throw an exception if number of elements is smaller 
    %            than the largest element number in to-refine list. 
    %            Correct code of method grid1D.h. 
     
     % Define the abstract prop 
     
    %% Properties with SetAccess = protected 
    %
    % * nPointsInElements   = 2
    % 
    properties(SetAccess = protected,...
            GetAccess = public)
        nPointsInElements = 2;
    end
    
    properties(Constant = true)
        spaceDimension = 1;
    end
    
    
    methods(Access = protected)
            
        function interval(obj,geo,hmax)
            %%
            %
            % * interval IN:self,double[,double] OUT:self
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
            %       grid1D.interval([0,1]),0.1)
            %
            
            
            % (c) 2013 Uwe Prüfert
           
            switch nargin
                case 2 % obj and geometry
                    if ~isa(geo,'double')
                        throw(obj.wrongClass)
                    end
                    obj.p = sort(geo(:))'; % row vector
                  
                    obj.t = [1:size(geo,2)-1 
                             2:size(geo,2)
                             ones(1,size(geo,2)-1)];
                    obj.e = [1 obj.nPoints
                                0 0
                                0 0
                                0 0
                                1 2
                                0 0]; 
                case 3                                    
                    try
                        geo = geo(:)';
                        if isa(hmax,'char')% hmax could be a char
                            hmax = eval(hmax); 
                        end                          
                        obj.interval(geo);                         
                        while obj.h > hmax
                            % refine                             
                            obj.refineMesh();                                                                            
                        end
                        
                    catch ME
                        ME.throwAsCaller;
                    end                    
                otherwise
                   obj.wrongNumberInputs.throw;
            end 
            obj.isExtended = false;
        end
    end
    %% Methods with Access = public
    %
    methods(Access = public)        
         % Make the following methods general in gridd?
        function diam = h(obj)
            %
            % diam = obj.h() 
            % Computes the meshsize h aka (in 2D)
            % the(outer) diameter of all triangles
            
            diam = obj.p(obj.t(end-1,:))-obj.p(obj.t(1,:));
             
        end
        
        
        
        function b = isBoundary(obj)
            %  b = isBoundary(gt) 
            % Method that indicates the boundary points.
            % b is a logical vector of length(#triangles)    %
                        
            b = false(1,size(obj.p,2));        
            b(obj.e(1,1))= true;
            b(obj.e(1,2))= true; 
        end  
        
        
        
        
        function refineMesh(obj,varargin)
            %%
            % 
            % * refineMesh IN:self[,double] OUT:self
            % 
            % Methods that  refinenes   grid.
            %
            % Call:
            %
            %         grid1D.refineMesh();
            %
            %         grid1D.refineMesh(elements);
            %
            % elements is an index vector of to refined elements. Its
            % length must be smaller than #elements.
            
            if obj.isExtended
                % Brute force re-initialization and re-extension
                x = obj.p(1:obj.ngpts);               
                obj2 = grid1D();
                obj2.interval(x);                 
                obj2.refineMesh(varargin{:});                
                obj2.extendMesh; 
                
                obj.p = obj2.p;
                obj.e = obj2.e;
                obj.t = obj2.t;

                obj.b = obj2.b;
                obj.ngpts = obj2.ngpts;
                obj.isExtended = obj2.isExtended;
                obj.nPointsInElements = obj2.nPointsInElements;                
            else                
                switch nargin
                    case 1                       
                        ind = 1:obj.nElements;
                    case 2                        
                        ind = varargin{1};
                        if max(ind)>obj.nElements
                            MException('GRID1D:REFINEMESH',...
                                ['The number of elements  and the',...
                                ' index given\nby the argument do not fit.']).throwAsCaller;
                        end
                    otherwise
                        obj.wrongNumberInputs.throw;            
                end
                 
                np = obj.nPoints;
                ne = obj.nElements;
                l = 1;
                for k = ind
                    obj.p(np+l) = 0.5*(obj.p(obj.t(1,k))+obj.p(obj.t(2,k)));
                    obj.t(1,ne+l) = np+l;
                    obj.t(2,ne+l) = obj.t(2,k);
                    obj.t(3,ne+l) = obj.t(3,k);
                    obj.t(2,k) = np+l;
                    l = l + 1  ;              
                end
            end
        end                
      
        function plot(obj,varargin)
            %%
            % 
            % * plot IN:self[,double[,char]] OUT:none
            %
            % plot method for class grid1D
            %
            % Call:
            %
            %       obj.plot()
            %
            %       obj.plot(y)
            %
            %       obj.plot(y,args)
            %
            % y can be a vector of length #elements or a vector of length
            % #points. In the first case, the function is assumed to be
            % piecewise constant, in the second case, it is assumed that
            % the function is piecewise linear.
            %
            % args can be an arbitrary number of pairs of arguments for
            % Matlab's common line object.
            %
            
            % (c) 2015 Uwe Prüfert
          
            if ~isempty(varargin) && isa(varargin{1},'matlab.graphics.axis.Axes')
                ax = newplot(varargin{1});
                if length(varargin)>1
                    varargin = varargin(2:end);
                else
                    varargin = {};
                end
            else
                ax = newplot();
            end

            if (nargin == 1) || (mod(length(varargin),2) == 0)
                % plot only mesh, maybe with pairs of options
                if obj.nPoints < 25
                    plot(obj.p,zeros(size(obj.p,2)),varargin{:});
                else                   
                    plot([min(obj.p) max(obj.p)],[1 1],varargin{:});
                end
            else
                y = varargin{1};
                if obj.nPoints < 200
                    switch obj.nPoints-length(y)
                        
                        case 0 % (continuous) linear FE                             
                            [x,indx] = sort(obj.p);
                            plot(x,y(indx),varargin{2:end});
                        case 1 % discontinuous P0 
                            
                            for k = 1:obj.nElements
                                line([obj.p(obj.t(1,k)),obj.p(obj.t(2,k))],...
                                    [y(k),y(k)],'LineStyle','-',varargin{1:end});
                            end
                        otherwise
                            obj.wrongNumberPoints.throwAsCaller;
                    end
                else
                    switch max(size(y))
                        case obj.nPoints
                            [x,indx] = sort(obj.p);
                            plot(x,y(indx),varargin{2:end});
                        case obj.nElements
                            x = obj.midpts;
                            plot(x,y,'.',varargin{2:end});
                        otherwise
                    end
                end
            end

            % switch length(varargin)
            %     case 1                     
            %         if obj.nPoints < 25
            %             plot(obj.p,zeros(size(obj.p,2)),'+');
            %         else
            %             plot([min(obj.p) max(obj.p)],[1 1]);
            %         end
            %     otherwise % with additional arguments
            %         if obj.nPoints < 200
            %             switch obj.nPoints-length(y)
            %                 case 0 % (continuous) linear FE                             
            %                     for k = 1:obj.nElements
            %                         line([obj.p(obj.t(1,k)),obj.p(obj.t(end-1,k))]',...
            %                             [y(obj.t(1,k)),y(obj.t(end-1,k))]',...
            %                             'LineStyle','-',varargin{1:end});
            %                     end
            %                 case 1 % discontinuous P0 
            % 
            %                     for k = 1:obj.nElements
            %                         line([obj.p(obj.t(1,k)),obj.p(obj.t(2,k))],...
            %                             [y(k),y(k)],'LineStyle','-',varargin{1:end});
            %                     end
            %                 otherwise
            %                     obj.wrongNumberPoints.throwAsCaller;
            %             end
            %         else % brute force plotting 
            %             switch max(size(y))
            %                 case obj.nPoints
            %                     [x,indx] = sort(obj.p);
            %                     plot(x,y(indx),varargin{:});
            %                 case obj.nElements
            %                     x = obj.midpts;
            %                     plot(x,y,'.',varargin{:});
            %                 otherwise
            %             end
            %         end
            % end
        end
        
        function obj2 = copy(obj1)
            %%
            %
            % * copy IN:self OUT:grid
            %
            % Hard copy method
            %
            % Call
            %
            %        obj2 = obj1.copy()
            %
            classArg1 = class(obj1);
            eval(['obj2 = ',classArg1,'();']);         
            obj2.p = obj1.p;
            obj2.e = obj1.e;
            obj2.t = obj1.t;
            obj2.b = obj1.b;
            obj2.ngpts = obj1.ngpts;
            obj2.isExtended = obj1.isExtended;
            obj2.nPointsInElements = obj1.nPointsInElements;
        end
           

        function bval = convCoefficientsMpt(obj,b)
            %%
            % * convCoefficientsMpt IN:self,double|char|inline|function_handle)
            % OUT:double
            % 
            % The method is used in the fem class methods to compte matrix 
            % entries for assembling C matrix.
            %
            % The method is used in the fem class method convection. The
            % argument must be an evaluable term or a double of length one,
            % #elements, or #points. 
            % The method return a
            % vector of lenght #elements.
            %
            % Call:
            %
            %       grid1d.convCoefficientsMpt(1)    
            %
            %       grid1d.convCoefficientsMpt('1')   
            %
            %       grid1d.convCoefficientsMpt(fun)
            %
            
            
            if isa(b,'function_handle') || isa(b,'inline')
                bval = feval(b,obj.x);
            elseif isa(b,'char')
                x = obj.point2Center(obj.p);  % if b depends from "x"...
                bval = eval(b).*ones(1,obj.nElements);
            elseif isa(b,'double')                
                switch length(b)
                    case obj.nPoints 
                        bval = obj.point2Center(b);
                    case obj.nElements                    
                        bval = b;    
                    case 1
                    % skalar                    
                        bval = b*ones(1,obj.nElements);
                    otherwise
                        obj.wrongFormat.throwAsCaller;
                end   
                
            else
                obj.wrongFormat.throwAsCaller;
            end    
        end
        
        function [cval,aval,fval] = aCoefficientsMpt(obj,c,a,f) 
            %%
            %
            % * aCoefficientsMpt IN self,c,a,f OUT:double,double,double
            % 
            % The method is used in the fem class methods to compte matrix 
            % entries for assembling K, M matrices and F vector.
            % The method conputes the coefficients for diffusion, linear 
            % part and source (c,a,f). c a f can be given as 
            % double, char, inline, or function_handle. All           
            % argument must be a evaluable terms or a doubles of length one,
            % #elements, or #points. Mixed formulations are allowed.
            % The method return three vectors of lenght #elements.
            %
            % Call:
            %
            %       [cval,aval,fval] = aCoefficientsMpt('1','0',10)
            % 
            %       [cval,aval,fval] = aCoefficientsMpt(1,1,10)
            % 
            %       [cval,aval,fval] = aCoefficientsMpt('1','sin(x)',@f)
            % 
          
            
            if isa(c,'function_handle') || isa(c,'inline')
                cval = feval(c,obj.point2Center(obj.p));
            elseif isa(c,'char')
                x = obj.point2Center(obj.p); %#ok<*NASGU>
                cval = eval(c).*ones(1,obj.nElements);
            elseif isa(c,'double')
                if length(c)==obj.nElements
                    cval = c(:)';
                elseif length(c)==obj.nPoints
                    cval = obj.point2Center(c);                    
                elseif isscalar(c)
                    cval = c*ones(1,obj.nElements);
                else
                    obj.wrongFormat.throwAsCaller;
                end
            elseif isa(c,'inline')
                cval = c(p);
            else
                obj.wrongFormat.throwAsCaller;
            end
            
            if isa(a,'function_handle') || isa(a,'inline')
                aval = feval(a,obj.point2Center(obj.p));                
            elseif isa(a,'char')
                aval = eval(a).*ones(1,obj.nElements);
            elseif isa(a,'double')
                if length(a)==obj.nElements
                    aval = a; 
                    aval = aval(:)';
                elseif length(a)==obj.nPoints
                    aval = obj.point2Center(a);
                elseif isscalar(a)
                    aval = a*ones(1,obj.nElements);
                else
                    obj.wrongFormat.throwAsCaller;
                end
            elseif isa(a,'inline')
                aval = a(p);
            else
                obj.wrongClass.throwAsCaller;
            end
            
            if isa(f,'function_handle') || isa(f,'inline')
                fval = feval(f,obj.point2Center(obj.p));
            elseif isa(f,'char')
                x = obj.point2Center(obj.p);
                fval = eval(f).*ones(1,obj.nElements);
            elseif isa(f,'double')
                if length(f)==obj.nElements                   
                    fval = f;
                elseif length(f)==obj.nPoints
                    fval = obj.point2Center(f);     
                elseif isscalar(f)
                    fval = f*ones(1,obj.nElements);
                else
                    obj.wrongFormat.throwAsCaller;
                end
            elseif isa(f,'inline')
                fval = f(p);
            else
                obj.wrongClass.throwAsCaller;
            end
        end
                       
        function element = pointToElementIndex(obj,pt)
            [~,indx] = sort(abs(obj.p-pt));
            element = indx(1);
        end
        
        function  dy = gradient(obj,f,varargin)
            dy = (f(2:end)-f(1:end-1))'...
                ./(obj.p(1,2:end)-obj.p(1,1:end-1));
        end
                
        function extendMesh(obj)
            %%
            %
            % * extendMesh IN:self, OUT:self
            %
            % Method that extends a 1D mesh for use Lagrange-2 finite
            % elements.
            %
            % Call:
            %
            %       grid1D.extendMesh();
            
            % (c) 2013 by Uwe Prüfert
            
            if obj.isExtended
                return
            end
            
            
            % midpoints Note the order of operations is importand to work.
                       
            % brute force re-ordering of points
            x = sort(obj.p);
            obj2 = Interval(x);                 
            obj.p = obj2.p;
            obj.t = obj2.t;
            obj.e = obj2.e;
             
            
            xad = 0.5*(obj.p(2:end)+obj.p(1:end-1));               
            t= [obj.t(1,:);...
                obj.nPoints+(1:length(xad));...
                obj.t(2,:);...
                obj.t(3,:)];
            obj.t = t;
            obj.ngpts = obj.nPoints;
            obj.p = [obj.p,xad];  
            obj.isExtended = true;    
            obj.nPointsInElements = 3;       
        end 
        
        function [sidelength,area] = sideLengthAndArea(obj)
        %Method to compute  side lengths and areas of elements.
        % Results are vectors of length nElements
        % (c) 2015 by Uwe Prüfert
        
            % Let them one 
            sidelength = ones(1,obj.nElements);
            % semi-last line minus first line of t are the left and right
            % points of the sub-intervals
            area =  obj.p(:,obj.t(end-1,:))-obj.p(:,obj.t(1,:)); 
            % Make it  compatible to extendedMeshs
            area = area(1:obj.nElements);            
        end  
    end  % public block ends
    
    methods(Access = private)
        function N = neighbours(obj,varargin)            
            % N = neigbours(gt[,indx]) 
            % computes to every triangle the 
            % index of neighbored triangles              
            nt = length(obj.t(1,:));            
            if nargin==2
                indx = varargin{1};
            else
                indx = 1:nt;
            end            
            N = sparse_null(3,nt);
            
            for k = 1:length(indx)
                nk = obj.ent(indx(k));
                [i] = find(nk==k);
                nk(i)=[];
                N(1:length(nk),indx(k))=nk;
            end
        end  
    end    
end
