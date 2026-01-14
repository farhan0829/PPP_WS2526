classdef(Abstract = true) grid2DR < gridd
    %grid2DR Class for 2D rectangle meshes
    %   Minimalistic implementation.
    
    properties(SetAccess = protected,...
            GetAccess = public)
        nPointsInElements = 4;
    end
    
    properties(Access = public)
        xgrid Interval
        ygrid Interval
    end
    
    properties(Constant = true)
        spaceDimension  = 2;
    end

    properties(Dependent = true,...
            SetAccess = protected)
         
        hmax
        hmin
        hmean
        a
    end

    methods
        function val = get.hmax(obj)
            val = max(max(obj.xgrid.h),max(obj.ygrid.h));
        end
        function val = get.hmin(obj)
            val = minx(min(obj.xgrid.h),min(obj.ygrid.h));
        end
        function val = get.hmean(obj)
            val = mean(mean(obj.xgrid.h),mean(obj.xgrid.h));
        end
        function val = get.a(obj)
            val = obj.xgrid.h'*obj.ygrid.h;
        end
    end
    
    methods(Access = public)
        function obj = rectangle(obj,x,y,h)
            %rectangle 
            switch nargin
                case 1
                    if isempty(obj.xgrid)||isempty(obj.ygrid)
                        obj.xgrid = Interval();
                        obj.ygrid = Interval();
                        % else: We use xgrid and ygrid.
                    end
                case 3
                    obj.xgrid = Interval(x);
                    obj.ygrid = Interval(y);
                case 4
                     obj.xgrid = Interval(x,h);
                     obj.ygrid = Interval(y,h);
            end
            nx = obj.xgrid.nPoints;
            ny = obj.ygrid.nPoints;
            
            x = sort(obj.xgrid.x);            
            y = sort(obj.ygrid.x);    
            [x,y] = meshgrid(x,y);
            obj.p = [reshape(x,1,nx*ny);reshape(y,1,nx*ny)];            
            
            obj.e = [[1:ny:(nx-1)*ny
                (1:ny:(nx-1)*ny)+ny
                zeros(2,nx-1)
                ones(1,nx-1)],...
                [ny*(nx-1)+1:ny*nx-1
                ny*(nx-1)+2:ny*nx
                zeros(2, ny-1)
                2*ones(1,ny-1)],...
                [nx*ny:-ny:ny+1
                    (nx*ny:-ny:ny+1)-ny
                    zeros(2,nx-1)
                3*ones(1,nx-1)],...
                [ny:-1:2
                ny-1:-1:1
                zeros(2,ny-1)
                4*ones(1,ny-1)]];
        
            idx = find(obj.e(5,:)==1);
            s1 = [[0, cumsum(sqrt(sum((obj.p(:,obj.e(2,idx(1:end-1)))-obj.p(:,obj.e(1,idx(1:end-1)))).^2,1)))]
            [cumsum(sqrt(sum((obj.p(:,obj.e(2,idx))-obj.p(:,obj.e(1,idx))).^2,1)))]];
            s1 = 1/max(max(s1))*s1;
             idx = find(obj.e(5,:)==2);
                s2 = [[0, cumsum(sqrt(sum((obj.p(:,obj.e(2,idx(1:end-1)))-obj.p(:,obj.e(1,idx(1:end-1)))).^2,1)))]
                [cumsum(sqrt(sum((obj.p(:,obj.e(2,idx))-obj.p(:,obj.e(1,idx))).^2,1)))]];
            s2 = 1/max(max(s2))*s2;
             idx = find(obj.e(5,:)==3);
                s3  = [[0, cumsum(sqrt(sum((obj.p(:,obj.e(2,idx(1:end-1)))-obj.p(:,obj.e(1,idx(1:end-1)))).^2,1)))]
                [cumsum(sqrt(sum((obj.p(:,obj.e(2,idx))-obj.p(:,obj.e(1,idx))).^2,1)))]];
            s3 = 1/max(max(s3))*s3;
             idx = find(obj.e(5,:)==4);
                s4 = [[0, cumsum(sqrt(sum((obj.p(:,obj.e(2,idx(1:end-1)))-obj.p(:,obj.e(1,idx(1:end-1)))).^2,1)))]
                [cumsum(sqrt(sum((obj.p(:,obj.e(2,idx))-obj.p(:,obj.e(1,idx))).^2,1)))]];
            s4 = 1/max(max(s4))*s4;
            
            obj.e(3:4,:)=[s1 s2 s3 s4];
            
            obj.t  = zeros(5,(nx-1)*(ny-1));
            
            for k = 1:nx-1
                obj.t(:,(k-1)*(ny-1)+(1:ny-1)) = (k-1)*(ny)+[ 1:ny-1
                 ny+1:2*ny-1
                 ny+2:2*ny
                 2:ny
                 zeros(size(1:ny-1))]; 
            end
            obj.t(5,:) = ones(1,obj.nElements);
        end
        
        function refineMesh(obj)
            if isempty(obj.xgrid)||isempty(obj.ygrid)
                MException('grid2DR:EMPTYMESH',...
                    'It make no sense to refine empty meshes.').throwAsCaller;
            end
            obj.xgrid.refineMesh;
            obj.ygrid.refineMesh;
            obj.rectangle;
        end
        
        function plot(obj,varargin)
            %%
            %
            % * plot IN:self[,double[,char]] OUT:none
            % 
            % Plot method for two dimensional grids. 
            %
            % Without arguments, the method plots the two dimensional grid
            % in the x-y-plane. Optional arguments are a double vector of
            % length #elements or of length #points and/or pairs of 
            % charachter array containing options that Matlabs patch class 
            % objects understand.
            %
            % Call:
            %
            %       obj.plot()
            % 
            % Plots the grid.
            %
            %       obj.plot(y)
            %
            % Here y must be a double vector of length #elements or of 
            % length #points. In the first case, the function y is assumed
            % to be piecewise constant,in the second case, it is assumed
            % that y ist piecewise linear.
            %
            %       obj.plot(arglist)
            %
            % arglist must  be paise parameter-value of options that 
            % Matlabs patch class objects understand.
            %
            % Example:
            %
            %       grid2D.plot('EdgeColor','k')
            %
            % plots the grid in black.
            % 
            %        obj.plot(y,'EdgeColor','k')
            % 
            % plots the data y over the grid as described above.
            %
            %       g.plot(sin(g.x.*g.y),'EdgeColor','k','FaceColor','none')
            %
            % plots the function 
            % $sin(x\cdot y)$
            % over the grid. The edge color is set to black, the faces are
            % deleted. 

           
             
            if ~isempty(varargin) && isa(varargin{1},'matlab.graphics.axis.Axes')
                ax = varargin{1} ;
                if length(varargin)>1
                    varargin = varargin(2:end);
                else
                    varargin = {};
                end
            else
                ax = newplot();
            end
             
             
            % Plot the  mesh triangles, no solution vector given.
            if (isempty(varargin)) || (mod(length(varargin),2) == 0)
                x = obj.p(1,:);
                y = obj.p(2,:);
                z = zeros(size(x));
                patch(ax,...
                    'faces',obj.t(1:4,:)',...
                    'vertices',[x(:),y(:),z(:)],...
                    'facecolor',[1 1 1],...
                    'edgecolor',get(ax,'DefaultSurfaceEdgeColor'),...
                    'parent',ax,...
                    varargin{:});
                view(ax,2);
            else
                x = obj.p(1,:);                
                y = obj.p(2,:);
                z = full(varargin{1});
                 
                
                % PO Plot 
                if length(z) == obj.nElements
                    x = reshape(obj.p(1,obj.t(1:4,:)),4,obj.nElements);
                    y = reshape(obj.p(2,obj.t(1:4,:)),4,obj.nElements);
                    z = [1;1;1;1]*z(:)';
                    cz = z;
                    patch(ax,...
                        x,y,z,cz,...
                        'lineStyle','none',...
                        varargin{2:end});
                     
                    grid(ax,'on');
                     
                % P1 or p2 plot
                elseif length(z) == obj.nPoints
                     
                     
                %  P1 plot
                   patch(ax,...
                       'faces',obj.t(1:4,:)','vertices',[x(:),y(:),z(:)],...
                        'facevertexcdata',z(:),...
                        'facecolor','interp',...
                        'edgecolor',get(ax,'DefaultSurfaceEdgeColor'),...
                        'parent',ax,...
                    varargin{2:end});
                     
                    grid(ax,'on');
                     

                else
                    MException('GRID2D:PLOT:WRONGELEMENT',...
                        'Elements can only be P0,P1 or P2.').throwAsCaller;
                end
            end
             
        end
        
        function identifyBoundarySegment(obj,nSegment)
            %% identifyBoundarySegment
            % 
            % IN:self[,double] OUT:none 
            %
            % Method that plot the grid and marks the boundary segments by colors.
            % The optional argument is the number of boundary segment to be marked.
            % If no argument given, all boundary segments were market by different colors.
            % If the segment not exists, nothing will be highlighted.
            %
            % Example:
            %
            %       grid2D.identifyBoundarySegment
            %
            % plots the grid and marks the four boundary segments in black,
            % red, green and blue. 
            %
            %       grid2D.identifyBoundarySegment(2)
            %
            % plots the grid and marks boundary segment no. two by red,
            % while all other boundaries are colored blue.
            
            
%             obj.plot;
            if nargin == 1 % on arguments    
                clf
                obj.plot;
                colors = [0 0 0
                          1 0 0
                          0 .6 0
                          0 0 1
                          0 .8 0.8
                          1 0 1
                          1 0.5 0.5
                          0.5 1 0.5
                          1 1 0.5
                          0.5 1 1                          
                          1 0.5 1
                          1 0.3 0.7
                          0.3 1 0.7
                          1 1 0.3
                          0.3 1 1                          
                          1 0.3 1];
                      
                if max(obj.e(5,:))>15
                    warning(['Number of boundary segments is to large,',...
                        'use option SegmentNumber to identify single ',...
                        'boundary segment.']);
                    return
                end
                for k = 1:length(obj.e(5,:))
                    line(obj.p(1,obj.e(1:2,k)),obj.p(2,obj.e(1:2,k)),...
                        'LineWidth',2,'color',colors(obj.e(5,k),:));
                    
                end   
                boundarysegments = unique(obj.e(5,:));
                for k =  boundarysegments    
                    annotation(...
                        'textbox',[0.85 .95-0.05*k 0.15 0.05],...
                        'String',['Segm. ',num2str(k)],...
                        'Color',colors(k,:),...
                        'LineStyle','none');
                end

            else
                clf
                obj.plot;
                indx = find(obj.e(5,:) == nSegment); 
                for k = 1:length(obj.e(5,:))
                    line(obj.p(1,obj.e(1:2,k)),obj.p(2,obj.e(1:2,k)),...
                        'LineWidth',1,'color','blue')
                end
                for k = 1:length(indx)
                    line(obj.p(1,obj.e(1:2,indx(k))),obj.p(2,obj.e(1:2,indx(k))),...
                        'LineWidth',3,'color','red')
                end
                if nSegment<=max(obj.e(5,:))
                    title(['Boundary segment no ',num2str(nSegment)]);
                else
                    cla
                    fprintf(['Sorry, there is not Boundary segment no ',...
                        num2str(nSegment),' in this geometry\n']);
                end
            end
        end
        
        function [cval,aval,fval] = aCoefficientsMpt(obj,c,a,f)
            %computes the value of the coefficients in the center of every triangle
            % 'symbolic' variables x, y  are neccesary for evaluation of string
            % objects like c = 'sin(x)' etc.
            % Restrictions
            % NOT jet implemented:
            % * cell array input
            % * full matrix c, c must be a 2 x 2 diagonal matrix.

            
           
            p = obj.p;  
            t = obj.t;

            switch class(c)
                case 'function_handle'
                    midp = obj.midpts;
                    x = midp(1,:);
                    y = midp(2,:);
                    cval = feval(c,x,y);
                    [rows,cols] = size(cval);
                    if rows == 1
                        % scalar
                        cval = [cval;cval];
                    end
                    if max(rows,cols)==obj.nPoints
                        cval = [obj.point2Center(cval(1,:));...
                            obj.point2Center(cval(2,:))];
                    end
                case 'double'
                    % four cases
                    % (i) scalar
                    % (ii) 2 x 2 matrix
                    % (iii) vector length np
                    % (iv) vector length nt
                    [rows,cols] = size(c);

                    switch max(rows,cols)
                        case 1
                            cval = c(ones(2,obj.nElements));
                        case 2
                            c1 = c(1,1); c2 = c(2,2);
                            cval(1,:) = c1(ones(1,obj.nElements));
                            cval(2,:) = c2(ones(1,obj.nElements));
                        case obj.nPoints
                            cval = obj.point2Center(c);
                            cval = [cval;cval];
                        case obj.nElements
                            % the good one, nothing to do
                            c = c(:);
                            cval = [c(:)';c(:)'];
                        otherwise
                             obj.wrongFormat.throwAsCaller
                    end
                case 'char'
                    % must be a single char symbolizing
                    % the coefficent function.
                    % For evaluating the coefficient function,
                    % we need x and y variables "hanging in the air".
                    % Hence, the next warning can be ignored!
                    midp = obj.midpts;
                    x = midp(1,:);
                    y = midp(2,:); %#ok<*NASGU>
                    try
                        cval = eval(c);
                    catch ME
                        ME.throwAsCaller;
                    end
                    [rows,cols] = size(cval);
                    switch max(rows,cols)
                        case 1
                            % c is a constant like 'pi'
                            cval = cval(ones(2,obj.nElements));
                        case obj.nElements
                            cval = [cval;cval];
                        otherwise
                            throw(obj.wrongFormat);
                    end
                case 'cell'
                    % must be a 2*2 cell array which refers to the entries
                    % of a 2*2
                    % diffussion matrix. Every Entry must be double or a single char symbolizing
                    % the part of the coefficent function like in the case
                    % 'char'.

                    [rows,cols] = size(c);
                    if (rows ~=2) || (cols ~=2)
                        throw(obj.wrongFormat);
                    end
                    cval = zeros(4,obj.nElements);
                    cvalcell = cell(2,2);


                    % For evaluating the coefficient function,
                    % we need x and y variables "hanging in the air".
                    % Hence, the next warning can be ignored!
                    midp = obj.midpts;
                    x = midp(1,:);
                    y = midp(2,:); %#ok<*NASGU>
                    try
                        for i = 1:2
                            for j = 1 : 2
                                switch class(c{i,j})
                                    case 'char'
                                        cvalcell{i,j} = eval(c{i,j});
                                    case 'double'
                                        cvalcell{i,j} = c{i,j};
                                    otherwise
                                        obj.wrongClass.throwAsCaller
                                end
                            end
                        end
                    catch ME
                        throw(ME);
                    end
                        for i = 1:2
                            for j = 1 : 2
                                linIndex = 2*(i-1)+j;
                                cvalij = cvalcell{i,j};
                                cvalij = cvalij(:);
                                [rows,cols] = size(cvalij);
                                switch max(rows,cols)
                                    case 1
                                        % c is a constant like 'pi'
                                        cval(linIndex,:) = cvalij(ones(1,obj.nElements));
                                    case  obj.nElements
                                        cval(linIndex,:) = cvalij';
                                    otherwise
                                        throw(obj.wrongFormat);
                                end
                            end
                        end
                otherwise
                    throw(obj.wrongClass)
            end

            % repeat code from case 'c'...
            switch class(a)
                case 'function_handle'
                    midp = obj.midpts;
                    x = midp(1,:);
                    y = midp(2,:);
                    aval = feval(a,x,y);
                    [rows,cols] = size(aval);
                    if max(rows,cols)==obj.nPoints
                        aval =  obj.point2Center(aval);
                        aval = aval(:)';
                    end
                case 'double'
                    [rows,cols] = size(a);
                    switch max(rows,cols)
                        case 1
                            aval = a(ones(1,obj.nElements));
                        case obj.nPoints
                            aval = obj.point2Center(a);
                            aval = aval(:)';
                            %
                        case obj.nElements
                            % the good one, nothing to do
                            aval = a(:)';
                        otherwise
                            obj.wrongFormat.throwAsCaller
                    end
                case 'char'
                    midp = obj.midpts;
                    x = midp(1,:);
                    y = midp(2,:);
                    try
                        aval = eval(a);
                    catch ME
                         ME.throwAsCaller;
                    end
                    [rows,cols] = size(aval);
                    switch max(rows,cols)
                        case 1
                            aval = aval(ones(1,obj.nElements));
                        case obj.nElements
                            % work already done
                        otherwise
                            finiteElements.wrongFormat.throwAsCaller;
                    end
                case 'cell'
                    obj.wrongClass.throwAsCaller;
                otherwise
                    obj.wrongClass.throwAsCaller;
            end

            switch class(f)
                case 'function_handle'
                    midp = obj.midpts;
                    x = midp(1,:);
                    y = midp(2,:);                     
                    fval = feval(f,x,y);
                    [rows,cols] = size(fval);
                    if max(rows,cols)==obj.nPoints
                        fval =  obj.point2Center(fval);
                    end
                case 'double'
                    [rows,cols] = size(f);
                    switch max(rows,cols)
                        case 1
                            fval = f(ones(1,obj.nElements));

                        case obj.nPoints
                            fval = obj.point2Center(f);

                        case obj.nElements
                            % the good one, nothing to do
                            fval = f;
                        otherwise
                            obj.wrongFormat.throw
                    end
%                     max(fval)
                case 'char'
                    midp = obj.midpts;
                    x = midp(1,:);
                    y = midp(2,:);
                    try
                        fval = eval(f);
                    catch ME
                        throw(ME);
                    end

                    [rows,cols] = size(fval);
                    switch max(rows,cols)
                        case 1
                           fval = fval(ones(1,obj.nElements));
                        case obj.nElements
                            % work already done
                        otherwise
                            throw(obj.wrongFormat);
                    end
                case 'cell'
                    error('Sorry, Cell array input not jet implemented!')
                otherwise
                    throw(obj.wrongClass)
            end
        end
        
        function [sidelength,area] = sideLengthAndArea(obj)
        end         
    end 
    methods(Static = true, Access = public)
        % here bcoefficient is STATIC
        function[qval,gval,hval,rval] = boundCoefficients(p,b,varargin)
            % compute the boundary coefficients
            %
            % q,g,h,r are everything that can be evaluated by eval of feval, The
            % independent variables must be named by x,y (Euklidian) or by s (arc-length
            % parametrization)
            % p must be grid.p(:,k) or grid.p
            % b must be one column of the boundary condition matrix.
            
            
            
            % Do not use gridd.x, p may be a subsequence of gridd.p
            if size(p,1) == 2
                x = p(1,:);
                y = p(2,:);
            else
                s = p;
            end
            if nargin==3
                t = varargin{1};
            end
            m = b(2);
            qval = 0;
            gval = 0;
            hval = 0;
            rval = 0;
            lengthq = b(3);
            lengthg = b(4);
            if m == 0 % only Neumann BCs
                try
                    qval = eval(char(b(5:5+lengthq-1)));
                    gval = eval(char(b(5+lengthq:5+lengthq+lengthg-1)));
                catch ME
                    ME.throw;
                end
            else % only Dirichlet BCs
                
                lengthh = b(5);
                lengthr = b(6);

                try
                    hval = eval(char(b(9:9+lengthh-1)));
                    if isscalar(hval)
                        try
                            hval = hval*ones(1,length(x));
                        catch
                            hval = hval*ones(1,length(s));
                        end
                    end
                    rval = eval(char(b(9+lengthh:9+lengthh+lengthr-1)));
                    if isscalar(rval)
                        try
                            rval = rval*ones(1,length(x));
                        catch
                            rval = rval*ones(1,length(s));
                        end
                    end
                catch ME
                    ME.throw;
                end
            end
        end
    end   
end

