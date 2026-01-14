classdef(Abstract) grid2D < gridd
    % grid2D class definition for 2D meshes
    % (c) 2025 Uwe Prüfert


    % NOTE grid2D is now abstract
    % All methods that shape geometry are protected.

    % Changes:
    % March 2024 redesign of abstract methods for  meshing.
    % All methods use now distmesh. freeGeometry is replaced by polygonGeometry
    % 2017/11/03 Bugfix in jiggleMesh: getMeshQuality(obj) => obj.meshQuality
    % Rename method  rectangle => rectangle.
    % Add identifySubdomain method.
    % 2016/08/31 Method gradient removed.
    % 2013/10/29
    % moveMesh, turnMesh methods added.
    % All pdetool dependent functionality removed.
    % Constructor: Only empty call allowed.
    % plot now adds a color bar when called with optional arguments.
    %
    % 2014/08/13 Fix a bug in identifyBoundarySegment when plotting
    % not ordered boundary segments.
    %
    % 2014/08/19 Fix some "strange behaviour" in plot. Use now low level
    % calls of patch.
    %
    % 2014/10/15 Fix the "convex set problem" in freeGeometry
    % U.P.
    % 2015/5/4 Add error estimator
    % Plot methods come now from plotUtils2d class


    % Define the abstract prop
    properties(SetAccess = protected,...
            GetAccess = public)
        nPointsInElements  = 3;
    end

    properties(Constant = true)
        spaceDimension = 2;
    end

    properties(Dependent = true,...
            SetAccess = protected)
        meshQuality
        hmax
        hmin
        hmean
        h
    end

    methods
        function val = get.meshQuality(obj)
            % angle = getInnerAngles(obj)

            % a,b,c are the rectangle of the edge lenghts
            a = sum((obj.p(:,obj.t(2,:))-obj.p(:,obj.t(1,:))).^2);
            b = sum((obj.p(:,obj.t(3,:))-obj.p(:,obj.t(2,:))).^2);
            c = sum((obj.p(:,obj.t(3,:))-obj.p(:,obj.t(1,:))).^2);

            val = sum([abs(acos((b+c-a)./(2*sqrt(b.*c)))-pi/3);...
                abs(acos((a+c-b)./(2*sqrt(a.*c)))-pi/3);...
                abs(acos((a+b-c)./(2*sqrt(a.*b)))-pi/3)]);
            val = 1 - 3*val/4/pi;
        end

        function val = get.hmax(obj)
            val = max(obj.triangleDiameters);
        end
        function val = get.hmin(obj)
            val = min(obj.triangleDiameters);
        end
        function val = get.hmean(obj)
            val = mean(obj.triangleDiameters);
        end
        function val = get.h(obj)
            val = obj.triangleDiameters;
        end
    end

    methods(Access = protected)
       % Basic geometries: free formed boundary by polygon, reactangle
       % and ellipse.

        function polygon(obj,polygon,hmax)
        % POLYGONG
        % Now using distmesh.
        % Note that ever segment of the polygon defines also a boundary segment
        % in the mesh. If you want to change this, redefine the nimbering in the
        % derived class.
        % DO NOT using this function to create circles etc. You will get a
        % huge number of boundry segments.

        % Error management
            switch nargin
                case 1
                    MException('GRID2D:POLYGONGEOMETRY:EMPTYCALL',...
                    ['Empty call is not allowed.',...
                    ' At least the polygomial must be given']).throwAsCaller;
                case 2
                    hmax = 0.5;
                case 3
                %
                otherwise
            end

            % polygon must be in the form x y, where x and y are column vectors.
            [m,n] = size(polygon);
            if m~=2 && n~=2
                MException('GRID2D:POLYGON:WRONGFORMAT',...
                ['Polygonial must be given as',...
               ' an n x 2 or 2 x n matrix']).throwAsCaller;
            end
            if max(m,n)<3
                MException('GRID2D:POLYGONGEOMETRY:NOTENOUGHPOINTS',...
                ['The polynon must contain'...
                 ' at least three points lying not on a line']).throwAsCaller;
            end
            if m == 2
                polygon = polygon';
            end
            % Call distmesh.
            % ho is now defined such that at leas one point is created in the domain.
            h0 = min(max(polygon(:,1))-min(polygon(:,1)),...
                  max(polygon(:,2))-min(polygon(:,2)))/2;

            [p,e,t] = distmesh2d(@dpoly,@huniform,...
                    h0,...
                    [min(polygon(:,1)) min(polygon(:,2))
                     max(polygon(:,1)) max(polygon(:,2))],...
                     polygon,polygon);

            % Bring p-e-t in OOPDE form
            obj.p = p';
            obj.e = e';
            obj.t = t';
            obj.t(4,:) = ones(1,size(obj.t,2));



            % Correct boundary segment numbering
            % We will have for every segment in the polygonial also a
            % boundary segment in the mesh.
            [m,n] = size(polygon);
            % force polygon to be a row matrix.
            if m > n
                polygon = polygon';
            end
            % Identify corner points from polygon
            cornerPoints =  zeros(1,size(polygon,2));
            for k = 1:max(size(polygon))
                n1 = obj.p(1,:)==polygon(1,k);
                n2 = obj.p(2,:)==polygon(2,k);
                cornerPoints(k) = find(n1.*n2 == 1);
            end
            % Keep number of corner points and extend cornerpoint by the
            % first one to wrap arround.
            NCP = length(cornerPoints);
            cornerPoints = [cornerPoints cornerPoints(1)];
            % We start with boundary segment 1
            BS = 1;
            % For number of segments in polygon we look for the elements
            % belonging to one segment
            for k1 = 1:NCP
                akt = cornerPoints(k1);
                while cornerPoints(k1+1)~=akt
                    k2 = find(akt==obj.e(1,:));
                    obj.e(end,k2) = BS;
                    akt = obj.e(2,k2);
                end
                BS = BS + 1;
            end

            % If hmax is given, refine uniformly
            h = obj.hmax;
            % Some messages...
            if h > hmax
                fprintf('Refine mesh.')
                toBeRefine = true;
            else
                toBeRefine = false;
            end
            while h > hmax
                obj.refineMesh();
                fprintf('.');
                h = obj.hmax;
            end
            if toBeRefine
                fprintf('done.\n');
            end

            % Try to jiggle the mesh for better quality.
            % TODO: Trigger this only if mesh quality is too bad
            fprintf('Jiggle mesh...');
            for k = 1:3
                obj.jiggleMesh;
            end
            fprintf('done.\n');
        end


        function ellipse(obj,a1,a2,a3,a4,a5)
            switch nargin
                case 1
                    A = 2;  B = 1;  hmax = 0.2; a = 0;  b = 0;
                case 2
                    A = 2;  B = 1;  hmax = a1;  a = 0;  b = 0;
                case 3
                    A = a1; B = a2; hmax = 0.2; a = 0;  b = 0;
                case 4
                    A = a1; B = a2; hmax = a3;  a = 0;  b = 0;
                case 5
                    A = a1; B = a2; hmax = 0.25;  a = a3; b = a4;
                case 6
                    A = a1; B = a2; a = a3; b = a4; hmax = a5;
                otherwise
                    obj.wrongNumberInputs.throw;
            end
            if A<=0 || B <= 0
                MException('ELLIPSE:wrongDefinedInput',...
                    'Error: A and B must be greater then zero.').throwAsCaller;
            end
            % Prepare the boundary description function
            fd=@(p) (p(:,1)-a).^2/A^2+(p(:,2)-b).^2/B^2-1;
            % Call distmesh
            [p,e,t] = distmesh2d(fd,@huniform,hmax,[-A+a,-B+b
                                                     A+a,B+b],[]);
            obj.p = p';
            obj.e = e';
            obj.t = t';
            obj.t(4,:) = ones(1,size(obj.t,2));
        end

        function rectangle(obj,varargin)
            % rectangle()
            % rectangle(h)
            % rectangle(xmin,xmax,ymin,ymax)
            % rectangle(xmin,xmax,ymin,ymax,h)
            % No options => unitsqaure with h = 0.1
            % hmax given but but no corner points => unitsquare with given
            % left lower and right upper point given,
            % but no hmax => rectangle with hmax = 0.1
            % The boundary lenght parameter runs in anti clockwise
            % order from 0 to 1 on everey boundary segment,
            % independetly from its length!
            % (c) 2013 by Uwe Prüfert



            % First handle all posible argument...
            switch nargin
                case 1
                    % no options,unitsqaure with h = 0.1
                    xmin = 0;
                    xmax = 1;
                    ymin = 0;
                    ymax = 1;
                    hmax = 0.1;
                case 2
                    % hmax but no coordinates, unitsquare with given
                    % hmax
                    xmin = 0;
                    xmax = 1;
                    ymin = 0;
                    ymax = 1;
                    hmax = varargin{1};
                case 5
                    % left lower and right upper point given,
                    % but no hmax, rectangle with hmax = 0.1
                    xmin = varargin{1};
                    xmax = varargin{2};
                    ymin = varargin{3};
                    ymax = varargin{4};
                    hmax = 0.1;
                case 6
                    xmin = varargin{1};
                    xmax = varargin{2};
                    ymin = varargin{3};
                    ymax = varargin{4};
                    hmax = varargin{5};
                otherwise
                     obj.wrongNumberInputs.throw;
            end

            fd=@(p) -min(min(min(-ymin+p(:,2),ymax-p(:,2)),-xmin+p(:,1)),xmax-p(:,1));
            [p,e,t]= distmesh2d(fd,@huniform,hmax,[xmin,ymin;...
                                                   xmax,ymax],[xmin,ymin;...
                                                               xmax,ymin;...
                                                               xmin,ymax;...
                                                               xmax,ymax]);

            obj.p = p';
            obj.e = e';
            obj.t = t';
            obj.t(4,:) = ones(1,size(obj.t,2));
        end

        % Helper methods for refinement

        function neighbourIndex = getNeighbourPointsIndex(obj,k,d)
            %neighbourIndex = obj.getNeighbourPointsIndex(k)
            % Computes the index of all neighboured points of
            % point k.
            %neighbourIndex = obj.getNeighbourPointsIndex(k,d)
            % Computes the index of all neighboured points of
            % point k in Domain d.

            [~,indx] = find(obj.t(1:3,:)==k);
            if (nargin == 3)
                indx = indx(obj.t(4,indx) == d);
            end
            neighbourIndex = unique(obj.t(1:3,indx))';
            neighbourIndex(neighbourIndex==k) = [];
        end

        function  refineRGB(obj,markedElements)
            %refineRGB: local refinement of finite element mesh by red-green-blue
            %           refinement, where marked elements are red-refined.
            %
            %Usage:
            %
            %obj.refineRGB(markedElements)
            %
            %    This methods bases on a code taken from the paper
            %    >> Efficient Implementation of Adaptive P1-FEM in Matlab <<
            %    by S. Funken, D. Praetorius, and P. Wissgott.
            % Original Authors:
            %
            %      10-07-08
            %  Adaption for grid2D class by Uwe Pruefert   2013
            %
            % ChangeLog
            % 03/13/2015
            % the function did not preserve Domain numbers. This is not
            % fixed


            % if object is extended, we only need the old grid structure
            if obj.isExtended
                coordinates = obj.p(:,1:obj.NrPO)';
                boundary    = obj.e(1:end-1,:)';
                elements    = [obj.t(1:3,:);obj.t(end,:)]';
            else
                coordinates = obj.p';
                boundary    = obj.e';
                elements    = obj.t';
            end



            nE = size(elements,1);
            % Sort elements such that first edge is longest
            dx = coordinates(elements(:,[2,3,1]),1)-coordinates(elements(:,1:3),1);
            dy = coordinates(elements(:,[2,3,1]),2)-coordinates(elements(:,1:3),2);
            [~,idxMax] = max(reshape(dx.^2+dy.^2,nE,3),[],2);
            idx = ( idxMax==2 );

            % Add 4th column for subdomain information
            elements(idx,:) = elements(idx,[2,3,1,4]);
            idx = ( idxMax==3 );
            elements(idx,:) = elements(idx,[3,1,2,4]);
            % Obtain geometric information on edges
            % Restrict to 1:3, ignore 4th column
            [edge2nodes,element2edges,boundary2edges{1}] ...
                = provideGeometricData(elements(:,1:3),boundary(:,1:2));
            % Mark edges for refinement
            edge2newNode = zeros(max(max(element2edges)),1);
            edge2newNode(element2edges(markedElements,:)) = 1;
            swap = 1;
            while ~isempty(swap)
                markedEdge = edge2newNode(element2edges);
                swap = find( ~markedEdge(:,1) & (markedEdge(:,2) | markedEdge(:,3)) );
                edge2newNode(element2edges(swap,1)) = 1;
            end
            % Generate new nodes

            edge2newNode((edge2newNode)>0) = size(coordinates,1) + ...
                (1:nnz(edge2newNode));
            idx = find(edge2newNode);
            coordinates(edge2newNode(idx),:) ...
                = (coordinates(edge2nodes(idx,1),:)+...
                coordinates(edge2nodes(idx,2),:))/2;

            %  Refine boundary elements
            % In constrast to the original code,
            % we must extend the boundary matrix
            % by adding informations on the number
            % of boundary and the start and end
            % parameter. Importand for evaluating
            % u = f(s) statements... U.P.
            if ~isempty(boundary)
                newNodes = edge2newNode(boundary2edges{1});
                markedEdges = find(newNodes);
                % We compute the   midpoint of
                % edge. It will be the end point
                % of the first end the starting
                % point of the second set of new
                % edges. U.P.
                newEdgeMidpoint = (boundary(markedEdges,4)+boundary(markedEdges,3))/2;

                if ~isempty(markedEdges)
                    boundary = [boundary(~newNodes,:); ...
                                boundary(markedEdges,1),newNodes(markedEdges),...
                                boundary(markedEdges,3),newEdgeMidpoint,boundary(markedEdges,5:end); ...
                                newNodes(markedEdges),boundary(markedEdges,2),...
                                newEdgeMidpoint,boundary(markedEdges,4:end)];
                end
            end
            newNodes = edge2newNode(element2edges);
            %  Determine type of refinement for each element
            markedEdges = (newNodes~=0);
            none = ~markedEdges(:,1);
            bisec1  = ( markedEdges(:,1) & ~markedEdges(:,2) & ~markedEdges(:,3) );
            bisec12 = ( markedEdges(:,1) &  markedEdges(:,2) & ~markedEdges(:,3) );
            bisec13 = ( markedEdges(:,1) & ~markedEdges(:,2) &  markedEdges(:,3) );
            red     = ( markedEdges(:,1) &  markedEdges(:,2) &  markedEdges(:,3) );
            %  Generate element numbering for refined mesh
            idx = ones(nE,1);
            idx(bisec1)  = 2; %*** green = newest vertex bisection of 1st edge
            idx(bisec12) = 3; %*** blue (right) = newest vertex bisection of 1st and 2nd edge
            idx(bisec13) = 3; %*** blue (left) = newest vertex bisection of 1st and 3rd edge
            idx(red)     = 4; %*** red refinement
            idx = [1;1+cumsum(idx)];

            % Generate new elements
            newElements = zeros(idx(end)-1,4);
            newElements(idx(none),:) = elements(none,:);

            % add iformation on subdomain number
            newElements([idx(bisec1),1+idx(bisec1)],:) ...
                = [elements(bisec1,3),elements(bisec1,1),newNodes(bisec1,1),elements(bisec1,4); ...
                   elements(bisec1,2),elements(bisec1,3),newNodes(bisec1,1),elements(bisec1,4)];

            newElements([idx(bisec12),1+idx(bisec12),2+idx(bisec12)],:) ...
                = [elements(bisec12,3),elements(bisec12,1),newNodes(bisec12,1),elements(bisec12,4); ...
                   newNodes(bisec12,1),elements(bisec12,2),newNodes(bisec12,2),elements(bisec12,4); ...
                   elements(bisec12,3),newNodes(bisec12,1),newNodes(bisec12,2),elements(bisec12,4)];

            newElements([idx(bisec13),1+idx(bisec13),2+idx(bisec13)],:) ...
                = [newNodes(bisec13,1),elements(bisec13,3),newNodes(bisec13,3),elements(bisec13,4); ...
                   elements(bisec13,1),newNodes(bisec13,1),newNodes(bisec13,3),elements(bisec13,4); ...
                   elements(bisec13,2),elements(bisec13,3),newNodes(bisec13,1),elements(bisec13,4)];

            newElements([idx(red),1+idx(red),2+idx(red),3+idx(red)],:) ...
                = [elements(red,1),newNodes(red,1),newNodes(red,3),elements(red,4); ...
                   newNodes(red,1),elements(red,2),newNodes(red,2),elements(red,4); ...
                   newNodes(red,3),newNodes(red,2),elements(red,3),elements(red,4); ...
                   newNodes(red,2),newNodes(red,3),newNodes(red,1),elements(red,4)];

            if obj.isExtended
                obj.isExtended = false;
                obj.p = coordinates';
                obj.e = boundary';
                obj.t = newElements';
                extendMesh(obj);

            else
                obj.p = coordinates';
                obj.e = boundary';
                obj.t = newElements';
            end


            % local function
            function [edge2nodes,element2edges,boundaries] ...
                    = provideGeometricData(elements,boundaryG)
                %provideGeometricData: returns geometric data for finite element mesh
                %
                %Usage:
                %
                %[edges2nodes,element2edges,dirichlet2edges,neumann2edges] ...
                %    = provideGeometricData(elements,dirichlet,neumann)
                %
                %Authors:
                %
                %    S. Funken, D. Praetorius, P. Wissgott  10-07-08

                nEE = size(elements,1);
                nB = nargin-1;
                %*** Node vectors of all edges (interior edges appear twice)
                I = elements(:);
                J = reshape(elements(:,[2,3,1]),3*nEE,1);
                %*** Symmetrize I and J (so far boundary edges appear only once)
                pointer = [1,3*nEE,zeros(1,nB)];


                I = [I;boundaryG(:,2)];
                J = [J;boundaryG(:,1)];

                pointer(3) = pointer(2) + size(boundaryG,1);

                %*** Create numbering of edges
                idxIJ = find(I < J);
                edgeNumber = zeros(length(I),1);
                edgeNumber(idxIJ) = 1:length(idxIJ);
                idxJI = find(I > J);
                number2edges = sparse(I(idxIJ),J(idxIJ),1:length(idxIJ));
                [foo{1:2},numberingIJ] = find(number2edges);
                [foo{1:2},idxJI2IJ] = find(sparse(J(idxJI),I(idxJI),idxJI) ); %#ok<ASGLU>

                edgeNumber(idxJI2IJ) = numberingIJ;
                % Provide element2edges and edge2nodes
                element2edges = reshape(edgeNumber(1:3*nEE),nEE,3);
                edge2nodes = [I(idxIJ),J(idxIJ)];
                % Provide boundary2edges

                boundaries = edgeNumber(pointer(2)+1:pointer(3));

            end
        end
    end

    methods(Static,Access = public)
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
                    qval = eval(char(b(5:5+lengthq-1)'));
                    gval = eval(char(b(5+lengthq:5+lengthq+lengthg-1)'));
                catch ME
                    MException(ME.identifier,ME.message).throwAsCaller;
                end
            else % only Dirichlet BCs

                lengthh = b(5);
                lengthr = b(6);

                try
                    hval = eval(char(b(9:9+lengthh-1)'));
                    if isscalar(hval)
                        try
                            hval = hval*ones(1,length(x));
                        catch
                            hval = hval*ones(1,length(s));
                        end
                    end
                    rval = eval(char(b(9+lengthh:9+lengthh+lengthr-1)'));
                    if isscalar(rval)
                        try
                            rval = rval*ones(1,length(x));
                        catch
                            rval = rval*ones(1,length(s));
                        end
                    end
                catch ME
                     MException(ME.identifier,ME.message).throwAsCaller;
                end
            end
        end
    end

    methods(Access = public)

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


            % New handling of inputs. Additional argument axes to plot into.
            % First two entries in varargin could be y and ax, last are optionsal arguments for patch


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



            % Plot the  mesh triangles, no solution vector given.
            if (nargin == 1) || (mod(length(varargin),2) == 0)
                x = obj.p(1,:);
                y = obj.p(2,:);
                z = zeros(size(x));
                patch(...
                    'faces',obj.t(1:3,:)',...
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
                    x = reshape(obj.p(1,obj.t(1:3,:)),3,obj.nElements);
                    y = reshape(obj.p(2,obj.t(1:3,:)),3,obj.nElements);
                    z = [1;1;1]*z(:)';
                    cz = z;
                    patch(x,y,z,cz,...
                        'lineStyle','none',...
                        varargin{2:end});
                    view(ax,3);
                    grid(ax,'on');
                    colorbar;
                % P1 or p2 plot
                elseif length(z) == obj.nPoints
                    % P2 plot
                    if obj.isExtended
                        patch('faces',obj.t([1,4,6],:)','vertices',[x(:),y(:),z(:)],...
                            'facevertexcdata',z(:),...
                            'facecolor','interp',...
                            'edgecolor',get(ax,'DefaultSurfaceEdgeColor'),...
                            'parent',ax,...
                            varargin{2:end});
                        patch('faces',obj.t([4,2,5],:)','vertices',[x(:),y(:),z(:)],...
                            'facevertexcdata',z(:),...
                            'facecolor','interp',...
                            'edgecolor',get(ax,'DefaultSurfaceEdgeColor'),...
                            'parent',ax,...
                            varargin{2:end});
                        patch('faces',obj.t([5,6,4],:)','vertices',[x(:),y(:),z(:)],...
                            'facevertexcdata',z(:),...
                            'facecolor','interp',...
                            'edgecolor',get(ax,'DefaultSurfaceEdgeColor'),...
                            'parent',ax,...
                            varargin{2:end});
                        patch('faces',obj.t([5,3,6],:)','vertices',[x(:),y(:),z(:)],...
                            'facevertexcdata',z(:),...
                            'facecolor','interp',...
                            'edgecolor',get(ax,'DefaultSurfaceEdgeColor'),...
                            'parent',ax,...
                            varargin{2:end});
                        view(ax,3);
                        grid(ax,'on');
                        colorbar;
                    else
                    %  P1 plot
                       patch('faces',obj.t(1:3,:)','vertices',[x(:),y(:),z(:)],...
                            'facevertexcdata',z(:),...
                            'facecolor','interp',...
                            'edgecolor',get(ax,'DefaultSurfaceEdgeColor'),...
                            'parent',ax,...
                            varargin{2:end});
                        view(ax,3);
                        grid(ax,'on');
                        colorbar;
                    end
                else
                    MException('GRID2D:PLOT:WRONGELEMENT',...
                        'Elements can only be P0,P1 or P2.').throwAsCaller;
                end
            end




            axis equal
         end

        function identifyBoundarySegment(obj,nSegment)
            %%
            %
            % * identifyBoundarySegment  IN:self[,double] OUT:none
            %
            % Method that plot the grid and makrs the boundary segments by colors.
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

            % (c) 2013 by Uwe Prüfert
            % Change log
            % 2014/08/13 Fix a bug in identifyBoundarySegment when plotting
            % not ordered boundary segments
            % U.P.

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
                        'String',['Segment ',num2str(k)],...
                        'Color',colors(k,:),...
                        'LineStyle','none');
                end

            else
                clf
                obj.plot;
                indx = [];
                for k = 1:length(nSegment)

                    indx = [indx,find(obj.e(5,:) == nSegment(k))];
                end
    %             line(obj.p(1,obj.e(1:2,:)),obj.p(2,obj.e(1:2,:)));
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

        function identifySubdomain(obj,sub)
            switch nargin
                case 1
                    obj.plot(obj.t(end,:));
                    colorbar off;
                    view(2);
                case 2
                    subs = unique(obj.t(end,:));
                    if sub<=max(subs)
                        obj.plot(double(obj.t(end,:)==sub));
                        view(2);colorbar off;
                    else
                        fprintf('Subdonain with number %d do not exists.\n',sub);
                    end
                otherwise
            end
        end

        function quiver(obj,dx,dy)
            %%
            %
            % * quiver IN:self,double,double OUT:none
            %
            % quiver plot method for grid2D.
            % Plots vector field. The arguments are double vectors dx and dy
            % of length #points or #elements.
            %
            % Example 1:
            %
            %       dfdx = @(x,y) 2*x;
            %       dfdy = @(x,y) pi*cos(pi*y);
            %       grid.quiver(dfdx(grid.x,grid.y),dfdy(grid.x,grid.y))
            %
            % evaluates gradient at grid points.
            %
            % Note that the gradient of a
            % finite element solutions is only defined in the interior of
            % the elements but not in the grid nodes or on the edges.
            %
            % Example 2:
            %
            %       pde.solve();
            %       [dx,dy] = pde.fem.gradient(pde.grid,pde.y);
            %       pde.grid.quiver(dx,dy)
            %
            % computes the gradient of a solution and plots its
            % vector-field.

            % Change log
            % Code revisited 3/3/2016

            obj.plot('EdgeColor',[0.75 0.75 0.7])
            hold on
            if length(dx) == obj.nPoints && length(dy) == obj.nPoints
                quiver(obj.x,obj.y,dx,dy);
            elseif length(dx) == obj.nElements && length(dy) == obj.nElements
                midpoints = obj.midpts;
                quiver(midpoints(1,:),midpoints(2,:),dx,dy);
            else
                grid2D.wrongFormat.throwAsCaller;
            end
            view(2);
            axis equal
            hold off
        end

        function extendMesh(obj)
        % extends the given mesh for P-2 elements
        % call with: yourobjectname.grid.extendMesh()
        %  Written by Marcel Heim for his Diploma thesis

            % if the mesh is already extended, return to calling function
            if obj.isExtended
                return
            end
            % get the number of old points before the mesh is extended
            obj.NrPO = size(obj.p,2);
            % get the number of triangles
            [~,NrTri] = size(obj.t);
            % get the number of edges
            [~,NrEdg] = size(obj.e);
            % prepare some index vectors with length NrTri and NrEdg
            idx = 1:NrTri;
            idx1 = 1:NrEdg;
            % get midpoints of each triangle edge
            xm1 = 0.5*(obj.p(1,obj.t(1,idx))+obj.p(1,obj.t(2,idx)));
            ym1 = 0.5*(obj.p(2,obj.t(1,idx))+obj.p(2,obj.t(2,idx)));
            xm2 = 0.5*(obj.p(1,obj.t(2,idx))+obj.p(1,obj.t(3,idx)));
            ym2 = 0.5*(obj.p(2,obj.t(2,idx))+obj.p(2,obj.t(3,idx)));
            xm3 = 0.5*(obj.p(1,obj.t(1,idx))+obj.p(1,obj.t(3,idx)));
            ym3 = 0.5*(obj.p(2,obj.t(1,idx))+obj.p(2,obj.t(3,idx)));
            % get the midpoint of each mesh edge of obj.e
            xme = 0.5*(obj.p(1,obj.e(1,idx1))+obj.p(1,obj.e(2,idx1)));
            yme = 0.5*(obj.p(2,obj.e(1,idx1))+obj.p(2,obj.e(2,idx1)));
            % sort the triangle values like this:
            % (xm1,ym1) of triangle 1
            %    ...
            % (xm3,ym3) of triangle 1
            % (xm1,ym1) of triangle 2 and so on
            tobesorted = [xm1;ym1;xm2;ym2;xm3;ym3];
            sorted = reshape(tobesorted(:),2,3*NrTri)';
            % get the points we have to add to the mesh
            % each interior triangulation mesh edge gets counted twice and
            % thus has to be omitted by unique()
            % idx2 provides the data which row of 'sorted' maps to which row
            % of idx2
            [newPoints,~,idx2] = unique(sorted,'rows');
            obj.p = [obj.p newPoints'];
            % reshape the idx2 vector into a 3 x NrTri matrix and add
            % the number of old points to the values ( compare with
            % new point structure above: [[1...NrPO] [newpoints]]
            obj.t = [obj.t(1:3,:);reshape(idx2+obj.NrPO,3,NrTri);obj.t(4,:)];
            idxMPe = zeros(1, NrEdg);
            for idx1=1:NrEdg
                 %   idxMPe(idx1) = findpoint(newPoints,[xme(idx1) yme(idx1)]);
                 % MEX File removed. Thanks to Hannes Uecker.
                [~,id] = ismember([xme(idx1) yme(idx1)],newPoints,'rows');
                idxMPe(idx1) = id;
            end
            %finally add the new edge indices as last row to obj.e
            obj.e = [obj.e;idxMPe+obj.NrPO];
            % object is now extended hence we set the value to true
            obj.isExtended = true;
            % and we get 6 points per triangle
            obj.nPointsInElements = 6;
        end
        %heim_neu

        function refineMesh(obj,varargin)
            % Method for refinement of grid
            % obj.refineMesh();
            % obj.refineMesh('triangle',triangles_to_refine);
            % obj.refineMesh('domain',domains_to_refine);
            % obj.refineMesh('locationAllPoints',restrictions_for_elements_based_on_location);
            % restrictions_for_elements_based_on_location have to be
            % satisfied for all points of a triangle which shall be refined
            % eg obj.refineMesh('locationAllPoints','x <= 5 & y >7');
            % obj.refineMesh('locationOnePoint',restrictions_for_elements_based_on_location);
            % restrictions_for_elements_based_on_location have to be
            % satisfied for one point of a triangle which shall be refined
            % and for compatibility with old versions (maybe think of romove it)
            % obj.refineMesh(triangles_to_refine);
            %
            % ChangeLog
            % 03/13/2015
            % added the functionallity corresponding to
            % obj.refineMesh('triangle',triangles_to_refine);
            % obj.refineMesh('domain',domains_to_refine);
            % obj.refineMesh('locationAllPoints',restrictions_for_elements_based_on_location);
            % obj.refineMesh('locationOnePoint',restrictions_for_elements_based_on_location);



            switch nargin
                case 1
                    toRefine = 1:obj.nElements;
                case 2
                    toRefine = varargin{1};
                    if ~(isvector(toRefine)&&length(toRefine)<=obj.nElements)
                        % error if more elements are in refinement list as in
                        % the mesh
                        obj.wrongFormat.throw;
                    end
                case 3
                    switch varargin{1}
                        case {'triangle', 'triangles'}
                            toRefine = varargin{2};
                            if ~(isvector(toRefine)&&length(toRefine)<=obj.nElements)
                                % error if more elements are in refinement list as in
                                % the mesh
                                obj.wrongFormat.throw;
                            end
                        case {'domain','domains'}
                            domainsToRefine = varargin{2};
                            %heim_neu
                            %toRefine = find(ismember(obj.t(4,:) , domainsToRefine));
                            toRefine = find(ismember(obj.t(7,:) , domainsToRefine));
                            %heim_neu
                            if isempty(toRefine)
                                warning(['No Triangles to Refine. ',...
                                    'Maybe the wrong domainnumber was specified. Check it!']);
                            end
                        case {'locationAllPoints','locationsAllPoints'}
                            % Stategy: Find for every triangel the points
                            % which satisfy the condition varargin{2}.
                            % Refine all triangles which have 3 Points
                            % satisfying the condition
                            condition = varargin{2};
                            satisfied = zeros(3,obj.nElements);
                            try
                                for ip = 1:3
                                    x = obj.p(1,obj.t(ip,:));
                                    y = obj.p(2,obj.t(ip,:));
                                    satisfied(ip,:) = eval(condition);
                                end
                            catch ME
                                throw(ME);
                            end
                            toRefine = find(sum(satisfied) == 3);
                        case {'locationOnePoint','locationsOnePoint'}
                            % Stategy: Find for every triangel the points
                            % which satisfy the condition varargin{2}.
                            % Refine all triangles which have 1 or more Points
                            % satisfying the condition
                            condition = varargin{2};
                            satisfied = zeros(3,obj.nElements);
                            try
                                for ip = 1:3
                                    x = obj.p(1,obj.t(ip,:));
                                    y = obj.p(2,obj.t(ip,:));
                                    satisfied(ip,:) = eval(condition);
                                end
                            catch ME
                                throw(ME);
                            end
                            toRefine = find(sum(satisfied) >= 1);
                        otherwise
                            wrongFormat.throw;
                    end
                otherwise
                    wrongNumberInputs.throw;
            end
            try
                obj.refineRGB(toRefine);
            catch ME

                MException('GRDI2D:REFINEMESH:INVALIDGRID',...
                    ['Error when refine the mesh. Check the structure',...
                    ' of your mesh. Message was: ',...
                    ME.message]).throwAsCaller;
            end

        end

        function jiggleMesh(obj,prop,value)
            switch nargin
                case 1
                    obj.jiggle;
                case 3
                    switch prop
                        case 'n'
                            for k = 1:value
                                obj.jiggle;
                            end
                        case 'quality'
                            k = 0;
                            maxjiggle = 15;
                            while min(obj.meshQuality)<value...
                                    &&k<maxjiggle
                               obj.jiggle;
                               k = k + 1;
                            end
                        otherwise
                            obj.wrongFormat.throwAsCaller;
                    end
                otherwise
                    obj.wrongNumberInputs.throw;
            end
        end

        % utilities

        function innerPoints = getInnerPoints(obj)
            %innerPoints = obj.getInnerPoints
            % Gives back the  points not on the boundary of a grid.
            % (c) 2013 by Uwe Prüfert
            indx = obj.getInnerPointsIndex;
            if ~isempty(indx)
                innerPoints = obj.p(:,indx);
            else
                warning('grid2D:emptyInnerPointsIndex',...
                    'Grid contains no inner points.');
            end
        end

        % coefficients

        function [bvalvec] = convCoefficientsMpt(obj,b)
            %computes the value of the coefficient b in the center of every triangle
            % coefficent can be a vector of dim 2 x 1, or a cell array of dim 2 x 1

            % 'symbolic' variables x, y  and t are neccesary for evaluation of string
            % objects like c = 'sin(x)' etc.

            if obj.isExtended
                p = obj.p(:,1:obj.NrPO);
            else
                p = obj.p;
            end

            t = obj.t;

            x = obj.x;
            x = obj.point2Center(x);
            y = obj.y;
            y = obj.point2Center(y);



            if ~((max(size(b))>=2)&&(min(size(b))>=1))
                MException('ccoefficients:wrongCoefficientDefinition',...
                    ' b must be a vector').throwAsCaller;

            end
            % two cases:
            % cell-array - entries are strings, or doubles
            if isa(b,'cell')
                b1 = b{1};
                b2 = b{2};

            elseif isa(b,'double')
                b1 = b(1,:);
                b2 = b(2,:);
            else
                MException('ccoefficients:wrongCoefficientDefinition',...
                    ' Wrong coefficient definition').throwAsCaller;

            end

            if isa(b1,'function_handle') || isa(b1,'inline')
                bval = feval(b1,x,y);
            elseif isa(b1,'char')
                % b1 is a scalar or a row vector of length nElements, see
                % definition of x and y
                bval = eval(b1).*ones(1,obj.nElements);
            elseif isa(b1,'numeric')
                if length(b1)==obj.nPoints
                    % c vektor and defined in p
                    bval = b1;
                elseif length(b1)==1
                    % skalar
                    bval = b1*ones(obj.nElements,1);
                elseif length(b1)==obj.nElements
                    bval = b1;
                else
                    MException('ccoefficients:wrongSize',...
                        'wrong sized b(1)').throwAsCaller;
                end
            else
                MException('ccoefficients:wrongSize',...
                    'wrong formated b(1)').throwAsCaller;
            end
            if length(bval) == obj.nElements
                % b(1) is a vektor and defined in center of mass of triagle
            else
                dimb = size(bval);
                if dimb(1) == 1
                    bval = bval';
                end
                bval = obj.point2Center(bval);
            end
            dimb = size(bval);
            if dimb(1) == 1
                bval = bval';
            end
            % first column of bvalvec
            bvalvec = bval;

            if isa(b2,'function_handle') || isa(b2,'inline')
                bval = feval(b2,x,y);
            elseif isa(b2,'char')
                bval = eval(b2).*ones(1,obj.nElements);
            elseif isa(b2,'numeric')
                if length(b2)==obj.nPoints
                    % c vektor and defined in p
                    bval = b2;
                elseif length(b2)==1
                    % skalar
                    bval = b2*ones(obj.nElements,1);
                elseif length(b2)==obj.nElements
                    bval = b2;
                else
                    Exception('ccoefficients:wrongSize',...
                        'wrong sized b(1)').throwAsCaller;

                end
            else
                MException('ccoefficients:wrongSize',...
                    'wrong formated b(1)').throwAsCaller;

            end

            if length(bval) == obj.nElements
                % b(1) is a vektor and defined in center of mass of triagle

            else
                dimb = size(bval);
                if dimb(1) == 1
                    bval = bval';
                end
                bval = obj.point2Center(bval);
            end
            dimb = size(bval);
            if dimb(1) == 1
                bval = bval';
            end
            bvalvec = [bvalvec,bval];
        end

        function [cval,aval,fval] = aCoefficientsMpt(obj,c,a,f)
            %computes the value of the coefficients in the center of every triangle
            % 'symbolic' variables x, y  are neccesary for evaluation of string
            % objects like c = 'sin(x)' etc.
            % Restrictions
            % NOT jet implemented:
            % * cell array input
            % * full matrix c, c must be a 2 x 2 diagonal matrix.

            %heim_neu
            if obj.isExtended
                p = obj.p(:,1:obj.NrPO);
            else
                p = obj.p;
            end
            %heim_neu
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
                        cval = [obj.point2Center(cval(1,:))
                            zeros(1,obj.nElements)
                            zeros(1,obj.nElements)
                            obj.point2Center(cval(2,:))];
                    end
                case 'double'
                    [rows,cols] = size(c);
                    %
                    % Three different cases:
                    % (i) scalar,
                    % (ii) diagonal matrix,
                    % (iii) full matrix.

                    % (ia) c is a scalar constant
                    if rows == 1 && cols == 1
                        cval = [c(ones(1,obj.nElements))
                            zeros(2,obj.nElements)
                            c(ones(1,obj.nElements))];

                    % (ib) c is a row vector if length nElements
                    elseif rows == 1 && cols == obj.nElements
                        cval = [c
                            zeros(2,obj.nElements)
                            c];

                    % (ic) c is a row vector if length nPoints
                    elseif rows == 1 && cols == obj.nPoints
                        cval = [obj.point2Center(c)
                            zeros(2,obj.nElements)
                            obj.point2Center(c)];

                    % (iia) c is a row vector if length 2*nElements
                    % It will be handled as a 2 x 2 diagonal matrix where
                    % c(1,1) and c(2,2) contain the value of a space
                    % dependent coefficient function at element's centers
                    % (and c(1,2) = c(2,1) = 0).
                    elseif rows == 1 && cols == 2*obj.nElements
                        cval = [c(1:obj.nElements)
                            zeros(2,obj.nElements)
                            c(1+obj.nElements:2*obj.nElements)];

                    % (iib) c is a row vector if length 2*nPoints
                    % It will be handled as a 2 x 2 diagonal matrix where
                    % c(1,1) and c(2,2) contain the value of a space
                    % dependent coefficient function at grid's nodes
                    % (and c(1,2) = c(2,1) = 0).
                    elseif rows == 1 && cols == 2*obj.nPoints
                        cval = [obj.point2Center(c(1:obj.nPoints))
                            zeros(2,obj.nElements)
                            obj.point2Center(c(1+obj.nPoints:2*obj.nPoints))];


                    % (iic) c is a 2 x nElements matrix.
                    % It will be handled as a 2 x 2 diagonal matrix where
                    % c(1,1) and c(2,2) contain the value of a space
                    % dependent coefficient function at element's centers
                    % (and c(1,2) = c(2,1) = 0).
                    elseif rows == 2 && cols == obj.nElements
                        cval = [c(1:obj.nElements)
                            zeros(2,obj.nElements)
                            c(2,1:obj.nElements)];

                        % (iid) c is a 2 x nPointss matrix.
                    % It will be handled as a 2 x 2 diagonal matrix where
                    % c(1,1) and c(2,2) contain the value of a space
                    % dependent coefficient function at mesh's nodes
                    % (and c(1,2) = c(2,1) = 0).
                    elseif rows == 2 && cols == obj.nPoints
                        cval = [obj.point2Center(c(1,:))
                            zeros(1,obj.nElements)
                            zeros(1,obj.nElements)
                            obj.point2Center(c(2,:))];


                    % (iiia) c is a 2 x 2 matrix
                    elseif rows == 2 && cols == 2
                        c1 = c(1,1); c2 = c(1,2);
                        c3 = c(2,1); c4 = c(2,2);
                        cval = [c1(ones(1,obj.nElements))
                            c2(ones(1,obj.nElements))
                            c3(ones(1,obj.nElements))
                            c4(ones(1,obj.nElements))];

                    % (iiib) c is a 2 x 2*nElements matrix.
                    % It will be handled as a 2 x 2 full matrix where
                    % c(i,j) contain the values of space dependent coefficient
                    % functions.
                    elseif rows == 2 && cols == 2*obj.nElements
                        cval = [c(1,1:obj.nElements)
                            c(1,1+obj.nElements:2*obj.nElements)
                            c(2,1:obj.nElements)
                            c(2,1+obj.nElements:2*obj.nElements)];

                    % (iiic) c is a 2 x 2*nPoints matrix.
                    % It will be handled as a 2 x 2 full matrix where
                    % c(i,j)  contain the value of a space dependent
                    % coefficient function at mesh points.
                    elseif rows == 2 && cols == 2*obj.nPoints
                        cval = [obj.point2Center(c(1,1:obj.nPoints))
                            obj.point2Center(c(1,1+obj.nPoints:2*obj.nPoints))
                            obj.point2Center(c(2,1:obj.nPoints))
                            obj.point2Center(c(2,1+obj.nPoints:2*obj.nPoints))];
                    else
                        errordlg('C has wrong format.','Error')
                        obj.wrongFormat.throwAsCaller
                    end%
                case {'char', 'string'}
                    % must be a single char/string symbolizing
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
                            obj.wrongFormat.throwAsCaller;
                    end
                case 'cell'
                    % must be a 2*2 cell array which refers to the entries
                    % of a 2*2
                    % diffussion matrix. Every Entry must be double or a single char symbolizing
                    % the part of the coefficent function like in the case
                    % 'char'.

                    [rows,cols] = size(c);
                    if (rows ~=2) || (cols ~=2)
                        obj.wrongFormat.throwAsCaller;
                    end
                    cval = zeros(4,obj.nElements);
                    cvalcell = cell(2,2);


                    % For evaluating the coefficient function,
                    % we need x and y variables "hanging in the air".
                    % The next warning can be ignored!
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
                    obj.wrongClass
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
                             MException('GRIDD:WRONGFORMAT',...
                             'Data for source term is missformed').throwAsCaller;
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
                            obj.wrongFormat.throwAsCaller
                    end
                case 'cell'
                    error('Sorry, Cell array input not jet implemented!')
                otherwise
                     obj.wrongClass.throwAsCaller;
            end
        end

        % utilities: diameter, center of triangle, has triangle a boundary edge etc.

        function diam = triangleDiameters(obj)
            % diam = obj.triangleDiameters()
            % Computes the (outer) diameter of all triangles
            d1 = (obj.p(:,obj.t(1,:))-obj.p(:,obj.t(2,:)));
            d2 = (obj.p(:,obj.t(1,:))-obj.p(:,obj.t(3,:)));
            d3 = (obj.p(:,obj.t(2,:))-obj.p(:,obj.t(3,:)));

            diam(1,:) = sqrt(sum(d1.*d1));
            diam(2,:) = sqrt(sum(d2.*d2));
            diam(3,:) = sqrt(sum(d3.*d3));
            diam = max(diam);
        end

        function b = isBoundaryTriangle(obj)
            %  b = isBoundaryTriangle(gt)
            % Method that indicates the boundary points.
            % b is a vector of length(#triangles)    %

            % Algorithmus: falls in der k-ten Spalte von N
            % Nullen stehen, hat das Dreieck weniger als drei
            % Nachbarn, ist also ein Dreieck mit Randkante
            N = obj.neighbours();
            b = false(1,size(obj.t,2));
            for k = 1:size(N,2)
                n = 0; % keine Randkanten
                for j = 1:3
                    if N(j,k)==0
                        n = n+1;
                    end
                end
                b(k) = (n>0);
            end
        end

        function indx = point2ElementIndex(obj,pt)
            % point2ElementIndex gives back the index of the triangle
            % containing the point pt pt must be a vector of length
            % two.
            % indx = point2ElementIndex(obj,pt)
            if ~isvector(pt)
                MException('GRID2D:NOVECTOR',...
                'The argument must be the coordinates of a single point').throwAsCaller;
            end

            p = obj.p;
            t = obj.t;
            p1 = p(:,(t(1,:)));
            p2 = p(:,(t(2,:)));
            p3 = p(:,(t(3,:)));
            x21 = p2(1,:)-p1(1,:);
            x31 = p3(1,:)-p1(1,:);
            y21 = p2(2,:)-p1(2,:);
            y31 = p3(2,:)-p1(2,:);
            J = x21.*y31-x31.*y21;

            ptx = pt(1)-p1(1,:);
            pty = pt(2)-p1(2,:);
            ptux = 1./J.*(y31.*ptx-x31.*pty);
            ptuy = 1./J.*(x21.*pty-y21.*ptx);

            indx = find((sum([ptux;ptuy])<=1) & (ptux>=0) & (ptuy>=0));
        end

        function [sidelength,area] = sideLengthAndArea(obj)
        %Method to compute  side lengths and areas of triangles.
            dx = zeros(3,obj.nElements);
            dy = zeros(3,obj.nElements);
            sidelength = zeros(3,obj.nElements);
            for k = 1:3
                k1 = rem(k,3)+1;
                k2 = rem(k1,3)+1;
                dx(k,:) = obj.p(1,obj.t(k1,:)) - obj.p(1,obj.t(k2,:));
                dy(k,:) = obj.p(2,obj.t(k1,:)) - obj.p(2,obj.t(k2,:));
                sidelength(k,:) = sqrt(dx(k,:).*dx(k,:) + dy(k,:).*dy(k,:));
            end
            area = 0.5*abs(dx(1,:).*dy(2,:) - dx(2,:).*dy(1,:));
        end

        function jiggle(obj,method)
            %obj.jiggle
            % Jiggles the mesh
            % (c) 2013 by Uwe Prüfert
            if nargin==1
                method = 'all';
            end

            innerPoints = obj.getInnerPointsIndex;
            switch method
                case 'mean'
                    q = obj.meshQuality;
                    [~,idx] = sort(q,'descend');
                    points2Jiggle = intersect(innerPoints,obj.t(1:3,idx))';
                case 'all'
                    points2Jiggle = innerPoints;
                otherwise
            end

            for k = points2Jiggle
                indx = obj.getNeighbourPointsIndex(k);
                % mean value
                obj.p(:,k) = sum(obj.p(:,indx),2)/length(indx);
            end
        end

        function b = isPointInDomain(obj,pt)
            % isPointInTriangle -  gives back a bool vector.
            % b ist true if the point pt is i the triangle of a given
            % mesh object.
            % b = obj.isPointInTriangle(pt)
            % Algortihm: Transforms the point into the unit triangle  and
            % decide based on the relation pt_x >= 0 pt_y >= 0 pt_x+pt_y <=1.
            p = obj.p;
            t = obj.t;
            p1 = p(:,(t(1,:)));
            p2 = p(:,(t(2,:)));
            p3 = p(:,(t(3,:)));
            x21 = p2(1,:)-p1(1,:);
            x31 = p3(1,:)-p1(1,:);
            y21 = p2(2,:)-p1(2,:);
            y31 = p3(2,:)-p1(2,:);
            J = x21.*y31-x31.*y21;

            b = isPointInAnyTriangle(pt,p1,x21,x31,y21,y31,J);


            function b = isPointInAnyTriangle(pt,p1,x21,x31,y21,y31,J)
                b = false(1,size(pt,2));
                nP = size(pt,2);
                % The not vectorized code, compare it with the C-code in
                % isPointInAnyTriangle.c
                % nE = size(p1,2);
                % for k = 1:nP
                %     for l = 1:nE
                %         ptx = pt(1,k)-p1(1,l);
                %         pty = pt(2,k)-p1(2,l);
                %         ptux = 1/J(l)*(y31(l)*ptx-x31(l)*pty);
                %         ptuy = 1/J(l)*(x21(l)*pty-y21(l)*ptx);
                %         if  (ptux+ptuy<=1 && ptux>=0 && ptuy>=0)
                %             b(k) = true;
                %             break
                %         end
                %     end
                % end
                for k = 1:nP
                    ptx = pt(1,k)-p1(1,:);
                    pty = pt(2,k)-p1(2,:);
                    ptux = 1./J.*(y31.*ptx-x31.*pty);
                    ptuy = 1./J.*(x21.*pty-y21.*ptx);
                    b(k) = any((sum([ptux;ptuy])<=1) & (ptux>=0) & (ptuy>=0));
                end
            end
        end

        function N = neighbours(obj,varargin)
            % N = neigbours(gt[,indx])
            % computes to every triangle the
            % index of neighbored triangles
            if nargin==2
                indx = varargin{1};
            else
                indx = 1:obj.nElements;
            end
            N = zeros(3,obj.nElements);
            for k = 1:length(indx)
                nk = obj.ent(indx(k));
                [i] = find(nk==k);
                nk(i)=[];
                N(1:length(nk),indx(k))=nk;
            end
        end

        function intl = ent(obj,it)
            % ent Indices of triangles neighboring
            % a given set of triangles.
            nt=size(obj.t,2);
            switch nargin
                case 1
                    % all neighbours
                    it = 1:obj.nElements;
                case 2
                    % okay
                otherwise
                    obj.wrongNumberInputs.throw;
            end
            it1=ones(1,obj.nElements);
            it1(it)=zeros(size(it));
            it1=find(it1);
            ip1= obj.t(1,it)';
            ip2= obj.t(2,it)';
            ip3= obj.t(3,it)';

            % Make a connectivity matrix.
            A=sparse(ip1,ip2,1,obj.nPoints,obj.nPoints);
            A=A+sparse(ip2,ip3,1,obj.nPoints,obj.nPoints);
            A=A+sparse(ip3,ip1,1,obj.nPoints,obj.nPoints);


            ntl=zeros(1,nnz(A-A')); % a slight overestimate
            nnt=0;

            for i=1:length(it1)
                if A( obj.t(2,it1(i)), obj.t(1,it1(i))) || ...
                        A( obj.t(3,it1(i)), obj.t(2,it1(i))) || ...
                        A( obj.t(1,it1(i)), obj.t(3,it1(i)))
                    nnt=nnt+1;
                    ntl(nnt)=it1(i);
                end
            end
            intl=sort([it ntl(1:nnt)]);
        end
    end



end
