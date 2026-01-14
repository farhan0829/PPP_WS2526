classdef Rectangle < grid2D
    %UNTITLED9 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = Rectangle(varargin)
            % square()
            % square(h)
            % square(xmin,xmax,ymin,ymax)
            % square(xmin,xmax,ymin,ymax,h)
            % No options => unitsqaure with h = 0.1
            % hmax given but but no corner points => unitsquare with given
            % left lower and right upper point given,
            % but no hmax => square with hmax = 0.1
            % The boundary lenght parameter runs in anti clockwise
            % order from 0 to 1 on everey boundary segment,
            % independetly from its length!
            % (c) 2013 by Uwe Prüfert

            % Code valid for MATLAB R>= 2013a by using new
            % delaunayTriangulation class. For older Matlab
            % Releases it uses the old
            % DelaunayTri Class

            % First handle all posible argument...
            switch nargin
                case 0
                    % default constructor
                   return
                case 1
                    % hmax but no coordinates, unitsquare with given
                    % hmax
                    xmin = 0;
                    xmax = 1;
                    ymin = 0;
                    ymax = 1;
                    h = varargin{1};
                case 4
                    % left lower and right upper point given,
                    % but no hmax, square with hmax = 0.1
                    xmin = varargin{1};
                    xmax = varargin{2};
                    ymin = varargin{3};
                    ymax = varargin{4};
                    h = 0.1;
                case 5
                    xmin = varargin{1};
                    xmax = varargin{2};
                    ymin = varargin{3};
                    ymax = varargin{4};
                    h = varargin{5};
                otherwise
                     obj.wrongNumberInputs.throw;
            end

            nx = max(2,round((xmax-xmin)/h)+1);
            ny = max(2,round((ymax-ymin)/h)+1);


            xmesh = linspace(xmin,xmax,nx);
            ymesh = linspace(ymin,ymax,ny);
            [X,Y] = meshgrid(xmesh,ymesh);
            P = [reshape(X,1,nx*ny);reshape(Y,1,nx*ny)];



             
            dt = delaunayTriangulation(P');
            obj.p = dt.Points';
            obj.e = dt.freeBoundary';
            obj.t = dt.ConnectivityList';


            % because we have an equidistant boundary, we can
            % compute the boundary "s" in this simple way...
            enx = linspace(0,1,nx);
            eny = linspace(0,1,ny);
            obj.e(3,1:nx-1) = enx(1:end-1);
            obj.e(3,nx:nx+ny-2) =  eny(1:end-1);
            obj.e(3,nx+ny-1:nx+ny+nx-3) = (enx(1:end-1));
            obj.e(3,nx+ny+nx-2:2*nx+2*ny-4) = (eny(1:end-1)) ;

            obj.e(4,1:nx-1) = enx(2:end);
            obj.e(4,nx:nx+ny-2) =  eny(2:end);
            obj.e(4,nx+ny-1:nx+ny+nx-3) =  (enx(2:end));
            obj.e(4,nx+ny+nx-2:2*nx+2*ny-4) =  (eny(2:end)) ;

            % set the boundary segment number...
            obj.e(5,1:nx-1) = ones(1,nx-1);
            obj.e(5,nx:nx+ny-2) = 2*ones(1,ny-1);
            obj.e(5,nx+ny-1:nx+ny+nx-3) = 3*ones(1,nx-1);
            obj.e(5,nx+ny+nx-2:2*nx+2*ny-4) = 4*ones(1,ny-1);


            obj.t(4,:) = ones(1,size(obj.t,2));
        end
    end    
end

