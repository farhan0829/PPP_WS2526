classdef RectangleR < grid2DR
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
     
    
    methods
        function obj = RectangleR(varargin)
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
                     %
            end

            nx = max(2,round((xmax-xmin)/h)+1);
            ny = max(2,round((ymax-ymin)/h)+1);


            xmesh = linspace(xmin,xmax,nx);
            ymesh = linspace(ymin,ymax,ny);
            obj.rectangle(xmesh,ymesh,h);
        end       
         
    end
end

