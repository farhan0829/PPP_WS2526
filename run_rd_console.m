function [u, v, geo, t] = run_rd_console(modelName, a, b, Du, Dv)
%RUN_RD_CONSOLE  Minimal console runner for Schnakenberg / Brusselator / Heat test.
%
% Usage (from MATLAB Command Window):
%   run_rd_console('schnakenberg', 0.1, 0.9, 1, 8);
%   run_rd_console('brusselator',  1.0, 3.0, 1, 8);
%   run_rd_console('heat',         0,   0,   1, 1);   % "all parameters 0" test
%
% Inputs:
%   modelName : 'schnakenberg' | 'brusselator' | 'heat'
%   a,b       : model parameters (ignored for 'heat')
%   Du,Dv     : diffusion coefficients (Dv ignored for 'heat' unless you want v to diffuse too)
%
% Requirements (on MATLAB path):
%   - Bilinear2D.m
%   - RectangleR.m (from your OOPDELight / FEM toolbox)
%
% This is intentionally "no-frills": fixed domain + timestep + final time.
% Edit the "Defaults" section below if you want different dt/tEnd/mesh.

    % -------------------------
    % Defaults (edit if needed)
    % -------------------------
    dt        = 0.01;     % time step
    tEnd      = 50;       % final time
    stepsPerFrame = 10;   % plot every N steps
    Lx = 100; Ly = 50;    % domain size (rectangle)
    h  = 2.0;             % mesh width (smaller = finer mesh)
    cmap = 'parula';      % plot colormap

    % ---------------
    % Input checking
    % ---------------
    if nargin < 1
        % Interactive console mode
        fprintf("Model options: schnakenberg | brusselator | heat\n");
        fprintf("Suggested ranges:\n");
        fprintf("  Schnakenberg: a in [0, 0.2], b in [0, 2]\n");
        fprintf("  Brusselator : A in [0.5, 2], B in [1, 5]\n");
        fprintf("  Diffusion   : Du,Dv in [0.1, 100]\n\n");

        modelName = input("Model name [schnakenberg]: ", "s");
        if isempty(modelName)
            modelName = "schnakenberg";
        end

        modelName = lower(string(modelName));
        if modelName == "heat"
            a = 0;
            b = 0;
        else
            a = promptNumeric("a (or A) [0.1]: ", 0.1);
            b = promptNumeric("b (or B) [0.9]: ", 0.9);
        end
        Du = promptNumeric("Du [1]: ", 1);
        Dv = promptNumeric("Dv [8]: ", 8);
    else
        if nargin < 5
            error("Usage: run_rd_console(modelName, a, b, Du, Dv)");
        end
        modelName = lower(string(modelName));
    end

    % -------------------------
    % Grid + FEM + matrices
    % -------------------------
    fem = Bilinear2D;
    geo = RectangleR(0, Lx, 0, Ly, h);

    % Neumann-ish default (Robin(0,1) is your toolbox's safe neutral option)
    geo.makeBoundaryMatrix(geo.robinBC(0,1));

    [K, M, F] = fem.assema(geo, 1, 1, 1);
    [Q, G, ~, ~] = fem.assemb(geo); %#ok<ASGLU>

    n = geo.nPoints;

    % -------------------------
    % Initial condition
    % -------------------------
    % Uniform + small noise
    rng(0);
    u = ones(n,1) + 0.01*(rand(n,1)-0.5);
    v = ones(n,1) + 0.01*(rand(n,1)-0.5);

    % For heat-test, start from a localized bump (so you clearly see diffusion)
    if modelName == "heat"
        cx = 0.5*(min(geo.x)+max(geo.x));
        cy = 0.5*(min(geo.y)+max(geo.y));
        bump = exp(-((geo.x-cx).^2 + (geo.y-cy).^2)/(2*(0.08*min(Lx,Ly))^2));
        u = bump(:);
        v = zeros(n,1);
    end

    % enforce column vectors (prevents vertcat size errors)
    u = u(:); v = v(:);

    % -------------------------
    % Model coefficients
    % -------------------------
    % We write the PDE as:
    %   u_t = Du * Lap(u) + f(u,v)
    %   v_t = Dv * Lap(v) + g(u,v)
    %
    % Implicit Euler on diffusion + linear reaction, explicit on nonlinear term (u^2 v).
    %
    % For Schnakenberg:
    %   f = a - u + u^2 v
    %   g = b - u^2 v
    %
    % For Brusselator:
    %   f = a - (b+1)u + u^2 v     where a = A, b = B
    %   g = b*u - u^2 v
    %
    % For Heat test:
    %   f = 0, g = 0   (pure diffusion)

    s = 0; % boundary source term off (keep 0)

    switch modelName
        case "schnakenberg"
            Au = a;               Av = b;
            Buu = -1;             Buv = 0;
            Bvu = 0;              Bvv = 0;
            Cu  =  1;             Cv  = -1;    % nonlinear term (+u^2 v, -u^2 v)

        case "brusselator"
            Au = a;               Av = 0;
            Buu = -(b+1);         Buv = 0;
            Bvu =  b;             Bvv = 0;
            Cu  =  1;             Cv  = -1;

        case "heat"
            Au = 0;               Av = 0;
            Buu = 0;              Buv = 0;
            Bvu = 0;              Bvv = 0;
            Cu  = 0;              Cv  = 0;

        otherwise
            error("Unknown modelName '%s'. Use 'schnakenberg', 'brusselator', or 'heat'.", modelName);
    end

    % -------------------------
    % Build constant system
    % -------------------------
    A = [ M - dt*(Buu*M - Du*K + s*Q),   -dt*Buv*M; ...
          -dt*Bvu*M,                     M - dt*(Bvv*M - Dv*K + s*Q) ];

    Cvec = dt * [ (Au*F + s*G); ...
                  (Av*F + s*G) ];

    % -------------------------
    % Plot initial frame
    % -------------------------
    fig = figure('Name', sprintf('RD console: %s', modelName));
    ax = axes('Parent', fig);
    setappdata(fig, 'stop', false);
    uicontrol('Parent', fig, ...
        'Style', 'pushbutton', ...
        'String', 'Stop', ...
        'Units', 'pixels', ...
        'Position', [10 10 60 25], ...
        'Callback', @(src, evt) setappdata(fig, 'stop', true)); %#ok<NASGU>
    colormap(ax, cmap);
    plotFrame(ax, geo, u, 0);

    % -------------------------
    % Time stepping
    % -------------------------
    t = 0;
    k = 0;
    while t < tEnd
        if ~isvalid(fig) || getappdata(fig, 'stop')
            break;
        end
        % Nonlinear term u^2 v (explicit)
        nl = u.^2 .* v;

        ll = [M*u; M*v];                 % column (2n x 1)
        NL = dt * [ Cu * (M*nl); ...
                    Cv * (M*nl) ];

        y = A \ (ll + Cvec + NL);

        u = y(1:n);
        v = y(n+1:end);

        t = t + dt;
        k = k + 1;

        if mod(k, stepsPerFrame) == 0
            plotFrame(ax, geo, u, t);
            drawnow limitrate
        end
    end

    plotFrame(ax, geo, u, t);
    drawnow;

end

% ===== local helper (OK inside function file) =====
function plotFrame(ax, geo, z, t)
    if ~isgraphics(ax, 'axes')
        ax = axes('Parent', gcf);
    end
    cla(ax);
    geo.plot(ax, z, 'EdgeColor', 'none');  % same style as simulatorGUI
    view(ax, 2);
    axis(ax, 'tight');
    xlabel(ax, 'x');
    ylabel(ax, 'y');
    title(ax, sprintf('t = %.3g', t));
    colorbar(ax);
end

function val = promptNumeric(msg, defaultVal)
    raw = input(msg, "s");
    if isempty(raw)
        val = defaultVal;
        return;
    end
    val = str2double(raw);
    if isnan(val)
        fprintf("Invalid input. Using default: %g\n", defaultVal);
        val = defaultVal;
    end
end
