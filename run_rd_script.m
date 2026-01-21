%RUN_RD_SCRIPT  Standalone script for Schnakenberg / Brusselator / Heat.
%
% Edit the parameters below, then run this script.
% If a == 0 and b == 0, it switches to the heat equation test.

% -------------------------
% User parameters
% -------------------------
modelName = "schnakenberg";   % "schnakenberg" | "brusselator" | "heat"
a = 0.1;                      % a (or A for Brusselator)
b = 0.9;                      % b (or B for Brusselator)
Du = 1;                       % diffusion for u
Dv = 8;                       % diffusion for v

% -------------------------
% Defaults (edit if needed)
% -------------------------
dt        = 0.01;     % time step
tEnd      = 50;       % final time
stepsPerFrame = 10;   % plot every N steps
Lx = 100; Ly = 50;    % domain size (rectangle)
h  = 2.0;             % mesh width (smaller = finer mesh)
cmap = 'parula';      % plot colormap

% -------------------------
% Auto-select heat test
% -------------------------
if a == 0 && b == 0
    modelName = "heat";
end

modelName = lower(string(modelName));

% -------------------------
% Grid + FEM + matrices
% -------------------------
fem = Bilinear2D;
geo = RectangleR(0, Lx, 0, Ly, h);

% Neumann-ish default (Robin(0,1) is a neutral option)
geo.makeBoundaryMatrix(geo.robinBC(0,1));

[K, M, F] = fem.assema(geo, 1, 1, 1);
[Q, G, ~, ~] = fem.assemb(geo); %#ok<ASGLU>

n = geo.nPoints;

% -------------------------
% Initial condition
% -------------------------
rng(0);
u = ones(n,1) + 0.01*(rand(n,1)-0.5);
v = ones(n,1) + 0.01*(rand(n,1)-0.5);

% For heat-test, start from a localized bump
if modelName == "heat"
    cx = 0.5*(min(geo.x)+max(geo.x));
    cy = 0.5*(min(geo.y)+max(geo.y));
    bump = exp(-((geo.x-cx).^2 + (geo.y-cy).^2)/(2*(0.08*min(Lx,Ly))^2));
    u = bump(:);
    v = zeros(n,1);
end

u = u(:); v = v(:);

% -------------------------
% Model coefficients
% -------------------------
s = 0; % boundary source term off (keep 0)

switch modelName
    case "schnakenberg"
        Au = a;               Av = b;
        Buu = -1;             Buv = 0;
        Bvu = 0;              Bvv = 0;
        Cu  =  1;             Cv  = -1;

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
fig = figure('Name', sprintf('RD script: %s', modelName));
ax = axes(fig);
colormap(ax, cmap);
plotFrame(ax, geo, u, 0);

% -------------------------
% Time stepping
% -------------------------
t = 0;
k = 0;
while t < tEnd
    nl = u.^2 .* v;

    ll = [M*u; M*v];
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

% ===== local helper =====
function plotFrame(ax, geo, z, t)
    cla(ax);
    geo.plot(ax, z, 'EdgeColor', 'none');
    view(ax, 2);
    axis(ax, 'tight');
    title(ax, sprintf('t = %.3g', t));
    colorbar(ax);
end
