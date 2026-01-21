classdef Simulator_GUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        HelpMenu                        matlab.ui.container.Menu
        ManualMenu                      matlab.ui.container.Menu
        InitialConditionsandControlsPanel  matlab.ui.container.Panel
        NoiseAmplitudeSliderLabel       matlab.ui.control.Label
        StopButton                      matlab.ui.control.Button
        ResetButton                     matlab.ui.control.Button
        StepButton                      matlab.ui.control.Button
        StartPauseButton                matlab.ui.control.StateButton
        NoiseAmplitudeSlider            matlab.ui.control.Slider
        BoundaryConditionsButtonGroup   matlab.ui.container.ButtonGroup
        PeriodicButton                  matlab.ui.control.RadioButton
        DirichletButton                 matlab.ui.control.RadioButton
        NeumannButton                   matlab.ui.control.RadioButton
        vBaseValueEditField             matlab.ui.control.NumericEditField
        vBaseValueEditFieldLabel        matlab.ui.control.Label
        uBaseValueEditField             matlab.ui.control.NumericEditField
        uBaseValueEditFieldLabel        matlab.ui.control.Label
        InitialConditionDropDown        matlab.ui.control.DropDown
        InitialConditionDropDownLabel   matlab.ui.control.Label
        AdvancedFeaturesPanel           matlab.ui.container.Panel
        PatternsButtonGroup             matlab.ui.container.ButtonGroup
        ComplexPatternsButton           matlab.ui.control.RadioButton
        StripesButton                   matlab.ui.control.RadioButton
        TuringSpotsButton               matlab.ui.control.RadioButton
        EnablerealtimeupdatesCheckBox   matlab.ui.control.CheckBox
        AutoscalecolorlimitsCheckBox    matlab.ui.control.CheckBox
        ColormapDropDown                matlab.ui.control.DropDown
        ColormapDropDownLabel           matlab.ui.control.Label
        PlottypeDropDown                matlab.ui.control.DropDown
        PlottypeDropDownLabel           matlab.ui.control.Label
        VariabletodisplayDropDown       matlab.ui.control.DropDown
        VariabletodisplayDropDownLabel  matlab.ui.control.Label
        ExportandStatusPanel            matlab.ui.container.Panel
        HistoryRecordingCheckBox        matlab.ui.control.CheckBox
        FrameRateSpinner                matlab.ui.control.Spinner
        FrameRateSpinnerLabel           matlab.ui.control.Label
        ExportAnimationButton           matlab.ui.control.Button
        SaveDataButton                  matlab.ui.control.Button
        CurrenttimeLabel                matlab.ui.control.Label
        StatusLabel                     matlab.ui.control.Label
        SimulationProgressGauge         matlab.ui.control.LinearGauge
        SimulationProgressGaugeLabel    matlab.ui.control.Label
        VisualizationPanel              matlab.ui.container.Panel
        UIAxes                          matlab.ui.control.UIAxes
        ModelandParametersPanel         matlab.ui.container.Panel
        bSliderLabel                    matlab.ui.control.Label
        aSliderLabel                    matlab.ui.control.Label
        aSlider                         matlab.ui.control.Slider
        UpdateEveryNstepsSpinner        matlab.ui.control.Spinner
        UpdateEveryNstepsSpinnerLabel   matlab.ui.control.Label
        TotalTimeEditField              matlab.ui.control.NumericEditField
        TotalTimeEditFieldLabel         matlab.ui.control.Label
        TimeStepdtEditField             matlab.ui.control.NumericEditField
        TimeStepdtEditFieldLabel        matlab.ui.control.Label
        GridPointsSpinner               matlab.ui.control.Spinner
        GridPointsSpinnerLabel          matlab.ui.control.Label
        DomainSizeEditField             matlab.ui.control.NumericEditField
        DomainSizeEditFieldLabel        matlab.ui.control.Label
        DvSpinner                       matlab.ui.control.Spinner
        DvSpinnerLabel                  matlab.ui.control.Label
        DuSpinner                       matlab.ui.control.Spinner
        bSlider                         matlab.ui.control.Slider
        DuSpinnerLabel                  matlab.ui.control.Label
        ModelTypeDropDown               matlab.ui.control.DropDown
        ModelTypeDropDownLabel          matlab.ui.control.Label
    end

    
    properties (Access = private)
        % ==== Backend (OOPDELight FEM) ====
        model
        fem
        geo
    
        % Cached matrices/vectors
        K; M; Q; G; F;
        H; R;          % for Dirichlet constraints
        A; Cvec;
        n
    
        % State
        u; v
        t double = 0
        k double = 0
    
        % Timer
        timerObj timer
    
    % History for export/save
          frameHistory = struct('cdata',{},'colormap',{});
          tHistory = [];
          uHistory = {};
          vHistory = {};
    end
    methods (Access = private)
        
        % =========================
        % Backend helpers (MUST exist)
        % =========================

        function clearHistory(app)
              app.frameHistory = struct('cdata',{},'colormap',{});
              app.tHistory = [];
              app.uHistory = {};
              app.vHistory = {};
          end

        function captureHistory(app)
              if ~app.HistoryRecordingCheckBox.Value || isempty(app.u) || isempty(app.v) || isempty(app.UIAxes)
                  return;
              end
              app.frameHistory(end+1) = getframe(app.UIAxes);
              app.tHistory(end+1) = app.t;
              app.uHistory{end+1} = app.u;
              app.vHistory{end+1} = app.v;
          end

        function params = getAllParameters(app)
              if app.PeriodicButton.Value
                  bc = 'Periodic';
              elseif app.DirichletButton.Value
                  bc = 'Dirichlet';
              elseif app.NeumannButton.Value
                  bc = 'Neumann';
              else
                  bc = '';
              end

              adv = '';
              if ~isempty(app.PatternsButtonGroup.SelectedObject)
                  adv = app.PatternsButtonGroup.SelectedObject.Text;
              end

              params = struct( ...
                  'Model', app.ModelTypeDropDown.Value, ...
                  'a', app.aSlider.Value(1), ...
                  'b', app.bSlider.Value(1), ...
                  'Du', app.DuSpinner.Value, ...
                  'Dv', app.DvSpinner.Value, ...
                  'DomainSize', app.DomainSizeEditField.Value, ...
                  'GridPoints', app.GridPointsSpinner.Value, ...
                  'dt', app.TimeStepdtEditField.Value, ...
                  'TotalTime', app.TotalTimeEditField.Value, ...
                  'UpdateEveryNSteps', app.UpdateEveryNstepsSpinner.Value, ...
                  'InitialCondition', app.InitialConditionDropDown.Value, ...
                  'NoiseAmplitude', app.NoiseAmplitudeSlider.Value, ...
                  'uBase', app.uBaseValueEditField.Value, ...
                  'vBase', app.vBaseValueEditField.Value, ...
                  'BoundaryCondition', bc, ...
                  'AdvancedPattern', adv, ...
                  'RealtimeUpdates', app.EnablerealtimeupdatesCheckBox.Value, ...
                  'AutoscaleColorLimits', app.AutoscalecolorlimitsCheckBox.Value, ...
                  'Colormap', app.ColormapDropDown.Value, ...
                  'PlotType', app.PlottypeDropDown.Value, ...
                  'VariableToDisplay', app.VariabletodisplayDropDown.Value, ...
                  'FrameRate', app.FrameRateSpinner.Value, ...
                  'HistoryRecording', app.HistoryRecordingCheckBox.Value ...
              );
          end

        function resetDashboardUI(app)
            % Stop any running sim
            if ~isempty(app.timerObj) && isvalid(app.timerObj)
                stop(app.timerObj);
            end
            app.StartPauseButton.Value = false;
        
            % ---- UI defaults (choose what you want as "Reset") ----
            app.ModelTypeDropDown.Items = {'Schnakenberg','Brusselator','Heat'};
            app.ModelTypeDropDown.Value = 'Schnakenberg';
        
            app.aSlider.Limits = [0.0 4.0];
            app.bSlider.Limits = [0.0 4.0];

            app.aSlider.Value  = 0.1;   % scalar
            app.bSlider.Value  = 0.9;   % scalar
            app.ComplexPatternsButton.Value = true;
            app.ModelTypeDropDown.Enable = 'on';
            app.aSlider.Enable = 'on';
            app.bSlider.Enable = 'on';
            app.DuSpinner.Enable = 'on';
            app.DvSpinner.Enable = 'on';
            app.DuSpinnerLabel.Enable = 'on';
            app.DvSpinnerLabel.Enable = 'on';
            app.aSliderLabel.Text = 'a';
            app.bSliderLabel.Text = 'b';
        
            app.DomainSizeEditField.Value = 50;
            app.GridPointsSpinner.Value  = 500;
            app.TimeStepdtEditField.Value = 0.1;
            app.TotalTimeEditField.Value  = 200;
            app.UpdateEveryNstepsSpinner.Value = 10;
        
            app.uBaseValueEditField.Value = 1;
            app.vBaseValueEditField.Value = 1;
        
            app.NoiseAmplitudeSlider.Value = 0.01;
            app.InitialConditionDropDown.Items = {'Uniform','Random','Localized Spot'};
            app.InitialConditionDropDown.Value = 'Random';
        
            app.VariabletodisplayDropDown.Items = {'u','v','Both (u+v)','Magnitude'};
            app.VariabletodisplayDropDown.Value = 'u';
        
            app.PlottypeDropDown.Items = {'2D Heatmap','3D Surface','Contour'};
            app.PlottypeDropDown.Value = '2D Heatmap';
        
            app.ColormapDropDown.Items = {'parula','jet','hot','cool','turbo','gray'};
            app.ColormapDropDown.Value = 'turbo';
        
            app.EnablerealtimeupdatesCheckBox.Value = true;
            app.AutoscalecolorlimitsCheckBox.Value = true;
        
            % Boundary default
            app.NeumannButton.Value   = false;
            app.DirichletButton.Value = false;
            app.PeriodicButton.Value  = true;
        
            % Fix Du enable 
            app.DuSpinner.Enable = 'on';
            app.DuSpinnerLabel.Enable = 'on';
        
            % Status + progress + axes clear
            app.t = 0; app.k = 0;
            app.StatusLabel.Text = 'Ready.';
            app.CurrenttimeLabel.Text = 't = 0';
            app.SimulationProgressGauge.Value = 0;
        
            cla(app.UIAxes);
            title(app.UIAxes,'');
            drawnow;
        end



        
        function resetSimulation(app)
      % ---- Fix broken dropdown strings coming from createComponents (grey) ----
      app.ModelTypeDropDown.Items = {'Schnakenberg','Brusselator','Heat'};
      if ~ismember(app.ModelTypeDropDown.Value, app.ModelTypeDropDown.Items)
          app.ModelTypeDropDown.Value = 'Schnakenberg';
      end

      % ---- Build model ----
      modelName = app.ModelTypeDropDown.Value;
      if strcmpi(modelName,'Heat')
          app.model = setModel();
      else
          app.model = setModel(modelName);   % uses your existing file setModel.m
      end

      % Read parameters (RangeSlider -> take first value)
      a = app.aSlider.Value(1);
      b = app.bSlider.Value(1);

      % Base coefficients per model
      if strcmpi(modelName,'Schnakenberg')
          setModelField('Au', a);
          setModelField('Av', b);
          setModelField('Buu', -1);
          setModelField('Buv', 0);
          setModelField('Bvu', 0);
          setModelField('Bvv', 0);
          setModelField('Cu', 1);
          setModelField('Cv', -1);
      elseif strcmpi(modelName,'Brusselator')
          setModelField('Au', a);
          setModelField('Av', 0);
          setModelField('Buu', -(b+1));
          setModelField('Buv', 0);
          setModelField('Bvu', b);
          setModelField('Bvv', 0);
          setModelField('Cu', 1);
          setModelField('Cv', -1);
      elseif strcmpi(modelName,'Heat')
          a = 0; b = 0;
          setModelField('Au', 0);
          setModelField('Av', 0);
          setModelField('Buu', 0);
          setModelField('Buv', 0);
          setModelField('Bvu', 0);
          setModelField('Bvv', 0);
          setModelField('Cu', 0);
          setModelField('Cv', 0);
       end

      % Push parameters into model (handle both struct/object safely)
      setModelField('dt', app.TimeStepdtEditField.Value);
      setModelField('Du', app.DuSpinner.Value);
      setModelField('Dv', app.DvSpinner.Value);

      % ---- Geometry + FEM ----
      if exist('Bilinear2D','class') ~= 8
          app.StatusLabel.Text = 'Missing OOPDELight classes on path.';
          return;
      end

      L  = app.DomainSizeEditField.Value;
      Lx = 2*L; Ly = 1*L;

      N = max(50, round(app.GridPointsSpinner.Value));
      h = max(1e-6, L / max(1,(sqrt(N)-1)));
      setModelField('h', h);

      app.fem = Bilinear2D;
      app.geo = RectangleR(0, Lx, 0, Ly, h);

      % n points
      if isprop(app.geo,'nPoints')
          app.n = app.geo.nPoints;
      else
          app.n = length(app.geo.x);
      end

      % Boundary matrix
      if app.DirichletButton.Value
          app.geo.makeBoundaryMatrix(app.geo.dirichletBC('0'));
      elseif app.PeriodicButton.Value
          app.geo.makeBoundaryMatrix(app.geo.robinBC(0,1));
      else
          app.geo.makeBoundaryMatrix(app.geo.robinBC(0,1)); % Neumann-ish safe default
      end

      % Assemble
      [app.K, app.M, app.F] = app.fem.assema(app.geo, 1, 1, 1);
      [app.Q, app.G, app.H, app.R] = app.fem.assemb(app.geo);

      % ---- Initial conditions ----
      u0 = app.uBaseValueEditField.Value * ones(app.n,1);
      v0 = app.vBaseValueEditField.Value * ones(app.n,1);

       x = app.geo.x; y = app.geo.y;
      % Heat-specific initial condition (visible diffusion)
  if strcmpi(modelName,'Heat')
      cx = 0.5*(min(x)+max(x));
      cy = 0.5*(min(y)+max(y));
      bump = exp(-((x-cx).^2 + (y-cy).^2)/(2*(0.08*L)^2));
      u0 = bump(:);
      v0 = zeros(app.n,1);
  end


      ic  = lower(string(app.InitialConditionDropDown.Value));
      amp = app.NoiseAmplitudeSlider.Value;

     
      if contains(ic,"random")
          u0 = u0 + amp*(rand(app.n,1)-0.5);
          v0 = v0 + amp*(rand(app.n,1)-0.5);
      elseif contains(ic,"spot")
          cx = 0.5*(min(x)+max(x));
          cy = 0.5*(min(y)+max(y));
          p = zeros(app.n,1);
          p(((x-cx).^2 + (y-cy).^2) <= (0.08*L)^2) = 1;
          u0 = u0 + p;
      end

      app.u = u0;
      app.v = v0;
      app.t = 0;
      app.k = 0;

      % ---- Build system matrix A and vector Cvec ----
      dt = getModelField('dt', app.TimeStepdtEditField.Value);
      Du = getModelField('Du', app.DuSpinner.Value);
      Dv = getModelField('Dv', app.DvSpinner.Value);

      Buu = getModelField('Buu', 0);
      Buv = getModelField('Buv', 0);
      Bvu = getModelField('Bvu', 0);
      Bvv = getModelField('Bvv', 0);
      Cu  = getModelField('Cu',  1);
      Cv  = getModelField('Cv', -1);

      Au  = getModelField('Au', a);
      Av  = getModelField('Av', b);

      s = 0;
        app.A = [ app.M - dt*(Buu*app.M - Du*app.K + s*app.Q),   -dt*Buv*app.M; -dt*Bvu*app.M,app.M - dt*(Bvv*app.M - Dv*app.K + s*app.Q) ];

      app.Cvec = dt * [ (Au*app.F + s*app.G);
                        (Av*app.F + s*app.G) ];

      % Dirichlet penalty
      if app.DirichletButton.Value && ~isempty(app.H)
          pen = 1e2;
          P = (app.H' * app.H);
          Z = sparse(size(P,1), size(P,2));
          app.A = app.A + pen * [P Z; Z P];
      end

      % Periodic penalty (ties opposite boundaries)
      if app.PeriodicButton.Value
          pen = 1e3;
          P = buildPeriodicPenalty(app.geo);
          Z = sparse(size(P,1), size(P,2));
          app.A = app.A + pen * [P Z; Z P];
      end

      app.StatusLabel.Text = 'Ready.';
      app.CurrenttimeLabel.Text = 't = 0';
      app.SimulationProgressGauge.Value = 0;
      clearHistory(app);
      renderFrame(app);

      % ---- small helpers (inside method is OK) ----
      function setModelField(fname, val)
          if isstruct(app.model)
              app.model.(fname) = val;
          else
              try app.model.(fname) = val; catch, end
          end
      end

      function v = getModelField(fname, defaultV)
          v = defaultV;
          if isstruct(app.model) && isfield(app.model,fname)
              v = app.model.(fname); return;
          end
          if ~isstruct(app.model)
              try v = app.model.(fname); catch, v = defaultV; end
          end
      end

      function P = buildPeriodicPenalty(geo)
          x = geo.x(:);
          y = geo.y(:);
          dx = max(x) - min(x);
          dy = max(y) - min(y);
          tol = max([dx, dy, 1]) * 1e-10;
    
          left = find(abs(x - min(x)) <= tol);
          right = find(abs(x - max(x)) <= tol);
          bottom = find(abs(y - min(y)) <= tol);
          top = find(abs(y - max(y)) <= tol);
    
          pairs = [pairByValue(y, left, right); pairByValue(x, bottom, top)];
          if isempty(pairs)
              P = sparse(geo.nPoints, geo.nPoints);
              return;
          end

          m = size(pairs,1);
          rows = (1:m)';
          n = geo.nPoints;
          C = sparse(rows, pairs(:,1), 1, m, n) + sparse(rows, pairs(:,2), -1, m, n);
          P = C' * C;

      function pairs = pairByValue(coord, idxA, idxB)
          if isempty(idxA) || isempty(idxB)
              pairs = zeros(0,2);
              return;
          end
          [~, ia] = sort(coord(idxA));
          [~, ib] = sort(coord(idxB));
          a = idxA(ia);
          b = idxB(ib);
          mloc = min(numel(a), numel(b));
          pairs = [a(1:mloc), b(1:mloc)];
      end
     end
        end
        function renderFrame(app)
      if isempty(app.geo) || isempty(app.u) || isempty(app.UIAxes) || ~isvalid(app.UIAxes)
          return;
      end

      var = lower(string(app.VariabletodisplayDropDown.Value));
      if var == "v"
          z = app.v;
      elseif contains(var,"both")
          z = app.u + app.v;
      elseif contains(var,"magnitude")
          z = sqrt(app.u.^2 + app.v.^2);
      else
          z = app.u;
      end

      cla(app.UIAxes);

      ptype = app.PlottypeDropDown.Value;
      switch ptype
          case '2D Heatmap'
              if ismethod(app.geo,'plot')
                  app.geo.plot(app.UIAxes, z, 'EdgeColor','none');
              else
                  scatter(app.UIAxes, app.geo.x, app.geo.y, 12, z, 'filled');
              end
              view(app.UIAxes, 2);

          case '3D Surface'
              if isprop(app.geo,'t') && ~isempty(app.geo.t)
                  trisurf(app.geo.t, app.geo.x, app.geo.y, z, ...
                      'Parent', app.UIAxes, 'EdgeColor','none');
                  view(app.UIAxes, 3);
              else
                  scatter(app.UIAxes, app.geo.x, app.geo.y, 12, z, 'filled');
                  view(app.UIAxes, 2);
              end

          case 'Contour'
              if isprop(app.geo,'t') && ~isempty(app.geo.t)
                  tri = triangulation(app.geo.t, app.geo.x, app.geo.y);
                  F = scatteredInterpolant(tri.Points(:,1), tri.Points(:,2), z, 'natural', 'none');
                  xmin = min(app.geo.x); xmax = max(app.geo.x);
                  ymin = min(app.geo.y); ymax = max(app.geo.y);
                  [X,Y] = meshgrid(linspace(xmin,xmax,150), linspace(ymin,ymax,150));
                  Z = F(X,Y);
                  contourf(app.UIAxes, X, Y, Z, 20, 'LineStyle','none');
                  view(app.UIAxes, 2);
              else
                  scatter(app.UIAxes, app.geo.x, app.geo.y, 12, z, 'filled');
                  view(app.UIAxes, 2);
              end

          otherwise
              if ismethod(app.geo,'plot')
                  app.geo.plot(app.UIAxes, z, 'EdgeColor','none');
              else
                  scatter(app.UIAxes, app.geo.x, app.geo.y, 12, z, 'filled');
              end
              view(app.UIAxes, 2);
      end

      axis(app.UIAxes,'tight');
      try
          colormap(app.UIAxes, app.ColormapDropDown.Value);
      catch
          colormap(app.UIAxes, 'parula');
      end
      colorbar(app.UIAxes);

      title(app.UIAxes, sprintf('t = %.4g', app.t));
      drawnow limitrate;
      captureHistory(app);
        end

        function onTick(app)
      if ~isvalid(app) || isempty(app.UIFigure) || ~isvalid(app.UIFigure)
          if ~isempty(app.timerObj) && isvalid(app.timerObj)
              stop(app.timerObj);
          end
          return;
      end

      tEnd = app.TotalTimeEditField.Value;
      if app.t >= tEnd
          StopButtonPushed(app, []);
          app.StatusLabel.Text = 'Finished.';
          return;
      end

      doOneBatch(app);

      if app.EnablerealtimeupdatesCheckBox.Value || app.HistoryRecordingCheckBox.Value
          renderFrame(app);
       
      end
        end
        function doOneBatch(app)
      nSteps = max(1, round(app.UpdateEveryNstepsSpinner.Value));
      if app.EnablerealtimeupdatesCheckBox.Value
          nSteps = 1; % keep UI responsive in realtime mode
      end
      dt = getDT();

      for ii = 1:nSteps
          nl = app.u.^2 .* app.v;

          rhs = [app.M*app.u;
                 app.M*app.v];

          Cu = getModelField('Cu',  1);
          Cv = getModelField('Cv', -1);

          NL = dt * [ Cu * (app.M*nl);
                      Cv * (app.M*nl) ];

          y = app.A \ (rhs + app.Cvec + NL);

          app.u = y(1:app.n);
          app.v = y(app.n+1:end);

          app.t = app.t + dt;
          app.k = app.k + 1;

          if app.t >= app.TotalTimeEditField.Value
              break;
          end
      end

      app.CurrenttimeLabel.Text = sprintf('t = %.4g', app.t);
       app.SimulationProgressGauge.Value = 100 * min(1, app.t/max(app.TotalTimeEditField.Value, eps));
      function dt = getDT()
          if isstruct(app.model) && isfield(app.model,'dt'), dt = app.model.dt;
          else
              try dt = app.model.dt; catch, dt = app.TimeStepdtEditField.Value; end
          end
      end
      function v = getModelField(fname, defaultV)
          v = defaultV;
          if isstruct(app.model) && isfield(app.model,fname)
              v = app.model.(fname); return;
          end
          if ~isstruct(app.model)
              try v = app.model.(fname); catch, v = defaultV; end
          end
      end

        end
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
    
    % 1) Make Du editable 
    app.DuSpinner.Enable = 'on';
    app.DuSpinnerLabel.Enable = 'on';

    % 2) Safe defaults / fix dropdowns 
    app.ModelTypeDropDown.Items = {'Schnakenberg','Brusselator','Heat'};
    if ~ismember(app.ModelTypeDropDown.Value, app.ModelTypeDropDown.Items)
        app.ModelTypeDropDown.Value = 'Schnakenberg';
    end

    app.InitialConditionDropDown.Items = {'Uniform','Random','Localized Spot'};
    if ~ismember(app.InitialConditionDropDown.Value, app.InitialConditionDropDown.Items)
        app.InitialConditionDropDown.Value = 'Random';
    end

    app.VariabletodisplayDropDown.Items = {'u','v','Both (u+v)','Magnitude'};
    if ~ismember(app.VariabletodisplayDropDown.Value, app.VariabletodisplayDropDown.Items)
        app.VariabletodisplayDropDown.Value = 'u';
    end

    app.PlottypeDropDown.Items = {'2D Heatmap','3D Surface','Contour','2D Quiver'};
    if ~ismember(app.PlottypeDropDown.Value, app.PlottypeDropDown.Items)
        app.PlottypeDropDown.Value = '2D Heatmap';
    end

    app.ColormapDropDown.Items = {'parula','jet','hot','cool','turbo','viridis','gray'};
    if ~ismember(app.ColormapDropDown.Value, app.ColormapDropDown.Items)
        app.ColormapDropDown.Value = 'turbo';
    end

    % 3) Status + gauge + axes
    app.StatusLabel.Text = 'Ready.';
    app.CurrenttimeLabel.Text = 't = 0';
    app.SimulationProgressGauge.Value = 0;
    cla(app.UIAxes);
    drawnow;

    % 4) Optional: start with a fresh initial frame
    resetSimulation(app);
       
        end

        % Callback function
        function ParameterbParameterbSliderValueChanged(app, event)
    % optional: do nothing or update label etc.
    value = app.bSlider.Value;
        end

        % Value changed function: ModelTypeDropDown
        function ModelTypeDropDownValueChanged(app, event)
            app.model = app.ModelTypeDropDown.Value;
      switch app.model
          case 'Schnakenberg'
              app.aSlider.Enable = 'on';
              app.bSlider.Enable = 'on';
              app.aSlider.Limits = [0.0 4.0];
              app.bSlider.Limits = [0.0 4.0];
              app.aSlider.Value  = 0.1;
              app.bSlider.Value  = 0.9;
              app.aSliderLabel.Text = 'a';
              app.bSliderLabel.Text = 'b';

          case 'Brusselator'
              app.aSlider.Enable = 'on';
              app.bSlider.Enable = 'on';
              app.aSlider.Limits = [0.0 5.0];
              app.bSlider.Limits = [0.0 15.0];
              app.aSlider.Value  = 1.0;
              app.bSlider.Value  = 3.0;
              app.aSliderLabel.Text = 'a';
              app.bSliderLabel.Text = 'b';

          case 'Heat'
            
              app.aSlider.Enable = 'on';
              app.bSlider.Enable = 'on';
              app.aSlider.Limits = [0 1];
              app.bSlider.Limits = [0 1];
              app.aSlider.Value  = 0;
              app.bSlider.Value  = 0;
              app.aSlider.Enable = 'off';
              app.bSlider.Enable = 'off';
              app.aSliderLabel.Text = 'a';
              app.bSliderLabel.Text = 'b';
                  otherwise
      end

        end

        % Callback function
        function ResumeButtonPushed(app, event)
            app.StartPauseButton.Value = true;
            StartPauseButtonValueChanged(app, event);
        end

        % Callback function
        function PauseButtonPushed(app, event)
            if ~isempty(app.timerObj) && isvalid(app.timerObj)
                stop(app.timerObj);
            end
            app.StatusLabel.Text = 'Paused.';
            app.StartPauseButton.Value = false;
        end

        % Button pushed function: StopButton
        function StopButtonPushed(app, event)
             if ~isempty(app.timerObj) && isvalid(app.timerObj)
                stop(app.timerObj);
            end
            app.StatusLabel.Text = 'Stopped.';
            app.StartPauseButton.Value = false;
        end

        % Button pushed function: ResetButton
        function ResetButtonPushed(app, event)
            resetDashboardUI(app);   % UI reset + stop timer + clear axes etc.
            resetSimulation(app);    % backend rebuild + initial condition + first frame
            app.StatusLabel.Text = 'Reset done.';
        end

        % Button pushed function: StepButton
        function StepButtonPushed(app, event)
            if isempty(app.model) || isempty(app.geo)
                resetSimulation(app);
            end
            doOneBatch(app);
            renderFrame(app);
        end

        % Value changed function: StartPauseButton
        function StartPauseButtonValueChanged(app, event)
           if app.StartPauseButton.Value
          % Always rebuild clean
          resetSimulation(app);

          % Force one visible step immediately
          doOneBatch(app);
          renderFrame(app);
          drawnow;

          % Restart timer cleanly
          if ~isempty(app.timerObj) && isvalid(app.timerObj)
              stop(app.timerObj);
              delete(app.timerObj);
          end

          app.timerObj = timer( ...
              'ExecutionMode','fixedSpacing', ...
              'Period', max(0.01, 1/max(1, app.FrameRateSpinner.Value)), ...
              'TimerFcn', @(~,~)onTick(app));

          start(app.timerObj);
          app.StatusLabel.Text = 'Running...';
      else
          PauseButtonPushed(app, event);
      end


        end

        % Value changed function: PlottypeDropDown
        function PlottypeDropDownValueChanged(app, event)
          renderFrame(app);
            
        end

        % Value changed function: bSlider
        function bSliderValueChanged(app, event)
            value = app.bSlider.Value;
            
        end

        % Selection changed function: PatternsButtonGroup
        function PatternsButtonGroupSelectionChanged(app, event)

      app.ModelTypeDropDown.Value = 'Schnakenberg';
      app.ModelTypeDropDown.Enable = 'off';
      ModelTypeDropDownValueChanged(app, []);

      selected = event.NewValue; % or app.PatternsButtonGroup.SelectedObject

      if selected == app.TuringSpotsButton
          app.aSlider.Value = 0.1;
          app.bSlider.Value = 2.0;
          app.aSlider.Enable = 'off';
          app.bSlider.Enable = 'off';
          
          app.DuSpinner.Value = 1.0;
          app.DvSpinner.Value = 100.0;
          app.DuSpinner.Enable = 'off';
          app.DvSpinner.Enable = 'off';
          app.DuSpinnerLabel.Enable = 'off';
          app.DvSpinnerLabel.Enable = 'off';

      elseif selected == app.StripesButton
          app.aSlider.Value = 0.1;
          app.bSlider.Value = 2.0;
          app.aSlider.Enable = 'off';
          app.bSlider.Enable = 'off';

          app.DuSpinner.Value = 1.0;
          app.DvSpinner.Value = 30.0;
          app.DuSpinner.Enable = 'off';
          app.DvSpinner.Enable = 'off';
          app.DuSpinnerLabel.Enable = 'off';
          app.DvSpinnerLabel.Enable = 'off';

   

      else
          % Complex patterns = editable
          app.ModelTypeDropDown.Enable = 'on';
          app.aSlider.Enable = 'on';
          app.bSlider.Enable = 'on';

          app.DuSpinner.Enable = 'on';
          app.DvSpinner.Enable = 'on';
          app.DuSpinnerLabel.Enable = 'on';
          app.DvSpinnerLabel.Enable = 'on';
      end

       % Periodic boundary for paper cases
  app.NeumannButton.Value   = false;
  app.DirichletButton.Value = false;
  app.PeriodicButton.Value  = true;

  % Initial condition = steady state + noise
  a = app.aSlider.Value;
  b = app.bSlider.Value;
  app.uBaseValueEditField.Value = a + b;
  app.vBaseValueEditField.Value = b / (a + b)^2;

  app.InitialConditionDropDown.Value = 'Localized Spot';
  app.NoiseAmplitudeSlider.Value = 0.01;


      resetSimulation(app);
   
            
        end

        % Value changed function: HistoryRecordingCheckBox
        function HistoryRecordingCheckBoxValueChanged(app, event)
           if app.HistoryRecordingCheckBox.Value
                  app.StatusLabel.Text = 'History recording enabled.';
              else
                  app.StatusLabel.Text = 'History recording disabled — Save/Export blocked.';
              end
            
        end

        % Button pushed function: SaveDataButton
        function SaveDataButtonPushed(app, event)
             if ~app.HistoryRecordingCheckBox.Value
                  uialert(app.UIFigure, 'Enable "History recording" before saving or exporting.', 'History Recording');
                  return;
              end

              if isempty(app.u) || isempty(app.v) || isempty(app.geo)
                  uialert(app.UIFigure, 'No data to save. Run the simulation first.', 'Save Data');
                  return;
              end

              defaultName = sprintf('rd_data_t%.3f.mat', app.t);
              [file, path] = uiputfile( ...
                  {'*.mat','MAT-file (full data)'; '*.csv','CSV (current x,y,u,v)'}, ...
                  'Save Data', defaultName);

              if isequal(file,0)
                  return;
              end

              fullpath = fullfile(path, file);
              [~,~,ext] = fileparts(fullpath);

              params = getAllParameters(app);
              history = struct('t', app.tHistory, 'u', {app.uHistory}, 'v', {app.vHistory});

              x = app.geo.x(:);
              y = app.geo.y(:);
              u = app.u(:);
              v = app.v(:);
              t = app.t;

              tri = [];
              if isprop(app.geo,'t') && ~isempty(app.geo.t)
                  tri = app.geo.t;
              end

              if strcmpi(ext, '.csv')
                  T = table(x, y, u, v);
                  writetable(T, fullpath);

                  sidecar = fullfile(path, [file(1:end-4) '_params.mat']);
                  save(sidecar, 'params', 'history', 'tri');
              else
                  save(fullpath, 'x', 'y', 'u', 'v', 't', 'tri', 'params', 'history');
              end

              app.StatusLabel.Text = sprintf('Saved data: %s', fullpath);
        end

        % Button pushed function: ExportAnimationButton
        function ExportAnimationButtonPushed(app, event)
             if ~app.HistoryRecordingCheckBox.Value
                  uialert(app.UIFigure, 'Enable "History recording" before saving or exporting.', 'History Recording');
                  return;
              end

              if isempty(app.frameHistory)
                  uialert(app.UIFigure, 'No frames recorded. Run the simulation first.', 'Export Animation');
                  return;
              end

              [file, path] = uiputfile( ...
                  {'*.avi','AVI'; '*.mp4','MP4 (MPEG-4)'}, ...
                  'Export Animation', 'rd_animation.avi');

              if isequal(file,0)
                  return;
              end

              fullpath = fullfile(path, file);
              [~,~,ext] = fileparts(fullpath);

              if strcmpi(ext, '.mp4')
                  vw = VideoWriter(fullpath, 'MPEG-4');
              else
                  vw = VideoWriter(fullpath);
              end
              vw.FrameRate = max(1, round(app.FrameRateSpinner.Value));
              open(vw);

              for i = 1:numel(app.frameHistory)
                  writeVideo(vw, app.frameHistory(i));
              end

              close(vw);
              app.StatusLabel.Text = sprintf('Exported animation: %s', fullpath);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1184 721];
            app.UIFigure.Name = 'MATLAB App';

            % Create HelpMenu
            app.HelpMenu = uimenu(app.UIFigure);
            app.HelpMenu.Text = 'Help';

            % Create ManualMenu
            app.ManualMenu = uimenu(app.HelpMenu);
            app.ManualMenu.Text = 'Manual';

            % Create ModelandParametersPanel
            app.ModelandParametersPanel = uipanel(app.UIFigure);
            app.ModelandParametersPanel.Title = 'Model and Parameters';
            app.ModelandParametersPanel.Position = [1 477 440 245];

            % Create ModelTypeDropDownLabel
            app.ModelTypeDropDownLabel = uilabel(app.ModelandParametersPanel);
            app.ModelTypeDropDownLabel.HorizontalAlignment = 'right';
            app.ModelTypeDropDownLabel.Position = [14 192 67 22];
            app.ModelTypeDropDownLabel.Text = 'Model Type';

            % Create ModelTypeDropDown
            app.ModelTypeDropDown = uidropdown(app.ModelandParametersPanel);
            app.ModelTypeDropDown.Items = {'Schnakenberg', 'Brusselator', 'Heat'};
            app.ModelTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @ModelTypeDropDownValueChanged, true);
            app.ModelTypeDropDown.Position = [94 192 111 22];
            app.ModelTypeDropDown.Value = 'Schnakenberg';

            % Create DuSpinnerLabel
            app.DuSpinnerLabel = uilabel(app.ModelandParametersPanel);
            app.DuSpinnerLabel.HorizontalAlignment = 'right';
            app.DuSpinnerLabel.Enable = 'off';
            app.DuSpinnerLabel.Position = [5 66 25 22];
            app.DuSpinnerLabel.Text = 'Du';

            % Create bSlider
            app.bSlider = uislider(app.ModelandParametersPanel);
            app.bSlider.Limits = [0 4];
            app.bSlider.ValueChangedFcn = createCallbackFcn(app, @bSliderValueChanged, true);
            app.bSlider.Position = [52 123 146 3];
            app.bSlider.Value = 0.5;

            % Create DuSpinner
            app.DuSpinner = uispinner(app.ModelandParametersPanel);
            app.DuSpinner.Step = 0.1;
            app.DuSpinner.Limits = [0.001 100];
            app.DuSpinner.Enable = 'off';
            app.DuSpinner.Position = [46 67 120 22];
            app.DuSpinner.Value = 0.1;

            % Create DvSpinnerLabel
            app.DvSpinnerLabel = uilabel(app.ModelandParametersPanel);
            app.DvSpinnerLabel.HorizontalAlignment = 'right';
            app.DvSpinnerLabel.Position = [5 27 25 22];
            app.DvSpinnerLabel.Text = 'Dv';

            % Create DvSpinner
            app.DvSpinner = uispinner(app.ModelandParametersPanel);
            app.DvSpinner.Limits = [0.1 100];
            app.DvSpinner.Position = [45 27 121 22];
            app.DvSpinner.Value = 0.1;

            % Create DomainSizeEditFieldLabel
            app.DomainSizeEditFieldLabel = uilabel(app.ModelandParametersPanel);
            app.DomainSizeEditFieldLabel.HorizontalAlignment = 'right';
            app.DomainSizeEditFieldLabel.Position = [240 192 72 22];
            app.DomainSizeEditFieldLabel.Text = 'Domain Size';

            % Create DomainSizeEditField
            app.DomainSizeEditField = uieditfield(app.ModelandParametersPanel, 'numeric');
            app.DomainSizeEditField.Limits = [1 100];
            app.DomainSizeEditField.Position = [326 192 100 22];
            app.DomainSizeEditField.Value = 50;

            % Create GridPointsSpinnerLabel
            app.GridPointsSpinnerLabel = uilabel(app.ModelandParametersPanel);
            app.GridPointsSpinnerLabel.HorizontalAlignment = 'right';
            app.GridPointsSpinnerLabel.Position = [249 152 65 22];
            app.GridPointsSpinnerLabel.Text = 'Grid Points';

            % Create GridPointsSpinner
            app.GridPointsSpinner = uispinner(app.ModelandParametersPanel);
            app.GridPointsSpinner.Step = 10;
            app.GridPointsSpinner.Limits = [50 500];
            app.GridPointsSpinner.Position = [328 152 100 22];
            app.GridPointsSpinner.Value = 500;

            % Create TimeStepdtEditFieldLabel
            app.TimeStepdtEditFieldLabel = uilabel(app.ModelandParametersPanel);
            app.TimeStepdtEditFieldLabel.HorizontalAlignment = 'right';
            app.TimeStepdtEditFieldLabel.Position = [232 112 80 22];
            app.TimeStepdtEditFieldLabel.Text = 'Time Step (dt)';

            % Create TimeStepdtEditField
            app.TimeStepdtEditField = uieditfield(app.ModelandParametersPanel, 'numeric');
            app.TimeStepdtEditField.Position = [327 112 100 22];
            app.TimeStepdtEditField.Value = 0.1;

            % Create TotalTimeEditFieldLabel
            app.TotalTimeEditFieldLabel = uilabel(app.ModelandParametersPanel);
            app.TotalTimeEditFieldLabel.HorizontalAlignment = 'right';
            app.TotalTimeEditFieldLabel.Position = [252 66 60 22];
            app.TotalTimeEditFieldLabel.Text = 'Total Time';

            % Create TotalTimeEditField
            app.TotalTimeEditField = uieditfield(app.ModelandParametersPanel, 'numeric');
            app.TotalTimeEditField.Position = [327 66 100 22];
            app.TotalTimeEditField.Value = 200;

            % Create UpdateEveryNstepsSpinnerLabel
            app.UpdateEveryNstepsSpinnerLabel = uilabel(app.ModelandParametersPanel);
            app.UpdateEveryNstepsSpinnerLabel.HorizontalAlignment = 'right';
            app.UpdateEveryNstepsSpinnerLabel.Position = [190 27 122 22];
            app.UpdateEveryNstepsSpinnerLabel.Text = 'Update Every N steps';

            % Create UpdateEveryNstepsSpinner
            app.UpdateEveryNstepsSpinner = uispinner(app.ModelandParametersPanel);
            app.UpdateEveryNstepsSpinner.Limits = [1 100];
            app.UpdateEveryNstepsSpinner.Position = [327 27 100 22];
            app.UpdateEveryNstepsSpinner.Value = 10;

            % Create aSlider
            app.aSlider = uislider(app.ModelandParametersPanel);
            app.aSlider.Limits = [0 4];
            app.aSlider.Position = [54 165 144 3];
            app.aSlider.Value = 1;

            % Create aSliderLabel
            app.aSliderLabel = uilabel(app.ModelandParametersPanel);
            app.aSliderLabel.HorizontalAlignment = 'right';
            app.aSliderLabel.Position = [2 156 25 22];
            app.aSliderLabel.Text = 'a';

            % Create bSliderLabel
            app.bSliderLabel = uilabel(app.ModelandParametersPanel);
            app.bSliderLabel.HorizontalAlignment = 'right';
            app.bSliderLabel.Position = [3 114 25 22];
            app.bSliderLabel.Text = 'b';

            % Create VisualizationPanel
            app.VisualizationPanel = uipanel(app.UIFigure);
            app.VisualizationPanel.Title = 'Visualization';
            app.VisualizationPanel.Position = [441 242 744 480];

            % Create UIAxes
            app.UIAxes = uiaxes(app.VisualizationPanel);
            title(app.UIAxes, 'Plot')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [16 24 692 416];

            % Create ExportandStatusPanel
            app.ExportandStatusPanel = uipanel(app.UIFigure);
            app.ExportandStatusPanel.Title = 'Export and Status';
            app.ExportandStatusPanel.Position = [439 1 746 243];

            % Create SimulationProgressGaugeLabel
            app.SimulationProgressGaugeLabel = uilabel(app.ExportandStatusPanel);
            app.SimulationProgressGaugeLabel.HorizontalAlignment = 'center';
            app.SimulationProgressGaugeLabel.Position = [22 127 112 22];
            app.SimulationProgressGaugeLabel.Text = 'Simulation Progress';

            % Create SimulationProgressGauge
            app.SimulationProgressGauge = uigauge(app.ExportandStatusPanel, 'linear');
            app.SimulationProgressGauge.Position = [17 164 119 40];

            % Create StatusLabel
            app.StatusLabel = uilabel(app.ExportandStatusPanel);
            app.StatusLabel.Position = [23 96 39 22];
            app.StatusLabel.Text = 'Status';

            % Create CurrenttimeLabel
            app.CurrenttimeLabel = uilabel(app.ExportandStatusPanel);
            app.CurrenttimeLabel.Position = [24 61 71 22];
            app.CurrenttimeLabel.Text = 'Current time';

            % Create SaveDataButton
            app.SaveDataButton = uibutton(app.ExportandStatusPanel, 'push');
            app.SaveDataButton.ButtonPushedFcn = createCallbackFcn(app, @SaveDataButtonPushed, true);
            app.SaveDataButton.Position = [176 181 100 23];
            app.SaveDataButton.Text = 'Save Data';

            % Create ExportAnimationButton
            app.ExportAnimationButton = uibutton(app.ExportandStatusPanel, 'push');
            app.ExportAnimationButton.ButtonPushedFcn = createCallbackFcn(app, @ExportAnimationButtonPushed, true);
            app.ExportAnimationButton.Position = [177 139 100 23];
            app.ExportAnimationButton.Text = 'Export Animation';

            % Create FrameRateSpinnerLabel
            app.FrameRateSpinnerLabel = uilabel(app.ExportandStatusPanel);
            app.FrameRateSpinnerLabel.HorizontalAlignment = 'right';
            app.FrameRateSpinnerLabel.Position = [284 182 67 22];
            app.FrameRateSpinnerLabel.Text = 'Frame Rate';

            % Create FrameRateSpinner
            app.FrameRateSpinner = uispinner(app.ExportandStatusPanel);
            app.FrameRateSpinner.Limits = [1 60];
            app.FrameRateSpinner.Position = [366 182 100 22];
            app.FrameRateSpinner.Value = 30;

            % Create HistoryRecordingCheckBox
            app.HistoryRecordingCheckBox = uicheckbox(app.ExportandStatusPanel);
            app.HistoryRecordingCheckBox.ValueChangedFcn = createCallbackFcn(app, @HistoryRecordingCheckBoxValueChanged, true);
            app.HistoryRecordingCheckBox.Text = 'History Recording';
            app.HistoryRecordingCheckBox.Position = [322 148 118 22];

            % Create AdvancedFeaturesPanel
            app.AdvancedFeaturesPanel = uipanel(app.UIFigure);
            app.AdvancedFeaturesPanel.Title = 'Advanced Features';
            app.AdvancedFeaturesPanel.Position = [0 1 441 234];

            % Create VariabletodisplayDropDownLabel
            app.VariabletodisplayDropDownLabel = uilabel(app.AdvancedFeaturesPanel);
            app.VariabletodisplayDropDownLabel.HorizontalAlignment = 'right';
            app.VariabletodisplayDropDownLabel.Position = [14 180 103 22];
            app.VariabletodisplayDropDownLabel.Text = 'Variable to display';

            % Create VariabletodisplayDropDown
            app.VariabletodisplayDropDown = uidropdown(app.AdvancedFeaturesPanel);
            app.VariabletodisplayDropDown.Items = {'u', 'v', 'Both'};
            app.VariabletodisplayDropDown.Position = [132 180 100 22];
            app.VariabletodisplayDropDown.Value = 'u';

            % Create PlottypeDropDownLabel
            app.PlottypeDropDownLabel = uilabel(app.AdvancedFeaturesPanel);
            app.PlottypeDropDownLabel.HorizontalAlignment = 'right';
            app.PlottypeDropDownLabel.Position = [17 141 53 22];
            app.PlottypeDropDownLabel.Text = 'Plot type';

            % Create PlottypeDropDown
            app.PlottypeDropDown = uidropdown(app.AdvancedFeaturesPanel);
            app.PlottypeDropDown.Items = {'2D Heatmap', '3D Surface', 'Contour', '2D Quiver'};
            app.PlottypeDropDown.ValueChangedFcn = createCallbackFcn(app, @PlottypeDropDownValueChanged, true);
            app.PlottypeDropDown.Position = [131 141 100 22];
            app.PlottypeDropDown.Value = '2D Heatmap';

            % Create ColormapDropDownLabel
            app.ColormapDropDownLabel = uilabel(app.AdvancedFeaturesPanel);
            app.ColormapDropDownLabel.HorizontalAlignment = 'right';
            app.ColormapDropDownLabel.Position = [17 101 58 22];
            app.ColormapDropDownLabel.Text = 'Colormap';

            % Create ColormapDropDown
            app.ColormapDropDown = uidropdown(app.AdvancedFeaturesPanel);
            app.ColormapDropDown.Items = {'parula', 'jet', 'hot', 'cool', 'turbo', 'viridis'};
            app.ColormapDropDown.Position = [131 101 100 22];
            app.ColormapDropDown.Value = 'parula';

            % Create AutoscalecolorlimitsCheckBox
            app.AutoscalecolorlimitsCheckBox = uicheckbox(app.AdvancedFeaturesPanel);
            app.AutoscalecolorlimitsCheckBox.Text = 'Autoscale color limits';
            app.AutoscalecolorlimitsCheckBox.Position = [23 61 137 22];
            app.AutoscalecolorlimitsCheckBox.Value = true;

            % Create EnablerealtimeupdatesCheckBox
            app.EnablerealtimeupdatesCheckBox = uicheckbox(app.AdvancedFeaturesPanel);
            app.EnablerealtimeupdatesCheckBox.Text = 'Enable real time updates';
            app.EnablerealtimeupdatesCheckBox.Position = [24 24 155 22];
            app.EnablerealtimeupdatesCheckBox.Value = true;

            % Create PatternsButtonGroup
            app.PatternsButtonGroup = uibuttongroup(app.AdvancedFeaturesPanel);
            app.PatternsButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @PatternsButtonGroupSelectionChanged, true);
            app.PatternsButtonGroup.Title = 'Patterns';
            app.PatternsButtonGroup.Position = [247 24 178 183];

            % Create TuringSpotsButton
            app.TuringSpotsButton = uiradiobutton(app.PatternsButtonGroup);
            app.TuringSpotsButton.Text = 'Turing Spots';
            app.TuringSpotsButton.Position = [11 117 90 22];

            % Create StripesButton
            app.StripesButton = uiradiobutton(app.PatternsButtonGroup);
            app.StripesButton.Text = 'Stripes';
            app.StripesButton.Position = [11 95 60 22];

            % Create ComplexPatternsButton
            app.ComplexPatternsButton = uiradiobutton(app.PatternsButtonGroup);
            app.ComplexPatternsButton.Text = 'Complex Patterns';
            app.ComplexPatternsButton.Position = [11 73 119 22];
            app.ComplexPatternsButton.Value = true;

            % Create InitialConditionsandControlsPanel
            app.InitialConditionsandControlsPanel = uipanel(app.UIFigure);
            app.InitialConditionsandControlsPanel.Title = 'Initial Conditions and Controls';
            app.InitialConditionsandControlsPanel.Position = [0 234 441 248];

            % Create InitialConditionDropDownLabel
            app.InitialConditionDropDownLabel = uilabel(app.InitialConditionsandControlsPanel);
            app.InitialConditionDropDownLabel.HorizontalAlignment = 'right';
            app.InitialConditionDropDownLabel.Position = [15 186 88 22];
            app.InitialConditionDropDownLabel.Text = 'Initial Condition';

            % Create InitialConditionDropDown
            app.InitialConditionDropDown = uidropdown(app.InitialConditionsandControlsPanel);
            app.InitialConditionDropDown.Items = {'Uniform', 'Random', 'Localized Spot', 'Custom'};
            app.InitialConditionDropDown.Position = [118 186 90 22];
            app.InitialConditionDropDown.Value = 'Localized Spot';

            % Create uBaseValueEditFieldLabel
            app.uBaseValueEditFieldLabel = uilabel(app.InitialConditionsandControlsPanel);
            app.uBaseValueEditFieldLabel.HorizontalAlignment = 'right';
            app.uBaseValueEditFieldLabel.Position = [13 146 79 22];
            app.uBaseValueEditFieldLabel.Text = 'u₀ Base Value';

            % Create uBaseValueEditField
            app.uBaseValueEditField = uieditfield(app.InitialConditionsandControlsPanel, 'numeric');
            app.uBaseValueEditField.Position = [107 146 100 22];
            app.uBaseValueEditField.Value = 1;

            % Create vBaseValueEditFieldLabel
            app.vBaseValueEditFieldLabel = uilabel(app.InitialConditionsandControlsPanel);
            app.vBaseValueEditFieldLabel.HorizontalAlignment = 'right';
            app.vBaseValueEditFieldLabel.Position = [15 107 78 22];
            app.vBaseValueEditFieldLabel.Text = 'v₀ Base Value';

            % Create vBaseValueEditField
            app.vBaseValueEditField = uieditfield(app.InitialConditionsandControlsPanel, 'numeric');
            app.vBaseValueEditField.Position = [108 107 100 22];
            app.vBaseValueEditField.Value = 1;

            % Create BoundaryConditionsButtonGroup
            app.BoundaryConditionsButtonGroup = uibuttongroup(app.InitialConditionsandControlsPanel);
            app.BoundaryConditionsButtonGroup.Title = 'Boundary Conditions';
            app.BoundaryConditionsButtonGroup.Position = [216 107 100 101];

            % Create NeumannButton
            app.NeumannButton = uiradiobutton(app.BoundaryConditionsButtonGroup);
            app.NeumannButton.Text = 'Neumann';
            app.NeumannButton.Position = [11 55 74 22];

            % Create DirichletButton
            app.DirichletButton = uiradiobutton(app.BoundaryConditionsButtonGroup);
            app.DirichletButton.Text = 'Dirichlet';
            app.DirichletButton.Position = [11 33 66 22];
            app.DirichletButton.Value = true;

            % Create PeriodicButton
            app.PeriodicButton = uiradiobutton(app.BoundaryConditionsButtonGroup);
            app.PeriodicButton.Text = 'Periodic';
            app.PeriodicButton.Position = [11 11 66 22];

            % Create NoiseAmplitudeSlider
            app.NoiseAmplitudeSlider = uislider(app.InitialConditionsandControlsPanel);
            app.NoiseAmplitudeSlider.Limits = [0 0.05];
            app.NoiseAmplitudeSlider.Step = 0.1;
            app.NoiseAmplitudeSlider.Position = [128 65 150 3];

            % Create StartPauseButton
            app.StartPauseButton = uibutton(app.InitialConditionsandControlsPanel, 'state');
            app.StartPauseButton.ValueChangedFcn = createCallbackFcn(app, @StartPauseButtonValueChanged, true);
            app.StartPauseButton.Text = 'Start/Pause';
            app.StartPauseButton.Position = [329 185 100 23];

            % Create StepButton
            app.StepButton = uibutton(app.InitialConditionsandControlsPanel, 'push');
            app.StepButton.ButtonPushedFcn = createCallbackFcn(app, @StepButtonPushed, true);
            app.StepButton.Position = [327 137 100 23];
            app.StepButton.Text = 'Step';

            % Create ResetButton
            app.ResetButton = uibutton(app.InitialConditionsandControlsPanel, 'push');
            app.ResetButton.ButtonPushedFcn = createCallbackFcn(app, @ResetButtonPushed, true);
            app.ResetButton.Position = [328 89 100 23];
            app.ResetButton.Text = 'Reset';

            % Create StopButton
            app.StopButton = uibutton(app.InitialConditionsandControlsPanel, 'push');
            app.StopButton.ButtonPushedFcn = createCallbackFcn(app, @StopButtonPushed, true);
            app.StopButton.Position = [328 38 100 23];
            app.StopButton.Text = 'Stop';

            % Create NoiseAmplitudeSliderLabel
            app.NoiseAmplitudeSliderLabel = uilabel(app.InitialConditionsandControlsPanel);
            app.NoiseAmplitudeSliderLabel.HorizontalAlignment = 'right';
            app.NoiseAmplitudeSliderLabel.Position = [13 56 93 22];
            app.NoiseAmplitudeSliderLabel.Text = 'Noise Amplitude';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Simulator_GUI

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end