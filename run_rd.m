function varargout = run_rd(mode)
%RUN_RD  Select GUI or console mode for the reaction-diffusion model.
% Usage:
%   run_rd                 % prompt for mode
%   run_rd("gui")          % launch App Designer GUI
%   run_rd("console")      % run console mode (interactive inputs)
%
% Outputs:
%   If mode is "gui", returns the app handle.
%   If mode is "console", returns [u, v, geo, t] like run_rd_console.

    if nargin < 1 || strlength(string(mode)) == 0
        mode = promptMode();
    else
        mode = lower(string(mode));
    end

    switch mode
        case {"gui", "g"}
            app = Simulator_GUI;
            if nargout > 0
                varargout{1} = app;
            else
                clear app
            end

        case {"console", "c", "non-gui", "nongui", "no-gui"}
            if nargout > 0
                [varargout{1:nargout}] = run_rd_console();
            else
                run_rd_console();
            end

        otherwise
            error("Unknown mode '%s'. Use 'gui' or 'console'.", mode);
    end
end

function mode = promptMode()
    while true
        raw = input("Select mode (gui/console) [gui]: ", "s");
        if isempty(raw)
            mode = "gui";
            return;
        end
        raw = lower(strtrim(raw));
        if any(raw == ["gui", "g", "console", "c", "non-gui", "nongui", "no-gui"])
            mode = string(raw);
            return;
        end
        fprintf("Please enter 'gui' or 'console'.\n");
    end
end
