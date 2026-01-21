function app = Farhan_app_2_cli
%FARHAN_APP_2_CLI  Command-line launcher for the App Designer GUI.
% This wraps the .mlapp so it can be started like a normal .m function.

    app = Farhan_app_2;

    if nargout == 0
        clear app
    end
end
